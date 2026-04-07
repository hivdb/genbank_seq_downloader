import contextlib
import csv
import gzip
import logging
import os
import re
import shutil
import argparse
from multiprocessing import Pool, Process, Queue
from pathlib import Path

import requests
import yaml
from tqdm import tqdm

GENBANK_FTP = "https://ftp.ncbi.nih.gov/genbank/"
DATE_RE = re.compile(r"(\d{2}-[A-Z]{3}-\d{4})$")
HEADER = ["accession", "taxonomy", "submission_date", "source_file"]


def load_yaml(file_path):
    return yaml.safe_load(open(file_path))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Download gbvrl files and extract GenBank metadata."
    )
    parser.add_argument(
        "--config",
        default="config.yml",
        help="Path to YAML config file (default: config.yml).",
    )
    parser.add_argument(
        "--query",
        action="append",
        default=[],
        help=(
            "Extraction query term from CLI. Use multiple --query flags for "
            "multiple values. Example: --query HIV-1 --query SARS2"
        ),
    )
    parser.add_argument(
        "--list-files",
        action="store_true",
        help="List all gbvrl files sorted by size (desc).",
    )
    select_group = parser.add_mutually_exclusive_group()
    select_group.add_argument(
        "--only-file",
        default=None,
        help="Download and process only this one file (e.g. gbvrl1.seq.gz).",
    )
    select_group.add_argument(
        "--only-index",
        type=int,
        default=None,
        help="Download and process only one file by rank from size-sorted list (1-based).",
    )
    return parser.parse_args()


def get_file_name_list(prefix):
    page = requests.get(GENBANK_FTP)
    page.raise_for_status()
    file_names = sorted(
        set(re.findall(rf"({re.escape(prefix)}\d+\.seq\.gz)", page.text)),
        key=lambda name: int(re.search(r"(\d+)", name).group(1)),
    )
    return file_names


def get_remote_file_size(file_name):
    response = requests.head(GENBANK_FTP + file_name, allow_redirects=True)
    response.raise_for_status()
    try:
        return int(response.headers.get("Content-Length", ""))
    except (TypeError, ValueError):
        return -1


def get_file_catalog(prefix, include_size=False):
    file_names = get_file_name_list(prefix)
    if not include_size:
        return [{"name": name, "size": -1} for name in file_names]
    return [{"name": name, "size": get_remote_file_size(name)} for name in file_names]


def sort_catalog_by_size_desc(catalog):
    return sorted(catalog, key=lambda x: (x["size"], x["name"]), reverse=True)


def print_catalog(catalog):
    print("rank\tfile_name\tsize_bytes\tsize_mb")
    for idx, item in enumerate(catalog, start=1):
        size = item["size"]
        size_mb = f"{(size / (1024 * 1024)):.2f}" if size >= 0 else "NA"
        print(f"{idx}\t{item['name']}\t{size}\t{size_mb}")
    print(f"Total: {len(catalog)} files")


def download_file(url, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        return
    part = Path(str(output_path) + ".part")
    if part.exists():
        part.unlink()
    try:
        with requests.get(url, stream=True) as response:
            response.raise_for_status()
            with open(part, "wb") as fd:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        fd.write(chunk)
        part.rename(output_path)
    except Exception:
        if part.exists():
            part.unlink()
        raise



def normalize_text(value):
    return (value or "").strip().lower()


def normalize_to_list(value):
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value]
    return [str(value)]



def taxonomy_matches(taxonomy, terms):
    tax = normalize_text(taxonomy)
    if not tax:
        return False
    for term in terms:
        t = normalize_text(term)
        if t and (tax == t or t in tax or tax in t):
            return True
    return False


def parse_record(record_lines):
    accession = ""
    taxonomy = ""
    submission_date = ""

    for line in record_lines:
        if line.startswith("ACCESSION") and not accession:
            parts = line.split()
            if len(parts) > 1:
                accession = parts[1]
        if line.startswith("  ORGANISM") and not taxonomy:
            taxonomy = line.replace("  ORGANISM", "", 1).strip()
        if line.startswith("LOCUS") and not submission_date:
            m = DATE_RE.search(line.strip())
            if m:
                submission_date = m.group(1)

    return accession, taxonomy, submission_date


def resolve_taxonomy_queries(config):
    raw_queries = normalize_to_list(config.get("query", []))
    if not raw_queries:
        return []

    common_lookup = {
        normalize_text(k): normalize_to_list(v)
        for k, v in config.get("common_name_table", {}).items()
    }

    taxonomies = []
    for query in raw_queries:
        q = normalize_text(query)
        if not q:
            continue
        mapped = common_lookup.get(q)
        if mapped:
            taxonomies.extend(mapped)
        else:
            taxonomies.append(query)

    # Preserve order, remove duplicates.
    seen = set()
    deduped = []
    for item in taxonomies:
        key = normalize_text(item)
        if key and key not in seen:
            seen.add(key)
            deduped.append(item)

    return deduped


def _csv_row_count(path):
    with open(path, newline="", encoding="utf-8") as f:
        return max(0, sum(1 for _ in f) - 1)  # minus header


def process_one_file(file_path, sars2_terms, extract_terms,
                     all_csv, extract_csv, write_records, extract_records_dir):
    stem = file_path.name.replace(".seq.gz", "")

    # If already complete, return cached counts.
    if all_csv.exists():
        total = _csv_row_count(all_csv)
        extracted = _csv_row_count(extract_csv) if (extract_csv and extract_csv.exists()) else 0
        return total, 0, extracted

    # Clean up any leftover .part files from a previous interrupted run.
    all_csv_part = Path(str(all_csv) + ".part")
    ext_csv_part = Path(str(extract_csv) + ".part") if extract_csv else None
    rec_part = (extract_records_dir / f"{stem}.seq.gz.part") if (write_records and extract_csv) else None
    for p in filter(None, [all_csv_part, ext_csv_part, rec_part]):
        if p.exists():
            p.unlink()

    all_csv_part.parent.mkdir(parents=True, exist_ok=True)
    if ext_csv_part:
        ext_csv_part.parent.mkdir(parents=True, exist_ok=True)
    if rec_part:
        extract_records_dir.mkdir(parents=True, exist_ok=True)

    total_records = 0
    skipped_sars2 = 0
    extracted_records = 0

    rec_ctx = gzip.open(rec_part, "wt", encoding="utf-8") if rec_part else contextlib.nullcontext()
    ext_ctx = open(ext_csv_part, "w", newline="", encoding="utf-8") if ext_csv_part else contextlib.nullcontext()

    with open(all_csv_part, "w", newline="", encoding="utf-8") as all_fd, \
         ext_ctx as extract_csv_fd, \
         rec_ctx as extract_records_fd:
        all_writer = csv.writer(all_fd)
        all_writer.writerow(HEADER)
        extract_writer = None
        if extract_csv_fd:
            extract_writer = csv.writer(extract_csv_fd)
            extract_writer.writerow(HEADER)

        with gzip.open(file_path, "rt", encoding="utf-8", errors="replace") as fd:
            record_lines = []
            for line in fd:
                record_lines.append(line)
                if line.startswith("//"):
                    total_records += 1
                    accession, taxonomy, submission_date = parse_record(record_lines)

                    if taxonomy_matches(taxonomy, sars2_terms):
                        skipped_sars2 += 1
                        record_lines = []
                        continue

                    row = [accession, taxonomy, submission_date, file_path.name]
                    all_writer.writerow(row)

                    if extract_writer and (not extract_terms or taxonomy_matches(taxonomy, extract_terms)):
                        extracted_records += 1
                        extract_writer.writerow(row)
                        if extract_records_fd:
                            extract_records_fd.write("".join(record_lines))

                    record_lines = []

    # Rename .part → final only on full success.
    all_csv_part.rename(all_csv)
    if ext_csv_part:
        ext_csv_part.rename(extract_csv)
    if rec_part:
        rec_part.rename(extract_records_dir / f"{stem}.seq.gz")

    return total_records, skipped_sars2, extracted_records


def download_worker(args):
    file_name, download_dir = args
    download_file(GENBANK_FTP + file_name, download_dir / file_name)
    return file_name


def producer_worker(selected_files, download_dir, tmp_dir, file_queue):
    for file_name in selected_files:
        stem = file_name.replace(".seq.gz", "")
        all_csv = tmp_dir / "all" / f"{stem}.csv"
        if not all_csv.exists():
            download_worker((file_name, download_dir))
        file_queue.put(file_name)  # blocks when queue is full
    file_queue.put(None)  # sentinel


def process_worker(args):
    file_name, download_dir, tmp_dir, sars2_terms, extract_terms, extract_enabled, write_records, extract_records_dir = args
    file_path = download_dir / file_name
    stem = file_name.replace(".seq.gz", "")
    all_csv = tmp_dir / "all" / f"{stem}.csv"
    extract_csv = tmp_dir / "extract" / f"{stem}.csv" if extract_enabled else None
    t, s, e = process_one_file(file_path, sars2_terms, extract_terms, all_csv, extract_csv, write_records, extract_records_dir)

    if file_path.exists():
        file_path.unlink()
    return file_name, t, s, e


def merge_csv_files(input_files, output_file):
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w", newline="", encoding="utf-8") as out_fd:
        writer = csv.writer(out_fd)
        header_written = False
        for csv_file in sorted(input_files):
            if not csv_file.exists():
                continue
            with open(csv_file, "r", newline="", encoding="utf-8") as in_fd:
                for i, row in enumerate(csv.reader(in_fd)):
                    if i == 0:
                        if not header_written:
                            writer.writerow(row)
                            header_written = True
                    else:
                        writer.writerow(row)
        if not header_written:
            writer.writerow(HEADER)


def extract_and_split(ctx):
    download_dir = ctx["download_dir"]
    output_dir = ctx["output_dir"]
    tmp_dir = output_dir / "_tmp"

    sars2_terms = normalize_to_list(ctx.get("sars2_taxonomy_terms", []))
    extract_terms = resolve_taxonomy_queries(ctx)
    extract_enabled = bool(ctx.get("extract", True))
    write_records = bool(ctx.get("write_selected_records", True))

    extract_records_dir = output_dir / "records"
    extract_path = output_dir / "metadata.csv"

    selected_files = ctx["selected_files"]
    n = len(selected_files)

    tmp_dir.mkdir(parents=True, exist_ok=True)

    cpu_count = os.cpu_count() or 1
    num_processes = ctx.get("num_processes") or cpu_count

    log_path = output_dir / "download_seq.log"
    log = logging.getLogger("download_seq")
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setFormatter(logging.Formatter("%(asctime)s  %(message)s", datefmt="%H:%M:%S"))
    log.addHandler(fh)

    print(f"Files: {n}  CPUs: {cpu_count}  Processes: {num_processes}")
    if extract_enabled:
        print(f"Extract metadata: {extract_path}")

    file_queue = Queue(maxsize=num_processes)
    prod = Process(target=producer_worker, args=(selected_files, download_dir, tmp_dir, file_queue), daemon=True)
    prod.start()

    results = []
    dl_bar  = tqdm(total=n, desc="Download", position=0, unit="file")
    proc_bar = tqdm(total=n, desc="Process ", position=1, unit="file")

    def flush_batch(pool, batch):
        pargs = [
            (f, download_dir, tmp_dir, sars2_terms, extract_terms,
             extract_enabled, write_records, extract_records_dir)
            for f in batch
        ]
        for file_name, t, s, e in pool.imap_unordered(process_worker, pargs):
            log.info(f"done: {file_name}  total={t}  sars2={s}  extracted={e}")
            proc_bar.update(1)
            results.append((file_name, t, s, e))

    with Pool(num_processes) as pool:
        batch = []
        while True:
            file_name = file_queue.get()
            if file_name is None:
                if batch:
                    flush_batch(pool, batch)
                break
            log.info(f"downloaded: {file_name}")
            dl_bar.update(1)
            batch.append(file_name)
            if len(batch) >= num_processes:
                flush_batch(pool, batch)
                batch = []

    dl_bar.close()
    proc_bar.close()
    prod.join()

    merge_csv_files((tmp_dir / "all").glob("*.csv"), output_dir / "all_records.csv")
    if extract_enabled:
        merge_csv_files((tmp_dir / "extract").glob("*.csv"), extract_path)

    shutil.rmtree(tmp_dir, ignore_errors=True)

    totals = [sum(r[i + 1] for r in results) for i in range(3)]
    print(f"Total processed: {totals[0]}  sars2 skipped: {totals[1]}  extracted: {totals[2]}")
    print(f"Total processed records: {totals[0]}")
    print(f"Total SARS-CoV-2 skipped: {totals[1]}")
    if extract_enabled:
        print(f"Total extracted records: {totals[2]}")


def work():
    args = parse_args()
    config = load_yaml(args.config)

    db_path = Path(config["db_path"]).expanduser().resolve()
    date_folder = str(config.get("date", "latest"))

    config["download_dir"] = db_path / config.get("download_folder", "downloads") / date_folder
    config["output_dir"] = db_path / config.get("processed_folder", "processed") / date_folder
    config["gb_file_prefix"] = config.get("gb_file_prefix", "gbvrl")
    config["write_selected_records"] = bool(config.get("write_selected_records", True))

    config["download_dir"].mkdir(parents=True, exist_ok=True)
    config["output_dir"].mkdir(parents=True, exist_ok=True)

    need_size_catalog = bool(args.list_files or args.only_index is not None)
    catalog = get_file_catalog(config["gb_file_prefix"], include_size=need_size_catalog)

    if args.only_file:
        names = {item["name"] for item in catalog}
        if args.only_file not in names:
            raise ValueError(f"--only-file not found: {args.only_file}")
        selected_files = [args.only_file]
    elif args.only_index is not None:
        sorted_catalog = sort_catalog_by_size_desc(catalog)
        if args.only_index < 1 or args.only_index > len(sorted_catalog):
            raise ValueError(
                f"--only-index out of range: {args.only_index}; valid 1..{len(sorted_catalog)}"
            )
        selected_files = [sorted_catalog[args.only_index - 1]["name"]]
    else:
        selected_files = [item["name"] for item in catalog]

    if args.list_files:
        sorted_catalog = sort_catalog_by_size_desc(
            catalog if need_size_catalog else get_file_catalog(config["gb_file_prefix"], include_size=True)
        )
        print_catalog(sorted_catalog)
        if args.only_file is None and args.only_index is None:
            print("No file selected for processing. Use --only-file or --only-index.")
            return

    config["selected_files"] = selected_files

    config["query"] = args.query

    extract_and_split(config)


if __name__ == "__main__":
    work()
