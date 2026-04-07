# GenBank viral downloader and extractor

## What this script does

1. Downloads all GenBank files matching `gbvrl*.seq.gz` from the NCBI FTP server.
2. Parses every GenBank record, **skipping SARS-CoV-2 records entirely** (they are not written to any output).
3. Writes a CSV of all remaining records with:
   - `accession`
   - `taxonomy` (from `ORGANISM` line)
   - `submission_date` (from `LOCUS` line)
   - `source_file`
4. Optionally extracts records for any virus by taxonomy or common name via `--query`.

## Setup

```bash
python3 -m venv env
./env/bin/pip install requests pyyaml tqdm
cp config.yml config.yml  # edit as needed
```

Dependencies: `requests`, `pyyaml`, `tqdm`.

## Configuration (`config.yml`)

| Key | Default | Description |
|-----|---------|-------------|
| `db_path` | `.` | Base directory for all downloads and outputs |
| `date` | (required) | Release label used as a subfolder (e.g. `2026-02-15`) |
| `download_folder` | `downloads` | Subfolder for raw `.seq.gz` files |
| `processed_folder` | `processed` | Subfolder for CSV outputs |
| `gb_file_prefix` | `gbvrl` | Prefix for GenBank file names |
| `sars2_taxonomy_terms` | (list) | Taxonomy strings used to identify and skip SARS-CoV-2 records |
| `write_selected_records` | `true` | Write `.seq.gz` record files for SARS2 and extracted selections |
| `extract.output_folder` | `extracted` | Subfolder for custom extraction results |
| `common_name_table` | (built-in) | Map of common names to taxonomy strings |

### Query resolution

`--query` terms are matched against `common_name_table` keys (case-insensitive). If found, the mapped taxonomy strings are used. If not found, the term is used as-is for direct taxonomy matching.

### Pre-configured common names

The included `config.yml` ships with a `common_name_table` covering many human-relevant viruses, including:

`HIV-1`, `HIV-2`, `SIV`, `SARS2`, `SARS`, `MERS`, `Influenza A/B/C`, `RSV`, `Hepatitis A–G`, `HPV`, `Dengue`, `Zika`, `Ebola`, `Mpox`, and many more.

Use these short names directly with `--query` (e.g. `--query HIV-1`) without needing to know the full taxonomy string.

## Output layout

```
<db_path>/
  <download_folder>/<date>/
    gbvrl1.seq.gz          # deleted after processing
    ...
  <processed_folder>/<date>/
    all_records.csv        # all records (SARS-CoV-2 excluded)
    metadata.csv           # extracted records (all non-SARS2 if no --query; query-filtered if --query)
    records/
      gbvrl1.seq.gz        # raw sequences for extracted records
      ...
```

## Usage

```
python download_seq.py [--config CONFIG] [--query TERM ...] [--list-files]
                       [--only-file FILE | --only-index N]
```

### Options

| Flag | Description |
|------|-------------|
| `--config CONFIG` | Path to YAML config file (default: `config.yml`) |
| `--query TERM` | Common name or taxonomy query. Repeat for multiple terms. |
| `--list-files` | List all `gbvrl*.seq.gz` files sorted by size (largest first) then exit |
| `--only-file FILE` | Download and process only this file (e.g. `gbvrl1.seq.gz`) |
| `--only-index N` | Download and process only the Nth file from the size-sorted list (1-based) |

`--only-file` and `--only-index` are mutually exclusive.

### Examples

```bash
# Download all files and extract HIV-1 records (common name lookup)
./env/bin/python download_seq.py --query HIV-1

# Also works with full taxonomy strings directly
./env/bin/python download_seq.py --query "Human immunodeficiency virus 1"

# Extract multiple viruses at once
./env/bin/python download_seq.py --query HIV-1 --query SARS2

# List available files with sizes
./env/bin/python download_seq.py --list-files

# Process only the largest file and extract HIV-1
./env/bin/python download_seq.py --only-index 1 --query HIV-1

# Process a specific file
./env/bin/python download_seq.py --only-file gbvrl1.seq.gz --query HIV-1

# Use an alternate config file
./env/bin/python download_seq.py --config my_config.yml --query Dengue
```

### Notes

- Already-downloaded `.seq.gz` files are skipped (re-runs are safe).
- If `--query` is omitted, custom extraction is disabled; SARS-CoV-2 splitting still runs.
- `--list-files` without `--only-file`/`--only-index` prints the catalog and exits without downloading.
