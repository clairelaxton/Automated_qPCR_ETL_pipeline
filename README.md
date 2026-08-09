# qPCR Pipeline

Automates processing raw qPCR exports (Bio-Rad or ABI), merging them with a plate map and
target mapping, running QC, and collating everything into a clean summary — so you don't
have to do it by hand in Excel.

## Quick start
1) Click  to <> Code, Download ZIP
2) Unzip file and move contents to Project Folder
3) Install Python 3.10+ if you don't already have it: https://realpython.com/installing-python/
4) Set up a project folder for your data — this is separate from the pipeline code and can
live anywhere:

```text
your-project-folder/
├── Automated_qPCR_ETL_pipeline/ <- unzipped contents of this repository
├── plate_maps.csv          <- your exported plate maps (see below)
├── target_mapping.csv      <- optional, your target assignments (see below)
└── raw/                    <- folder with your qPCR .csv exports
    ├── Run1.csv
    ├── Run2.csv
    └── ...
```
### Option 1
1) Install dependencies (once, from this repo):

```bash
pip install -r requirements.txt
```
2) Run the pipeline against it:

```bash
python3 ./your/path/to/Automated_qPCR_ETL_pipeline/qpcr_pipeline_main.py /path/to/your-project-folder
```

(Or `cd` into your project folder first and just run `python3 /path/to/qpcr_pipeline_main.py`
with no argument — either order works, since either the argument or the current directory
tells it where your data is.)

### Option 2 
#### This is the recommended method if you're running into permissions issues with Python
Use the provided convenience script, which creates a virtual environment, installs
dependencies, and runs the pipeline in one step:
```bash
./your/path/to/run_pipeline.sh /path/to/your-project-folder
```

If you run into an issue saying that you don't have permissions to run this script, you can change your permissions by running:
```bash
chmod +x ./your/path/to/run_pipeline.sh
```
## plate_maps.csv format

One or more plates, stacked vertically in a single CSV, each shaped like this:

```csv
ECHH_MYLC1_112024,1,2,3,4,5,6,7,8,9,10,11,12
A,sample1,sample1,sample2,sample2,...,NTC,NTC
B,sample6,sample6,sample7,sample7,...,ECHH POS 50000,ECHH POS 50000
...
H,...

ECHH_MYLC2_112024,1,2,3,4,5,6,7,8,9,10,11,12
A,...
```

- **Run name** (top-left cell of each plate) **must match its raw CSV file's name exactly**
  (without `.csv`) — this is how the pipeline knows which raw export belongs to which plate.
  If it doesn't find a match, it warns you and, if the mismatch looks like a typo or a
  missing suffix, suggests which raw file it probably meant.
- **Sample names** should be unique per plate, one per well position (duplicate a sample's
  name across its replicate wells).
- **Separate plates with at least one blank line.**
- **Controls**: include `POS`/`POSITIVE`/`PC` (positive) or `NTC`/`NEC`/`NEG`/`NEGATIVE`
  (negative) somewhere in the sample name to have it recognized as a control — numbered
  variants like `NTC1`, `POS2` work too.

A few common formatting errors are bypassed automatically (you'll still see a warning
in the log, but the data is parsed correctly either way):
- A stray leading blank column shifting the whole plate one or more columns to the right.
- Row-letter cells (the `A`–`H` column) left blank — the row's position in the grid is used
  instead.

**Do not reuse the same run name for two different plate sections in the file.** Only the first
occurrence will be kept, and it's reported as an error.

## target_mapping.csv format (optional)

Maps each (run, fluorophore) pair to a target name:

```csv
RunPrefix,Fluorophore,Target
ECHH,FAM,EBV
ECHH,ROX,HHV6
RP-SC2,FAM,SC2
DEFAULT,ANY,Unknown
```

- **RunPrefix** matches the *start* of a run name (the longest matching prefix wins) — useful
  if you run different assays that reuse the same fluorophore for different targets. If you
  don't need that, just use `DEFAULT` for every row.
- **Fluorophore** must match the fluorophore/reporter name in your exported qPCR files
  (case-insensitive). `ANY` matches anything not otherwise mapped for that RunPrefix.
- **Target** is the name assigned wherever that pair appears. Each fluorophore can only map to
  one target per RunPrefix — if the same fluorophore needs different targets in different
  samples, distinguish those samples by name instead.
- If this file is missing entirely, every target defaults to `Unknown`.

The pipeline checks this file for common mistakes and warns you about them, usually
suggesting the likely correct spelling:
- A RunPrefix that doesn't match any of your actual run names (typo, or the wrong prefix).
- A (RunPrefix, Fluorophore) pair that shows up in your data but isn't covered here.
- Every target resolving to `Unknown` (a sign the file hasn't actually been filled in yet).
- Two rows mapping the same (RunPrefix, Fluorophore) to different targets (only the last one
  wins — this is almost always an accidental duplicate row).

## Internal control correction (optional)

If you run an internal control (e.g. RNase P) alongside your targets, the pipeline can
normalize `SQ Mean` against it — add `Correction Factor` and `Corrected SQ Mean` columns to
the clean summary. It's off unless you ask for it:

```bash
python3 qpcr_pipeline_main.py /path/to/your-project-folder --internal-control RNaseP
```

`RNaseP` here must match a `Target` value your `target_mapping.csv` actually produces
(case-insensitive) — it doesn't have to be RNase P specifically, it's whatever you've mapped
your internal control fluorophore to.

For each sample, the matching internal control reading is looked up in tiers: same sample +
same "batch" first, then the same sample in any batch, and if neither exists, that row is
left uncorrected (`Correction Factor = 1.0`) and flagged `missing_internal_control`. A
"batch" is derived from the run name — by default the middle underscore-delimited token
(`ECHH_MYLC1_112024` → batch `MYLC1`); if your run names don't follow that convention, a
developer can set `internal_control_batch_regex` in `PipelineConfig` to a custom pattern with
one capturing group.

If `TARGET_NAME` doesn't match anything in your data, the pipeline stops with an error
listing the target names it actually found — check `target_mapping.csv` and the flag itself
for spelling or case.

## Outputs

Files are named after your project folder.

- In the project folder:
  - `{folder}_qpcr_combined_results.csv` — every qPCR result, merged with the plate and target maps.
- In `results/`:
  - `{folder}_qpcr_clean_summary.xlsx` — one row per Sample/Target, QC'd and collapsed across
    runs. Includes `Correction Factor`/`Corrected SQ Mean` if internal control correction was enabled.
  - `{folder}_qpcr_exclusions.xlsx` — everything excluded during QC (inconclusive calls, ties,
    weak single positives, ...), which samples need to be repeated, and the controls that were
    set aside before collapsing.
- In `logs/`:
  - `qpcr_pipeline_<timestamp>.log` — a full, timestamped record of the run, including every
    warning described above.

## Notes

- QC thresholds and other defaults live in `src/config.py` (`PipelineConfig`).
- If a run is in `plate_maps.csv` but has no matching file in `raw/`, that plate's samples
  are skipped entirely and a warning is logged (see plate_maps.csv format, above).
