# Automated_qPCR_ETL_pipeline README
This pipeline automates processing raw qPCR exports (Bio-Rad or ABI), merging them with a plate map and target mapping, running QC, and collating everything into a clean summary — so you don't have to do it by hand in Excel.

## Citation:
Please cite the Preprint: https://www.medrxiv.org/content/10.64898/2026.05.19.26353495v1.
A detailed description of how the pipeline processes your qPCR data can be found in the Supplementary Methods section of this paper too, and the PDF
is included in this Repository.

### Important notes:
- This pipeline does not calculate quantity from standard curves or relative gene expression e.g. by ΔΔCq. It simply maps Samples to
  plate maps and fluorphore targets, runs basic QC to filter out technically erroneous results, and aggregates technical replicates.
- QC thresholds and other defaults live in `src/config.py` (`PipelineConfig`).

# Quick Start
## Download and install dependancies
1) Click <> Code, Download ZIP.
2) Install Python 3.10+ if you don't already have it: https://realpython.com/installing-python/
3) Set up a project folder for your data, this can live anywhere, but you'll need to use the path, or cd inside it to run the code.

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
If you're having issues installing requirements, you can install each package individually:
```bash
pip3 install pandas
pip3 install numpy
pip3 install openpyxl
```
If you're still having issues, try Option 2.

#### Run the pipeline
2) Run the pipeline against it:

```bash
python3 ./your/path/to/Automated_qPCR_ETL_pipeline/qpcr_pipeline_main.py /path/to/your-project-folder
```

### Option 2 
#### This is recommended if you're running into permissions issues with Python
Use the provided convenience script, which creates a virtual environment, installs
dependencies, and runs the pipeline in one step:
```bash
./your/path/to/run_pipeline.sh /path/to/your-project-folder
```

If you run into an issue saying that you don't have permissions to run this script, you can change your permissions by running:
```bash
chmod +x ./your/path/to/run_pipeline.sh
```
Note: you cannot use modifiers (e.g. internal control toggled on using this method)

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
A template is included in the repository that can be used for 96-well plates.
### Some important notes:
- **Run name** (top-left cell of each plate) **must match its raw CSV file's name exactly**
  (without `.csv`) — this is how the pipeline knows which raw export belongs to which plate.
  If it doesn't find a match, it warns you and, if the mismatch looks like a typo or a
  missing suffix, suggests which raw file it probably meant.
      _We recommend the naming convention: RunPrefix_ProjectName+RunNumber_Date (e.g. ECHH_MYLC1_112024)_
- **Sample names** must be unique per sample, duplicates will be treated as replicates and aggregated by the pipeline.
- **Separate plates with at least one empty line.**
- **Controls**: must include `POS`/`POSITIVE`/`PC` (positive) or `NTC`/`NEC`/`NEG`/`NEGATIVE`
  (negative) somewhere in the sample name to have it recognized as a control.

A few common formatting errors are bypassed automatically (you'll still see a warning
in the log, but the data is parsed correctly either way):
- A stray leading blank column shifting the whole plate one or more columns to the right.
- Row-letter cells (the `A`–`H` column) left blank — the row's position in the grid is used
  instead.

**Do not reuse the same run name for two different plate sections in the file.** Only the first
occurrence will be kept, and it's reported as an error.

## target_mapping.csv format

Maps each (run, fluorophore) pair to a target name, necassary for multiplex experiments:

```csv
RunPrefix,Fluorophore,Target
ECHH,FAM,EBV
ECHH,ROX,HHV6
RP-SC2,FAM,SC2
DEFAULT,ANY,Unknown
```
There is a template in the repository that demonstrates how this can be used for multiple multiplex assays.
### Some important notes:
- **RunPrefix** matches the *start* of a run name (E.g. ECHH_MYLC1_112024, ECHH is the run prefix, denoting that specific multiplex assay).
  This allows you to run the code on multiple different multiplex qPCR plates from different assays that reuse the same fluorophore for
  different targets. If you don't need that, use `DEFAULT`.
- **Fluorophore** must match the fluorophore/reporter name in your exported qPCR files (case-insensitive).
- **Target** is the genetic target assigned to a particular fluorophore (or RunPrefix/fluorophore combination). Each fluorophore
  can only map to *one target per RunPrefix*.
- If the target_mapping.csv file is missing entirely, every target name defaults to `Unknown`, this would be fine if you were only using
  a single fluorophore for all your runs (e.g. SYBR).

The pipeline is not configured to run multiple different assays that use the same fluorophores *within the same plate*, as the RunPrefix (which
denotes the entire plate) is what is used to distinguish this. If this is your goal, you'll need to denote your targets in the Sample name itself in
your plate map.

The pipeline checks the target_mapping.csv file for common mistakes and warns you about them, usually suggesting the likely correct spelling:
- A RunPrefix that doesn't match any of your actual run names (typo, or the wrong prefix).
- A (RunPrefix, Fluorophore) pair that shows up in your data but isn't covered here.
- Every target resolving to `Unknown` (a sign the file hasn't actually been filled in yet).
- Two rows mapping the same (RunPrefix, Fluorophore) to different targets (only the last one
  wins — this is almost always an accidental duplicate row).

# Modifiers
Modifiers can be called for, but you will need to use Option 1 to install the dependancies and run the code, _not Option 2_ (the bash script).

## Internal control correction (optional)
If you run an internal control (e.g. RNaseP) alongside your targets, the pipeline can
normalize `SQ Mean` (starting quantity) against it — add `Correction Factor` and `Corrected SQ Mean` columns to
the clean summary. It's switched off unless you call for it using the ```--internal control``` modifier:

```bash
python3 qpcr_pipeline_main.py /path/to/your-project-folder --internal-control RNaseP
```

`RNaseP` (or other internal control, IC) here must match a `Target` value your `target_mapping.csv` actually produces (case-insensitive).

For each sample, the matching internal control reading is looked up in tiers: same sample + same batch (by Run names e.g. ECHH_MYLC1_112024 and RP
SC2_MYLC1_112024 first (assuming you ran both plates on the same day), then the same sample in any plate with a matching RunPrefix and Sample (in
case you didn't). If neither exists, that row is left uncorrected (`Correction Factor = 1.0`) and flagged `missing_internal_control`.

If the internal control `TARGET_NAME` doesn't match anything in your data, the pipeline stops with an error
listing the target names it actually found — check `target_mapping.csv` and the flag itself for spelling or case.

# Outputs

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

- If a run is in `plate_maps.csv` but has no matching file in `raw/`, that plate's samples
  are skipped entirely and a warning is logged (see plate_maps.csv format, above).
