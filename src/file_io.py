"""
File I/O for the qPCR pipeline: locating and parsing raw instrument CSVs
(Bio-Rad and ABI exports), and transforming the plate map CSV into a
long-format Run/Well/Sample table.
"""
import pandas as pd
import re
from pathlib import Path
import logging

from .config import (
    ABI_INDICATORS, COLUMN_ALIASES, HEADER_DETECTION_GROUPS,
    COL, EMPTY_WELL_SAMPLE, PipelineConfig,
)


logger = logging.getLogger(__name__)


def split_header_line(line: str) -> list[str]:
    """
    Split a header candidate line on comma, tab, or semicolon.

    Args:
        line: Header line to split

    Returns:
        List of stripped, non-empty tokens (column names)
    """
    parts = re.split(r'[,\t;]', line.strip())
    return [p.strip() for p in parts if p.strip()]


def normalize_name(name: str) -> str:
    """
    Normalize a column name for case/whitespace-insensitive matching.

    Args:
        name: Column name to normalize

    Returns:
        Lowercase name with internal whitespace collapsed to single spaces
    """
    if not isinstance(name, str):
        name = str(name)
    name = name.strip()
    name = re.sub(r'\s+', ' ', name)
    return name.lower()


def normalize_well(well: str) -> str:
    """
    Normalize well IDs to the canonical 'A01' format.

    Args:
        well: Well ID in various formats (e.g., 'a1', 'A1', 'A01')

    Returns:
        Normalized well ID (e.g., 'A01'), '' for empty/NaN input, or the
        original value unchanged if it doesn't match a row-letter+column pattern
    """
    if pd.isna(well) or well == '':
        return ''

    well = str(well).strip().upper()
    match = re.match(r'^([A-H])0?(\d{1,2})$', well)

    if match:
        row = match.group(1)
        col = int(match.group(2))
        return f"{row}{col:02d}"

    return well


# Normalized (lowercase) alias sets, one per required column group, used for
# header detection. Built once at import time from config.HEADER_DETECTION_GROUPS
# so the alias spellings only have to be maintained in one place.
_REQUIRED_ALIAS_SETS = {
    group: {normalize_name(alias) for alias in aliases}
    for group, aliases in HEADER_DETECTION_GROUPS.items()
}


def find_qpcr_header_row(
    file_path: Path,
    max_scan_lines: int = 200,
    config: PipelineConfig | None = None,
) -> tuple[int, str]:
    """
    Scan a qPCR CSV for its header row and identify the instrument format.

    Args:
        file_path: Path to qPCR CSV file
        max_scan_lines: Maximum number of lines to scan
        config: Pipeline configuration (optional; overrides max_scan_lines)

    Returns:
        Tuple of (header_row_index, profile_str) where profile is 'biorad' or 'abi'

    Raises:
        ValueError: If no line within max_scan_lines has all required columns
    """
    if config is not None:
        max_scan_lines = config.max_scan_lines

    with file_path.open('r', encoding='utf-8', errors='ignore') as f:
        for idx, line in enumerate(f):
            if idx > max_scan_lines:
                break

            cols = split_header_line(line)
            if not cols:
                continue

            norm_set = {normalize_name(c) for c in cols}

            # A header row must contain at least one alias for every
            # required column.
            has_all_required = all(
                aliases & norm_set for aliases in _REQUIRED_ALIAS_SETS.values()
            )

            if has_all_required:
                profile = 'abi' if (norm_set & ABI_INDICATORS) else 'biorad'
                return idx, profile

    raise ValueError(f"Could not detect a header row in {file_path.name} (scanned first {max_scan_lines} lines).")


def load_qpcr_file(file_path: Path, config: PipelineConfig | None = None) -> pd.DataFrame:
    """
    Load a qPCR CSV, locating its header row and renaming instrument-specific
    columns to their canonical names (see config.COLUMN_ALIASES).

    Args:
        file_path: Path to qPCR CSV file
        config: Pipeline configuration (optional)

    Returns:
        DataFrame with canonical columns, normalized well IDs, and numeric
        Cq/quantity columns coerced to numbers

    Raises:
        ValueError: If no Well column is present after renaming
    """
    header_row, profile = find_qpcr_header_row(file_path, config=config)
    logger.info(f"Detected header at row {header_row} ({profile}) for file: {file_path.name}")

    # Read CSV with automatic delimiter detection
    df = pd.read_csv(file_path, header=header_row, sep=None, engine='python')

    # Clean up raw column names (strip + collapse spaces)
    df.columns = [re.sub(r'\s+', ' ', str(c)).strip() for c in df.columns]

    # Build a lowercase lookup for case-insensitive matching
    lower_cols = {c.lower(): c for c in df.columns}

    def present(*aliases):
        """Return the actual DF column name if any alias is present; else None."""
        for a in aliases:
            if a.lower() in lower_cols:
                return lower_cols[a.lower()]
        return None

    # Rename every recognized alias to its canonical column name.
    rename_map = {}
    for canonical, aliases in COLUMN_ALIASES.items():
        actual_col = present(*aliases)
        if actual_col:
            rename_map[actual_col] = canonical
    df = df.rename(columns=rename_map)

    # Ensure 'Well' exists after normalization
    if COL.WELL not in df.columns:
        logger.error(f"'{COL.WELL}' column not found after normalization in {file_path.name}")
        raise ValueError(f"Cannot proceed without a '{COL.WELL}' column.")

    # Normalize well IDs
    df[COL.WELL] = df[COL.WELL].apply(normalize_well)

    # Convert numeric columns
    numeric_cols = [COL.CQ, COL.CQ_MEAN_RAW, COL.CQ_SD_RAW, COL.SQ_RAW]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # Create Cq Std. Dev if missing (for ABI files)
    if COL.CQ_SD_RAW not in df.columns:
        df[COL.CQ_SD_RAW] = pd.NA

    return df


def _find_plate_start(row: pd.Series, ncols: int) -> int | None:
    """
    Return the column index of the run-name cell if this row looks like a
    plate-start row: some cell holding the literal value "1", immediately
    preceded by a non-empty (run-name) cell.

    Normally that's column 0 (name) / column 1 ("1"). Scanning every position
    instead of assuming column 0 means a plate map shifted right — e.g. by a
    stray leading blank column — is still recognized, just starting later.
    """
    for label_col in range(ncols - 1):
        name_cell = row.iloc[label_col]
        if not name_cell:
            continue
        try:
            if int(row.iloc[label_col + 1]) == 1:
                return label_col
        except (ValueError, TypeError):
            continue
    return None


def _extract_plate_wells(
    plate_df: pd.DataFrame,
    run_name: str,
    grid_start: int,
    label_col: int,
    config: PipelineConfig,
) -> list[dict[str, str]]:
    """
    Read one plate's row x column grid (rows A-H, columns 1-config.plate_cols,
    offset by label_col) starting at grid_start, returning one Run/Well/Sample
    dict per well.

    The row letter (A-H) is derived from its position in the grid rather than
    read from the plate map's own row-label cell, so a blank or missing label
    doesn't lose data — but a *present, disagreeing* label is still worth a
    warning, since that usually means something is misaligned.
    """
    ncols = plate_df.shape[1]
    rows = []

    for offset in range(config.plate_rows):
        r = grid_start + offset
        row_letter = chr(ord('A') + offset)

        actual_label = plate_df.iat[r, label_col] if label_col < ncols else ''
        if actual_label and actual_label.upper() != row_letter:
            logger.warning(
                f"⚠️ Run '{run_name}', plate map row {r}: expected row label '{row_letter}' "
                f"(based on position) but found '{actual_label}'. Using '{row_letter}' — "
                "please check this plate map section for missing or shifted rows."
            )

        for c in range(1, config.plate_cols + 1):
            col_idx = label_col + c
            sample = plate_df.iat[r, col_idx] if col_idx < ncols else ''
            sample = sample.strip() if sample else EMPTY_WELL_SAMPLE
            well = f"{row_letter}{c:02d}"
            rows.append({COL.RUN: run_name, COL.WELL: well, COL.SAMPLE: sample})

    return rows


def load_plate_map(plate_map_path: Path, config: PipelineConfig) -> pd.DataFrame:
    """
    Parse the plate map CSV (one or more grids stacked vertically, each
    preceded by a "<run name>,1,2,...,12" header row) into a long-format
    table of one row per well.

    Mechanical problems (a shifted grid, blank row labels) are recovered from
    automatically. A run name reused by more than one plate can't be safely
    guessed at, so it's loudly flagged instead: only the first occurrence is
    kept, since merging both would silently double-count every well shared
    between them against whichever single raw CSV matches that run name.

    Args:
        plate_map_path: Path to plate map CSV
        config: Pipeline configuration (plate dimensions)

    Returns:
        DataFrame with columns Run, Well, Sample (Sample is EMPTY_WELL_SAMPLE
        for any well left blank in the plate map)
    """
    plate_df = pd.read_csv(plate_map_path, header=None, dtype=str).fillna('')

    # Strip whitespace
    try:
        plate_df = plate_df.map(lambda x: x.strip() if isinstance(x, str) else '')
    except AttributeError:
        # Fallback for older pandas
        plate_df = plate_df.applymap(lambda x: x.strip() if isinstance(x, str) else '')

    transformed_rows = []
    nrows, ncols = plate_df.shape
    seen_runs: dict[str, int] = {}  # run name -> row index of its first plate
    i = 0

    while i < nrows:
        label_col = _find_plate_start(plate_df.iloc[i], ncols)

        if label_col is not None:
            run_name = plate_df.iat[i, label_col]
            grid_start = i + 1

            if label_col != 0:
                logger.warning(
                    f"⚠️ Plate map row {i} for run '{run_name}' is shifted {label_col} "
                    "column(s) to the right of where it's normally expected — "
                    "auto-correcting, but please check for a stray leading blank column."
                )

            if grid_start + config.plate_rows > nrows:
                logger.warning(f"Incomplete plate for run '{run_name}' starting at row {i}")
                break

            if run_name in seen_runs:
                logger.error(
                    f"❌ Duplicate run name '{run_name}' in plate_maps.csv: seen again at row {i} "
                    f"(first seen at row {seen_runs[run_name]}). Only the first occurrence will be "
                    "used — a raw CSV file can only match one of them, and keeping both would silently "
                    "double-count every shared well. Please rename one of these plate map sections "
                    "and re-run."
                )
            else:
                seen_runs[run_name] = i
                transformed_rows.extend(
                    _extract_plate_wells(plate_df, run_name, grid_start, label_col, config)
                )

            i = grid_start + config.plate_rows
            continue

        i += 1

    plate_map = pd.DataFrame(transformed_rows)
    logger.info(f"✅ Plate map transformed: {len(plate_map)} wells processed.")

    return plate_map
