"""
Per-file qPCR processing: load one instrument export, attach sample names
from the plate map, compute replicate statistics, and call amplification
status per well.
"""
import pandas as pd
from pathlib import Path
import logging

from .file_io import load_qpcr_file
from .config import PipelineConfig, EMPTY_WELL_SAMPLE, COL, STATUS
from .statistics import compute_replicate_stats


logger = logging.getLogger(__name__)


def process_qpcr_file(
    file_path: Path,
    plate_map: pd.DataFrame,
    config: PipelineConfig
) -> tuple[pd.DataFrame, bool]:
    """
    Load one qPCR file, merge in its plate-map sample names, and compute
    replicate statistics and amplification status per well.

    Args:
        file_path: Path to qPCR CSV file (its stem is treated as the run name)
        plate_map: Long-format plate map DataFrame (Run, Well, Sample)
        config: Pipeline configuration

    Returns:
        Tuple of (processed DataFrame, success boolean). On any failure
        (unreadable file, no matching plate-map rows, missing required
        columns) returns (empty DataFrame, False) so the caller can skip it.
    """
    run_name = file_path.stem
    logger.info(f"📄 Processing {file_path.name} (Run: {run_name})")

    try:
        # Load qPCR file
        qpcr_data = load_qpcr_file(file_path, config)
    except Exception as e:
        logger.error(f"❌ Error reading {file_path.name}: {e}")
        return pd.DataFrame(), False

    # Drop instrument's Sample column (use plate map instead)
    if COL.SAMPLE in qpcr_data.columns:
        qpcr_data = qpcr_data.drop(columns=[COL.SAMPLE])

    # Filter plate map for this run
    run_plate_map = plate_map[plate_map[COL.RUN] == run_name]
    if run_plate_map.empty:
        logger.warning(
            f"Run '{run_name}' not found in plate map. "
            f"All rows will be labeled '{EMPTY_WELL_SAMPLE}' and then dropped."
        )

    # Merge with plate map
    merged = qpcr_data.merge(
        run_plate_map[[COL.WELL, COL.SAMPLE]],
        on=COL.WELL,
        how='left'
    )

    # Replace missing samples with the empty-well sentinel
    merged[COL.SAMPLE] = merged[COL.SAMPLE].fillna(EMPTY_WELL_SAMPLE)

    # Log and drop empty wells
    empty_count = int((merged[COL.SAMPLE] == EMPTY_WELL_SAMPLE).sum())
    rows_before = len(merged)
    logger.info(
        f"🧪 Ignoring {empty_count} empty wells "
        f"out of {rows_before} rows in {file_path.name}"
    )

    merged = merged[merged[COL.SAMPLE] != EMPTY_WELL_SAMPLE]
    if merged.empty:
        logger.warning(
            f"⚠️ All rows were empty after plate-map merge for run '{run_name}'. "
            "Skipping file."
        )
        return pd.DataFrame(), False

    # Validate required columns
    if COL.FLUOR not in merged.columns:
        logger.error(
            f"'{COL.FLUOR}' column missing in {file_path.name} after normalization. Skipping."
        )
        return pd.DataFrame(), False

    required_numeric = [COL.SQ_RAW, COL.CQ]
    missing = [col for col in required_numeric if col not in merged.columns]
    if missing:
        logger.error(
            f"Required numeric columns missing in {file_path.name}: {missing}. Skipping."
        )
        return pd.DataFrame(), False

    # Ensure numeric types
    merged[COL.CQ] = pd.to_numeric(merged[COL.CQ], errors='coerce')
    merged[COL.SQ_RAW] = pd.to_numeric(merged[COL.SQ_RAW], errors='coerce')

    # Compute replicate statistics
    merged = compute_replicate_stats(merged, [COL.SAMPLE, COL.FLUOR])

    # Add run name
    merged[COL.RUN] = run_name

    # Build output
    output = merged[[
        COL.SAMPLE, COL.FLUOR, COL.WELL, COL.RUN,
        COL.CQ, COL.SQ_RAW,
        COL.MEAN_CQ, COL.CQ_SD, COL.SQ_MEAN,
    ]].copy()

    # Rename Fluor to Fluorophore
    output = output.rename(columns={COL.FLUOR: COL.FLUOROPHORE})

    # Add amplification status
    output[COL.AMPLIFICATION_STATUS] = output[COL.CQ].apply(
        lambda x: STATUS.AMP if pd.notna(x) and x > 0 else STATUS.NOAMP
    )

    return output, True


def process_all_qpcr_files(
    qpcr_files: list[Path],
    plate_map: pd.DataFrame,
    config: PipelineConfig
) -> tuple[pd.DataFrame, list[str]]:
    """
    Process every qPCR file and concatenate the successful ones.

    Args:
        qpcr_files: List of qPCR file paths
        plate_map: Plate map DataFrame
        config: Pipeline configuration

    Returns:
        Tuple of (combined DataFrame across all files, list of skipped file names)
    """
    all_data = []
    skipped_files = []

    for file_path in qpcr_files:
        output, success = process_qpcr_file(file_path, plate_map, config)

        if success:
            all_data.append(output)
        else:
            skipped_files.append(file_path.name)

    if not all_data:
        logger.error("⚠️ No valid qPCR data was processed.")
        return pd.DataFrame(), skipped_files

    combined = pd.concat(all_data, ignore_index=True)
    logger.info(
        f"📈 Processed {len(all_data)} file(s); total rows: {len(combined)}"
    )

    if skipped_files:
        logger.warning(f"⚠️ Skipped {len(skipped_files)} file(s): {', '.join(skipped_files)}")

    return combined, skipped_files
