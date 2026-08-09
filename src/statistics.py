"""
Numeric and aggregation helpers for qPCR data: replicate-level statistics
(mean/SD ignoring non-amplifying zero/negative Cq values), amplification
status aggregation, and the Run/Sample/Fluorophore summary table.
"""
import logging
import pandas as pd

from .config import COL, STATUS


logger = logging.getLogger(__name__)


def mean_positive_only(series: pd.Series) -> float:
    """
    Calculate the mean of positive (>0) values only, treating non-positive
    values (e.g. a non-amplifying replicate's Cq of 0) as not contributing.

    Args:
        series: Series of numeric values

    Returns:
        Mean of positive values, or 0.0 if none are positive
    """
    vals = pd.to_numeric(series, errors='coerce')
    vals = vals[vals > 0]
    return float(vals.mean()) if not vals.empty else 0.0


def sd_positive_only(series: pd.Series) -> float:
    """
    Calculate the population standard deviation (ddof=0) of positive (>0)
    values only.

    Args:
        series: Series of numeric values

    Returns:
        Standard deviation of positive values, or 0.0 if fewer than 2 are positive
    """
    vals = pd.to_numeric(series, errors='coerce')
    vals = vals[vals > 0]
    return float(vals.std(ddof=0)) if len(vals) >= 2 else 0.0


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    """
    Calculate a replicate-count-weighted mean, falling back to a simple mean
    if no weights are usable (e.g. all zero or NaN).

    Args:
        values: Series of values to average
        weights: Series of weights (typically replicate counts)

    Returns:
        Weighted mean as float
    """
    v = pd.to_numeric(values, errors="coerce")
    w = pd.to_numeric(weights, errors="coerce").fillna(0)
    mask = (w > 0) & v.notna()

    if mask.any():
        return float((v[mask] * w[mask]).sum() / w[mask].sum())

    return float(v.mean()) if v.notna().any() else 0.0


def aggregate_amp_status(series: pd.Series) -> str:
    """
    Aggregate per-replicate amplification status into one call for the group.

    Args:
        series: Series of STATUS.AMP/STATUS.NOAMP values, one per replicate

    Returns:
        STATUS.AMP or STATUS.NOAMP if all replicates agree, STATUS.INCONCLUSIVE
        if they disagree, or STATUS.NODATA if the group has no values
    """
    vals = series.dropna().astype(str)
    if vals.empty:
        return STATUS.NODATA

    n_amp = (vals == STATUS.AMP).sum()
    n_noamp = (vals == STATUS.NOAMP).sum()

    if n_amp > 0 and n_noamp > 0:
        return STATUS.INCONCLUSIVE

    return STATUS.AMP if n_amp == len(vals) else STATUS.NOAMP


def _round_stat_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Round the three computed statistics columns (Mean Cq, Cq SD, SQ Mean) to 3 decimal places."""
    for col in (COL.MEAN_CQ, COL.CQ_SD, COL.SQ_MEAN):
        df[col] = df[col].round(3)
    return df


def compute_replicate_stats(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    """
    Add per-group replicate statistics to a raw (one row per well) DataFrame.

    Args:
        df: DataFrame with raw Cq and starting-quantity values
        group_cols: Columns identifying a replicate group (e.g. Sample, Fluor)

    Returns:
        Copy of df with added COL.MEAN_CQ, COL.CQ_SD, COL.SQ_MEAN columns,
        each broadcast across every row of its group and rounded to 3 places
    """
    df = df.copy()

    # Compute statistics by group
    df[COL.MEAN_CQ] = df.groupby(group_cols)[COL.CQ].transform(mean_positive_only)
    df[COL.CQ_SD] = df.groupby(group_cols)[COL.CQ].transform(sd_positive_only)
    df[COL.SQ_MEAN] = df.groupby(group_cols)[COL.SQ_RAW].transform(mean_positive_only)

    return _round_stat_columns(df)


def create_summary_table(combined: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse per-well rows from every run into one row per
    Run/Sample/Fluorophore, with aggregated Cq/quantity statistics, the set
    of contributing wells, a replicate count, and an aggregated amplification
    status.

    Args:
        combined: Combined qPCR data across all runs (one row per well)

    Returns:
        Summary DataFrame: one row per Run/Sample/Fluorophore
    """
    group_cols = [COL.RUN, COL.SAMPLE, COL.FLUOROPHORE]

    # Well positions per group
    well_positions = (
        combined
        .groupby(group_cols, dropna=False)[COL.WELL]
        .apply(lambda s: ','.join(sorted(map(str, set(s)))))
        .reset_index(name=COL.WELL)
    )

    # Replicate count per group
    rep_counts = (
        combined
        .groupby(group_cols, dropna=False)
        .size()
        .reset_index(name=COL.REPLICATES)
    )

    # Aggregate stats from raw values (named aggregation: output_name=(column, func))
    agg_stats = (
        combined
        .groupby(group_cols, dropna=False)
        .agg(**{
            'Mean_Cq': (COL.CQ, mean_positive_only),
            'Cq_SD': (COL.CQ, sd_positive_only),
            'SQ_Mean': (COL.SQ_RAW, mean_positive_only),
            COL.AMP_STATUS: (COL.AMPLIFICATION_STATUS, aggregate_amp_status),
        })
        .reset_index()
        .rename(columns={'Mean_Cq': COL.MEAN_CQ, 'Cq_SD': COL.CQ_SD, 'SQ_Mean': COL.SQ_MEAN})
    )

    # Merge components
    summary = (
        agg_stats
        .merge(well_positions, on=group_cols, how='left')
        .merge(rep_counts, on=group_cols, how='left')
    )
    summary = _round_stat_columns(summary)

    # Order columns
    summary = summary[
        [COL.RUN, COL.SAMPLE, COL.FLUOROPHORE, COL.WELL,
         COL.REPLICATES, COL.MEAN_CQ, COL.CQ_SD, COL.SQ_MEAN,
         COL.AMP_STATUS]
    ]

    logger.info(f"✅ Summary created ({len(summary)} rows).")

    return summary
