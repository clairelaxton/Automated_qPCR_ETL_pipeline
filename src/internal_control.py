"""
Optional internal control (e.g. RNase P) correction (Step 3C): normalizes
SQ Mean against a user-chosen internal control target, so downstream
comparisons account for sample-to-sample input/extraction variability.

Disabled unless config.internal_control_target is set. The correction factor
for a (Sample, Fluorophore-derived Target) row is looked up by tier:
  1) internal control Cq for the same Sample + Batch (derived from Run)
  2) fallback: internal control Cq for the same Sample, any batch
  3) still missing: Correction Factor = 1.0 (no correction), flagged
"""
import re
import numpy as np
import pandas as pd
import logging

from .config import PipelineConfig, COL
from .quality_control import is_negative_control, is_positive_control
from .target_mapping import closest_match


logger = logging.getLogger(__name__)


def _is_control(sample_name: str, config: PipelineConfig) -> bool:
    """True for negative/positive control rows, which are excluded from internal-control lookups."""
    return is_negative_control(sample_name, config) or is_positive_control(sample_name, config)


def _extract_batch_tag(run_name: str, batch_regex: str | None) -> str | float:
    """
    Derive a "Batch" tag from a run name, used to pair an assay plate with
    its matching internal-control plate.

    Default: run names are "ASSAY_BATCH_DATE" (underscore-delimited); the
    batch tag is the middle token. If your run names follow a different
    convention, set config.internal_control_batch_regex to a pattern with
    one capturing group for the batch tag instead.
    """
    s = str(run_name)

    if batch_regex:
        m = re.search(batch_regex, s, flags=re.IGNORECASE)
        return m.group(1) if m else np.nan

    parts = s.split("_")
    return parts[1] if len(parts) >= 3 else np.nan


def _resolve_target_name(requested: str, summary: pd.DataFrame) -> str:
    """Case-insensitively match the requested internal control target against Targets actually present in summary."""
    available = sorted(summary[COL.TARGET].dropna().unique())
    matches = [t for t in available if t.upper() == requested.upper()]
    if matches:
        return matches[0]

    suggestion = closest_match(requested, available)
    hint = f" Did you mean '{suggestion}'?" if suggestion else ""
    raise ValueError(
        f"Internal control target '{requested}' was not found in your data.{hint} "
        f"Available targets: {', '.join(available) if available else '(none)'}. "
        "Check target_mapping.csv and --internal-control for spelling/case."
    )


def apply_internal_control_correction(summary: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    """
    Add Correction Factor / Corrected SQ Mean columns normalized against
    config.internal_control_target, if one is configured; otherwise a no-op.

    Args:
        summary: Summary DataFrame (post target-assignment, pre-collapse)
        config: Pipeline configuration

    Returns:
        summary unchanged if the feature is disabled; otherwise a copy with
        Correction Factor, Corrected SQ Mean, and IC_Missing columns added

    Raises:
        ValueError: If internal_control_target doesn't match any Target in
            the data, or no valid internal-control Cq values exist to build
            a baseline from
    """
    if not config.internal_control_target:
        return summary

    if summary.empty:
        return summary

    target_name = _resolve_target_name(config.internal_control_target, summary)
    summary = summary.copy()

    ic_mask = summary[COL.TARGET] == target_name
    is_sample_row = ~summary[COL.SAMPLE].apply(lambda s: _is_control(s, config))

    summary[COL.BATCH] = summary[COL.RUN].apply(
        lambda run: _extract_batch_tag(run, config.internal_control_batch_regex)
    )

    # Tiered internal-control Cq lookups: prefer Sample+Batch, fall back to Sample-only
    ic_rows = summary.loc[ic_mask & is_sample_row, [COL.SAMPLE, COL.BATCH, COL.MEAN_CQ]]
    ic_per_sample_batch = (
        ic_rows.groupby([COL.SAMPLE, COL.BATCH], as_index=False)
        .agg(IC_Cq_batch=(COL.MEAN_CQ, "mean"))
    )
    ic_per_sample_any = (
        ic_rows.groupby(COL.SAMPLE, as_index=False)
        .agg(IC_Cq_any=(COL.MEAN_CQ, "mean"))
    )

    # Global baseline: median Cq across all non-control internal-control rows
    baseline_vals = pd.to_numeric(ic_rows[COL.MEAN_CQ], errors="coerce")
    baseline_vals = baseline_vals[baseline_vals > 0]
    if baseline_vals.empty:
        raise ValueError(
            f"No valid internal control Cq values found for target '{target_name}'. "
            "Check that this target actually amplified in your data."
        )
    global_baseline = float(baseline_vals.median())
    logger.info(f"Internal control ('{target_name}') global baseline (median Cq) = {global_baseline:.3f}")

    summary = summary.merge(ic_per_sample_batch, on=[COL.SAMPLE, COL.BATCH], how="left")
    summary = summary.merge(ic_per_sample_any, on=COL.SAMPLE, how="left")

    summary[COL.IC_CQ] = summary["IC_Cq_batch"].fillna(summary["IC_Cq_any"])
    summary[COL.IC_SOURCE] = np.select(
        [summary["IC_Cq_batch"].notna(), summary["IC_Cq_batch"].isna() & summary["IC_Cq_any"].notna()],
        ["sample+batch", "sample-only"],
        default="missing",
    )
    summary = summary.drop(columns=["IC_Cq_batch", "IC_Cq_any"])

    summary[COL.IC_CQ] = pd.to_numeric(summary[COL.IC_CQ], errors="coerce")
    summary[COL.IC_GLOBAL_BASELINE] = global_baseline

    summary[COL.CORRECTION_FACTOR] = summary[COL.IC_CQ] / summary[COL.IC_GLOBAL_BASELINE]
    summary[COL.CORRECTION_FACTOR] = np.where(
        summary[COL.IC_CQ].notna() & (summary[COL.IC_CQ] > 0),
        summary[COL.CORRECTION_FACTOR],
        1.0,
    )
    summary[COL.CORRECTED_SQ_MEAN] = summary[COL.SQ_MEAN] * summary[COL.CORRECTION_FACTOR]
    summary[COL.IC_MISSING] = summary[COL.IC_CQ].isna()

    n_batch = int((summary[COL.IC_SOURCE] == "sample+batch").sum())
    n_sample = int((summary[COL.IC_SOURCE] == "sample-only").sum())
    n_missing = int((summary[COL.IC_SOURCE] == "missing").sum())
    log_fn = logger.warning if n_missing else logger.info
    log_fn(
        "Internal control matching counts — "
        f"sample+batch: {n_batch}, sample-only: {n_sample}, not matched: {n_missing}"
    )

    return summary.round({
        COL.IC_CQ: 3,
        COL.IC_GLOBAL_BASELINE: 3,
        COL.CORRECTION_FACTOR: 4,
        COL.CORRECTED_SQ_MEAN: 3,
    })
