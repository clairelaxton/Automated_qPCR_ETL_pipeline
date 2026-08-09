"""
Quality control for qPCR data: per-row flagging plus inter-run collapsing.

These two concerns live in one module because they're two steps of a single
QC phase (see QPCRPipeline.run(): flag_samples() then collapse_interrun()),
and collapsing calls the control-detection helpers defined here directly.
"""
import re
import logging
import pandas as pd
import numpy as np

from .config import (
    PipelineConfig, EXCLUSIONS_COLUMNS, COL, STATUS, FLAG, REASON, REPEAT,
)
from .statistics import weighted_mean


logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Per-row flagging
# ---------------------------------------------------------------------------

def _name_matches_keywords(name: str, keywords: set[str]) -> bool:
    """
    True if any alphanumeric token in name matches a keyword exactly, or
    matches a keyword plus trailing digits (e.g. 'NTC1', 'POS3'). Tokens are
    split on runs of non-alphanumeric characters, e.g. 'NEG CTRL-1' -> ['NEG', 'CTRL', '1'].
    """
    tokens = (t for t in re.split(r"[^A-Z0-9]+", str(name).upper()) if t)
    return any(
        token == kw or (token.startswith(kw) and token[len(kw):].isdigit())
        for token in tokens
        for kw in keywords
    )


def is_negative_control(name: str, config: PipelineConfig) -> bool:
    """Return True if name matches one of config's negative control keywords (NTC, NEC, ...)."""
    return _name_matches_keywords(name, config.negative_control_keywords)


def is_positive_control(name: str, config: PipelineConfig) -> bool:
    """Return True if name matches one of config's positive control keywords (POS, PC, ...)."""
    return _name_matches_keywords(name, config.positive_control_keywords)


def collect_flags(row: pd.Series, config: PipelineConfig) -> str:
    """
    Compute the comma-separated QC flags for one summary row.

    Negative/positive controls are checked first and flagged only for
    unexpected amplification behavior (a control returns early with just
    that verdict). Everything else is checked for single-replicate calls,
    mixed-replicate (inconclusive) results, unexpectedly low Cq, high
    replicate variability, and a missing internal control match — any
    number of which can apply at once.

    Args:
        row: A row from the summary table (Sample, Mean Cq, Cq SD,
            Replicates, AmpStatus_Aggregated, and IC_Missing if internal
            control correction ran)
        config: Pipeline configuration (thresholds and control keywords)

    Returns:
        Comma-and-space-joined flag labels, or '' if none apply
    """
    flags = []

    # Extract needed fields with safe defaults
    name = str(row.get(COL.SAMPLE, '')).upper()
    cq = row.get(COL.MEAN_CQ, np.nan)
    sd = row.get(COL.CQ_SD, np.nan)
    reps = int(row.get(COL.REPLICATES, 0)) if pd.notna(row.get(COL.REPLICATES, np.nan)) else 0
    agg = str(row.get(COL.AMP_STATUS, '')).strip()

    # ---- Control-specific checks (handled and returned early) ----
    if is_negative_control(name, config):
        if agg == STATUS.AMP or (pd.notna(cq) and cq > 0):
            if pd.notna(cq) and cq > 0:
                bucket = "<35" if cq < config.pos_cq_strong else f">={int(config.pos_cq_strong)}"
                flags.append(f"{FLAG.NTC_NEC_AMP}_Cq{bucket}")
            else:
                flags.append(FLAG.NTC_NEC_AMP)
        return ", ".join(flags)

    if is_positive_control(name, config):
        if agg != STATUS.AMP or not (pd.notna(cq) and cq > 0):
            flags.append(FLAG.POS_CONTROL_INCONCLUSIVE if agg == STATUS.INCONCLUSIVE else FLAG.POS_CONTROL_NO_AMP)
        return ", ".join(flags)

    # ---- Single-replicate cases ----
    if reps == 1:
        if agg == STATUS.AMP:
            flags.append(FLAG.SINGLE_REPLICATE_AMP)
        elif agg == STATUS.NOAMP:
            flags.append(FLAG.SINGLE_REPLICATE_NOAMP)

    # ---- Inconclusive groups: mixed Amp/NoAmp ----
    if agg == STATUS.INCONCLUSIVE:
        if pd.notna(cq) and cq > 0:
            flags.append(FLAG.INCONCLUSIVE_LOW_CQ if cq < config.pos_cq_strong else FLAG.INCONCLUSIVE_HIGH_CQ)
        else:
            flags.append(FLAG.INCONCLUSIVE)

    # ---- Low Cq (unexpectedly low) ----
    if pd.notna(cq) and (config.low_cq_min < cq < config.low_cq_max):
        flags.append(FLAG.LOW_CQ)

    # ---- High variation across replicates ----
    if pd.notna(sd) and sd >= config.high_sd:
        flags.append(FLAG.HIGH_SD)

    # ---- Missing internal control (only present if that correction ran) ----
    if bool(row.get(COL.IC_MISSING, False)):
        flags.append(FLAG.MISSING_INTERNAL_CONTROL)

    return ", ".join(flags)


def contains_any_bad_flag(flags_cell: str, bad_flags: set[str]) -> bool:
    """Return True if flags_cell (a comma-separated Flags value) contains any flag in bad_flags."""
    if not isinstance(flags_cell, str) or not flags_cell.strip():
        return False
    parts = {p.strip() for p in flags_cell.split(",") if p.strip()}
    return any(bad in parts for bad in bad_flags)


def flag_samples(summary: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    """
    Add a Flags column to the summary table by running collect_flags() over every row.

    Args:
        summary: Summary DataFrame (one row per Run/Sample/Fluorophore)
        config: Pipeline configuration

    Returns:
        Copy of summary with a 'Flags' column added
    """
    summary = summary.copy()
    summary[COL.FLAGS] = summary.apply(lambda row: collect_flags(row, config), axis=1)
    n_flagged = (summary[COL.FLAGS].astype(str).str.len() > 0).sum()
    logger.info(f"🚩 Flagging complete. Rows with any flags: {n_flagged} / {len(summary)}")
    return summary


# ---------------------------------------------------------------------------
# Inter-run collapsing
# ---------------------------------------------------------------------------

def _append_excluded(
    rows: pd.DataFrame,
    reason: str,
    repeat_needed: bool,
    sink: list[pd.DataFrame],
    sample: str,
    target: str
) -> None:
    """
    Record rows dropped during collapsing into the exclusions accumulator.

    Args:
        rows: Rows being excluded (no-op if empty)
        reason: One of the REASON.* exclusion reasons
        repeat_needed: Whether these rows should be flagged for a repeat
        sink: List accumulating one exclusions DataFrame per call
        sample: Sample name (shared by every row)
        target: Target name (shared by every row)
    """
    if rows.empty:
        return

    sink.append(pd.DataFrame({
        COL.SAMPLE: sample,
        COL.TARGET: target,
        COL.RUN: rows[COL.RUN].astype(str),
        COL.WELL: rows[COL.WELL].astype(str),
        COL.REPLICATES: pd.to_numeric(rows[COL.REPLICATES], errors="coerce").fillna(0).astype(int),
        COL.AMP_STATUS: rows[COL.AMP_STATUS].astype(str),
        COL.FLAGS: rows[COL.FLAGS].astype(str),
        COL.EXCLUDE_REASON: reason,
        COL.REPEAT_NEEDED: REPEAT.YES if repeat_needed else REPEAT.NO,
    }))


def _resolve_majority(
    grp: pd.DataFrame,
    good_rows: pd.DataFrame,
    winning_status: str,
    reason: str,
    excluded_accum: list[pd.DataFrame],
    sample: str,
    target: str,
) -> pd.DataFrame:
    """
    Keep the good rows agreeing with winning_status and record every other
    row in grp (both the disagreeing good rows and any bad-flagged ones) as
    excluded under reason.
    """
    kept = good_rows.loc[good_rows[COL.AMP_STATUS] == winning_status].copy()
    dropped = grp.loc[grp[COL.AMP_STATUS] != winning_status].copy()
    _append_excluded(dropped, reason=reason, repeat_needed=False, sink=excluded_accum, sample=sample, target=target)
    return kept


def _aggregate_kept_rows(kept: pd.DataFrame, sample: str, target: str, final_call: str) -> pd.DataFrame:
    """Combine the rows kept for one (Sample, Target) into a single, replicate-weighted result row."""
    reps = pd.to_numeric(kept[COL.REPLICATES], errors="coerce").fillna(0)
    fields = {
        COL.SAMPLE: sample,
        COL.TARGET: target,
        COL.RUN: ",".join(sorted(set(map(str, kept[COL.RUN])))),
        COL.WELL: ",".join(sorted(set(map(str, kept[COL.WELL])))),
        COL.REPLICATES: int(reps.sum()),
        COL.AMP_STATUS: final_call,
        COL.MEAN_CQ: weighted_mean(kept[COL.MEAN_CQ], reps),
        COL.CQ_SD: weighted_mean(kept[COL.CQ_SD], reps),
        COL.SQ_MEAN: weighted_mean(kept[COL.SQ_MEAN], reps),
        COL.FLAGS: ",".join(sorted({
            f for f in ",".join(kept[COL.FLAGS].astype(str)).split(",") if f.strip()
        })),
    }

    # Only present when internal control correction ran; carry them forward if so.
    for col in (COL.CORRECTION_FACTOR, COL.CORRECTED_SQ_MEAN):
        if col in kept.columns:
            fields[col] = weighted_mean(kept[col], reps)

    out = pd.Series(fields)
    for col in (COL.MEAN_CQ, COL.CQ_SD, COL.SQ_MEAN, COL.CORRECTION_FACTOR, COL.CORRECTED_SQ_MEAN):
        if col in out.index and pd.notna(out[col]):
            out[col] = round(float(out[col]), 3)

    return out.to_frame().T


def _collapse_interrun(
    grp: pd.DataFrame,
    excluded_accum: list[pd.DataFrame],
    config: PipelineConfig,
    sample: str | None = None,
    target: str | None = None,
) -> pd.DataFrame:
    """
    Collapse every run's row for one (Sample, Target) into a single result
    row by majority vote: bad-flagged rows are set aside first, then the
    remaining rows' amplification calls decide the outcome — a clear
    majority wins (the minority is dropped), a tie or an all-Inconclusive
    group excludes everything and marks it for a repeat.

    Args:
        grp: Rows for one (Sample, Target) across all runs
        excluded_accum: List to append any excluded rows to
        config: Pipeline configuration (bad_flags)
        sample: Sample name (read from grp if not given)
        target: Target name (read from grp if not given)

    Returns:
        Single-row DataFrame with the collapsed result, or an empty
        DataFrame if everything was excluded
    """
    if sample is None:
        if COL.SAMPLE not in grp.columns:
            raise KeyError(COL.SAMPLE)
        sample = grp[COL.SAMPLE].iloc[0]
    if target is None:
        if COL.TARGET not in grp.columns:
            raise KeyError(COL.TARGET)
        target = grp[COL.TARGET].iloc[0]

    # Separate bad-flag rows
    mask_bad = grp[COL.FLAGS].apply(lambda x: contains_any_bad_flag(x, config.bad_flags))
    good_rows = grp.loc[~mask_bad].copy()

    # Count statuses among good rows
    status = good_rows[COL.AMP_STATUS].astype(str)
    n_amp = (status == STATUS.AMP).sum()
    n_noamp = (status == STATUS.NOAMP).sum()

    # Decide by majority
    if n_amp == 0 and n_noamp == 0:
        # All Inconclusive (or nothing survived flag-filtering): exclude all, repeat needed
        _append_excluded(grp, reason=REASON.ALL_INCONCLUSIVE, repeat_needed=True,
                          sink=excluded_accum, sample=sample, target=target)
        return pd.DataFrame()

    if n_amp > n_noamp:
        final_call = STATUS.AMP
        kept = _resolve_majority(grp, good_rows, STATUS.AMP, REASON.MAJORITY_AMP, excluded_accum, sample, target)
    elif n_noamp > n_amp:
        final_call = STATUS.NOAMP
        kept = _resolve_majority(grp, good_rows, STATUS.NOAMP, REASON.MAJORITY_NOAMP, excluded_accum, sample, target)
    else:
        # Tie: exclude all, repeat needed
        _append_excluded(grp, reason=REASON.TIE, repeat_needed=True,
                          sink=excluded_accum, sample=sample, target=target)
        return pd.DataFrame()

    return _aggregate_kept_rows(kept, sample, target, final_call)


def collapse_interrun(
    summary: pd.DataFrame,
    config: PipelineConfig
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Set controls aside, then collapse every remaining (Sample, Target)
    group's per-run rows into one row using majority-vote rules (see
    _collapse_interrun).

    Args:
        summary: Summary DataFrame, potentially with multiple runs per
            (Sample, Target)
        config: Pipeline configuration

    Returns:
        Tuple of (final_summary, excluded_df, controls_excluded)
    """
    # Remove controls before collapsing
    controls_mask = summary[COL.SAMPLE].apply(
        lambda s: is_negative_control(s, config) or is_positive_control(s, config)
    )
    controls_excluded = summary[controls_mask].copy()
    summary_no_controls = summary[~controls_mask].copy()

    logger.info(
        f"🧪 Removed {len(controls_excluded)} control rows "
        "(NTC/NEC/POS) before inter-run collapsing."
    )

    # Collapse by (Sample, Target)
    excluded_rows = []

    collapsed_rows = []
    for (sample, target), grp in summary_no_controls.groupby([COL.SAMPLE, COL.TARGET]):
        collapsed = _collapse_interrun(
            grp=grp,
            excluded_accum=excluded_rows,
            config=config,
            sample=sample,
            target=target,
        )
        if not collapsed.empty:
            collapsed_rows.append(collapsed)

    final_summary = (
        pd.concat(collapsed_rows, ignore_index=True)
        if collapsed_rows
        else pd.DataFrame()
    )

    excluded_df = (
        pd.concat(excluded_rows, ignore_index=True)
        if excluded_rows
        else pd.DataFrame(columns=list(EXCLUSIONS_COLUMNS))
    )

    logger.info(
        f"🧮 Inter-run collapse complete. "
        f"Final rows: {len(final_summary)}; Excluded rows: {len(excluded_df)}"
    )

    return final_summary, excluded_df, controls_excluded
