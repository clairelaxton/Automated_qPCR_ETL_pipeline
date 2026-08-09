"""
Target assignment: maps each (Run, Fluorophore) pair to a target name using
target_mapping.csv, matching runs by their longest configured prefix and
falling back to a DEFAULT/ANY entry (or "Unknown") when nothing matches.
"""
import re
import difflib
import pandas as pd
import logging

from .config import (
    PipelineConfig, COL,
    DEFAULT_RUN_PREFIX, ANY_FLUOROPHORE, UNKNOWN_TARGET,
)


logger = logging.getLogger(__name__)


def closest_match(name: str, candidates: list[str]) -> str | None:
    """
    Best-guess suggestion for a likely-mistyped name, for use in warnings.

    Prefers a candidate that's a prefix/suffix of name (or vice versa) — the
    common case of a truncated or extended identifier, e.g. a run name
    missing its date suffix — falling back to general fuzzy similarity.
    Returns None if nothing looks close enough to suggest.
    """
    prefix_matches = [c for c in candidates if c.startswith(name) or name.startswith(c)]
    if prefix_matches:
        return min(prefix_matches, key=len)
    close = difflib.get_close_matches(name, candidates, n=1)
    return close[0] if close else None


def load_target_mapping(config: PipelineConfig) -> dict[str, dict[str, str]]:
    """
    Load target_mapping.csv into a lookup dict, or fall back to a single
    DEFAULT/ANY -> Unknown entry if the file doesn't exist.

    Blank rows (e.g. trailing rows left over from editing in Excel) are
    dropped silently. A (RunPrefix, Fluorophore) pair mapped to more than one
    Target is a real ambiguity that can't be silently resolved, so it's
    logged as a warning — whichever row appears last in the file wins.

    Args:
        config: Pipeline configuration

    Returns:
        Dictionary mapping RunPrefix -> {Fluorophore -> Target}
    """
    target_map_file = config.target_mapping_path

    if target_map_file.exists():
        target_df = pd.read_csv(target_map_file, dtype=str).fillna("")
        logger.info(f"✅ Loaded target mapping from {target_map_file.name}")
    else:
        logger.warning(f"⚠️ No {target_map_file.name} found — using built-in defaults.")
        target_df = pd.DataFrame({
            COL.RUN_PREFIX: [DEFAULT_RUN_PREFIX],
            COL.FLUOROPHORE: [ANY_FLUOROPHORE],
            COL.TARGET: [UNKNOWN_TARGET],
        })

    # Normalize mapping table
    target_df[COL.RUN_PREFIX] = target_df[COL.RUN_PREFIX].str.upper().str.strip()
    target_df[COL.FLUOROPHORE] = target_df[COL.FLUOROPHORE].str.upper().str.strip()
    target_df[COL.TARGET] = target_df[COL.TARGET].str.strip()

    # Silently drop blank rows (RunPrefix left empty)
    target_df = target_df[target_df[COL.RUN_PREFIX] != ""]

    _warn_about_conflicting_rows(target_df)

    # Build lookup dictionary
    target_lookup = (
        target_df
        .groupby(COL.RUN_PREFIX, dropna=False)
        .apply(lambda g: dict(zip(g[COL.FLUOROPHORE], g[COL.TARGET])))
        .to_dict()
    )

    return target_lookup


def _warn_about_conflicting_rows(target_df: pd.DataFrame) -> None:
    """Warn about (RunPrefix, Fluorophore) pairs mapped to more than one Target — an unresolvable ambiguity."""
    dupe_mask = target_df.duplicated(subset=[COL.RUN_PREFIX, COL.FLUOROPHORE], keep=False)
    for (rp, fluor), grp in target_df.loc[dupe_mask].groupby([COL.RUN_PREFIX, COL.FLUOROPHORE]):
        targets = list(grp[COL.TARGET].unique())
        if len(targets) > 1:
            logger.warning(
                f"⚠️ target_mapping.csv maps ({rp}, {fluor}) to multiple different targets "
                f"({', '.join(targets)}) — only the last one ('{targets[-1]}') will be used. "
                "Please remove the duplicate row."
            )


def get_run_prefix(run_name: str, target_lookup: dict[str, dict[str, str]]) -> str:
    """
    Find the longest configured RunPrefix that run_name starts with.

    Args:
        run_name: Run name to match
        target_lookup: Target lookup dictionary

    Returns:
        The longest matching RunPrefix, or DEFAULT_RUN_PREFIX if none match
    """
    run_upper = str(run_name).upper()
    prefixes = [
        rp for rp in target_lookup.keys()
        if rp != DEFAULT_RUN_PREFIX and run_upper.startswith(rp)
    ]
    return max(prefixes, key=len) if prefixes else DEFAULT_RUN_PREFIX


def assign_target(row: pd.Series, target_lookup: dict[str, dict[str, str]]) -> str:
    """
    Assign a target name for one row, trying its RunPrefix's exact
    Fluorophore mapping, then that prefix's ANY entry, then the DEFAULT
    prefix's exact and ANY entries, before giving up.

    Args:
        row: Row with 'Run' and 'Fluorophore' columns
        target_lookup: Target lookup dictionary

    Returns:
        Target name, or UNKNOWN_TARGET if nothing matches
    """
    run = str(row.get(COL.RUN, "")).upper()
    fluor = str(row.get(COL.FLUOROPHORE, "")).upper().strip()

    rp = get_run_prefix(run, target_lookup)

    # Try exact RunPrefix mapping
    mapping = target_lookup.get(rp, {})
    if fluor in mapping:
        return mapping[fluor]
    if ANY_FLUOROPHORE in mapping:
        return mapping[ANY_FLUOROPHORE]

    # Try DEFAULT fallback
    default_mapping = target_lookup.get(DEFAULT_RUN_PREFIX, {})
    if fluor in default_mapping:
        return default_mapping[fluor]
    if ANY_FLUOROPHORE in default_mapping:
        return default_mapping[ANY_FLUOROPHORE]

    return UNKNOWN_TARGET


def check_unmatched_run_prefixes(summary: pd.DataFrame, target_lookup: dict[str, dict[str, str]]) -> None:
    """
    Warn about RunPrefix values in target_mapping.csv that never match any
    observed run name — almost always a typo, since that RunPrefix's rules
    are then silently dead (every run that should have used them falls back
    to DEFAULT instead, usually landing on 'Unknown').

    Args:
        summary: Summary DataFrame
        target_lookup: Target lookup dictionary
    """
    if summary.empty:
        return

    observed_runs = [str(r) for r in summary[COL.RUN].dropna().unique()]
    configured_prefixes = [rp for rp in target_lookup if rp != DEFAULT_RUN_PREFIX]
    unmatched = [
        rp for rp in configured_prefixes
        if not any(run.upper().startswith(rp) for run in observed_runs)
    ]
    if not unmatched:
        return

    # Compare against each run's leading alphanumeric segment (its likely
    # intended prefix) rather than the full run name, for a fair fuzzy match.
    observed_prefix_candidates = sorted({
        m.group(0) for m in (re.match(r'^[A-Z0-9]+', run.upper()) for run in observed_runs) if m
    })

    msg_lines = [
        "The following RunPrefix value(s) in target_mapping.csv don't match any run "
        "name in your data, so their rows are never used (everything falls back to "
        "DEFAULT instead):"
    ]
    for rp in unmatched:
        suggestion = closest_match(rp, observed_prefix_candidates)
        hint = f" — did you mean '{suggestion}'?" if suggestion else ""
        msg_lines.append(f"  - {rp}{hint}")
    msg_lines.append("→ Check target_mapping.csv for a typo in RunPrefix.")
    logger.warning("\n".join(msg_lines))


def check_missing_mappings(summary: pd.DataFrame, target_lookup: dict[str, dict[str, str]]) -> None:
    """
    Log a warning listing any (RunPrefix, Fluorophore) pair observed in the
    data that resolves to UNKNOWN_TARGET, so it's easy to see what
    target_mapping.csv doesn't cover. Reuses assign_target() itself (rather
    than re-deriving the same fallback chain) so this can never disagree
    with what the pipeline actually assigns.

    Args:
        summary: Summary DataFrame
        target_lookup: Target lookup dictionary
    """
    if summary.empty:
        return

    observed_pairs = (
        summary[[COL.RUN, COL.FLUOROPHORE]]
        .dropna()
        .drop_duplicates()
    )

    missing_pairs = [
        (get_run_prefix(run, target_lookup), str(fluor).upper().strip())
        for run, fluor in observed_pairs.itertuples(index=False, name=None)
        if assign_target(pd.Series({COL.RUN: run, COL.FLUOROPHORE: fluor}), target_lookup) == UNKNOWN_TARGET
    ]

    if not missing_pairs:
        logger.info("🔎 Target mapping QC: all observed (RunPrefix, Fluorophore) pairs are covered.")
        return

    msg_lines = [
        "The following (RunPrefix, Fluorophore) pairs are not mapped in target_mapping.csv:"
    ]
    for rp, flu in sorted(set(missing_pairs)):
        # Configured fluorophores for this prefix (or DEFAULT) that ended up unused —
        # a likely typo target, e.g. 'HOX' configured but 'ROX' is what's in the data.
        configured_fluors = list(target_lookup.get(rp, target_lookup.get(DEFAULT_RUN_PREFIX, {})).keys())
        suggestion = closest_match(flu, configured_fluors)
        hint = f" — did you mean '{suggestion}'?" if suggestion else ""
        msg_lines.append(f"  - {rp}, {flu}{hint}")
    msg_lines += [
        "→ To fix: open target_mapping.csv and add rows for each pair, e.g.:",
        "    RunPrefix,Fluorophore,Target",
        "    RP-SC2,FAM,SC2",
        "    ECHH,ROX,HHV6",
        "  (You can also use: DEFAULT,ANY,Unknown as a generic fallback.)"
    ]
    logger.warning("\n".join(msg_lines))


def apply_target_mapping(summary: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    """
    Load target_mapping.csv, warn about any configuration issues, and add a
    Target column to the summary table.

    Args:
        summary: Summary DataFrame
        config: Pipeline configuration

    Returns:
        Copy of summary with a 'Target' column added
    """
    target_lookup = load_target_mapping(config)

    check_unmatched_run_prefixes(summary, target_lookup)
    check_missing_mappings(summary, target_lookup)

    # Apply mapping
    summary = summary.copy()
    summary[COL.TARGET] = summary.apply(
        lambda row: assign_target(row, target_lookup),
        axis=1
    )

    if not summary.empty and (summary[COL.TARGET] == UNKNOWN_TARGET).all():
        logger.warning(
            "⚠️ Every sample resolved to Target='Unknown'. This usually means "
            "target_mapping.csv hasn't been filled in with real target names yet — "
            "please check it has rows beyond the DEFAULT/ANY fallback."
        )

    logger.info("🧬 Target labels assigned using target_mapping.csv")

    return summary
