"""
Configuration for the qPCR pipeline: the PipelineConfig settings dataclass,
plus every shared string constant (column names, status/flag values, output
file and sheet names, instrument column aliases). Centralizing these means a
value only needs to change in one place, and every other module imports the
namespace (COL, STATUS, FLAG, ...) it needs instead of retyping literals.
"""
from pathlib import Path
from dataclasses import dataclass


LOGS_DIR_NAME = "logs"
RESULTS_DIR_NAME = "results"

DEFAULT_NEGATIVE_CONTROL_KEYWORDS = ("NTC", "NEC", "NEG", "NEGATIVE", "NEGCTRL", "NEGCONTROL")
DEFAULT_POSITIVE_CONTROL_KEYWORDS = ("POS", "POSITIVE", "PC", "POSCTRL", "POSCONTROL")


class COL:
    """
    Canonical column names used from load_qpcr_file() onward. Raw instrument
    files use their own vocabulary, handled separately by COLUMN_ALIASES.
    """
    SAMPLE = "Sample"
    TARGET = "Target"
    RUN = "Run"
    RUNS = "Runs"
    WELL = "Well"
    FLUOR = "Fluor"
    FLUOROPHORE = "Fluorophore"
    CQ = "Cq"
    CQ_MEAN_RAW = "Cq Mean"
    CQ_SD_RAW = "Cq Std. Dev"
    SQ_RAW = "Starting Quantity (SQ)"
    MEAN_CQ = "Mean Cq"
    CQ_SD = "Cq SD"
    SQ_MEAN = "SQ Mean"
    REPLICATES = "Replicates"
    AMP_STATUS = "AmpStatus_Aggregated"
    AMPLIFICATION_STATUS = "Amplification_Status"
    FLAGS = "Flags"
    CORRECTION_FACTOR = "Correction Factor"
    CORRECTED_SQ_MEAN = "Corrected SQ Mean"
    EXCLUDE_REASON = "Exclude_Reason"
    REPEAT_NEEDED = "Repeat_Needed"
    RUN_PREFIX = "RunPrefix"

    # Internal control correction (Step 3C). BATCH/IC_* are working columns
    # dropped implicitly at export time (CLEAN_SUMMARY_COLUMNS doesn't list
    # them); CORRECTION_FACTOR/CORRECTED_SQ_MEAN above are the two that carry
    # forward into the clean summary.
    BATCH = "Batch"
    IC_CQ = "IC Cq"
    IC_GLOBAL_BASELINE = "IC Global Baseline"
    IC_SOURCE = "IC Source"
    IC_MISSING = "IC_Missing"


class STATUS:
    """Amplification call values, stored in COL.AMP_STATUS."""
    AMP = "Amp"
    NOAMP = "NoAmp"
    INCONCLUSIVE = "Inconclusive"
    NODATA = "NoData"


class REPEAT:
    """Values stored in COL.REPEAT_NEEDED."""
    YES = "YES"
    NO = "NO"


class FLAG:
    """QC flag labels produced by quality_control.collect_flags."""
    LOW_CQ = "low_Cq"
    HIGH_SD = "SD≥2"
    SINGLE_REPLICATE_AMP = "single_replicate_Amp"
    SINGLE_REPLICATE_NOAMP = "single_replicate_NoAmp"
    INCONCLUSIVE = "inconclusive"
    INCONCLUSIVE_LOW_CQ = "inconclusive_Cq<35"
    INCONCLUSIVE_HIGH_CQ = "inconclusive_Cq>35"
    NTC_NEC_AMP = "NTC/NEC_amp"
    POS_CONTROL_INCONCLUSIVE = "POS_control_inconclusive"
    POS_CONTROL_NO_AMP = "POS_control_no_amp"
    MISSING_INTERNAL_CONTROL = "missing_internal_control"


class REASON:
    """Exclusion reasons recorded by quality_control.collapse_interrun."""
    ALL_INCONCLUSIVE = "all_inconclusive_or_no_good_rows"
    MAJORITY_AMP = "majority_Amp_drop_inconclusive"
    MAJORITY_NOAMP = "majority_NoAmp_drop_Inconclusive"
    TIE = "tie_Amp_vs_NoAmp"


class SUFFIX:
    """Output filename suffixes, appended to the working directory's folder name."""
    COMBINED_RESULTS = "_qpcr_combined_results.csv"
    CLEAN_SUMMARY = "_qpcr_clean_summary.xlsx"
    EXCLUSIONS = "_qpcr_exclusions.xlsx"


class SHEET:
    """Sheet names in the exclusions workbook."""
    EXCLUDED_SAMPLES = "Excluded Samples"
    SAMPLES_TO_REPEAT = "Samples to be repeated"
    CONTROLS_EXCLUDED = "Controls (excluded)"


# Sentinels for unmapped/empty values
EMPTY_WELL_SAMPLE = "empty"
UNKNOWN_TARGET = "Unknown"
DEFAULT_RUN_PREFIX = "DEFAULT"
ANY_FLUOROPHORE = "ANY"

# Reason recorded in the exclusions workbook for weak single positives
WEAK_SINGLE_POSITIVE_REASON = "weak_single_amp_Cq>35"

EXCLUSIONS_COLUMNS = (
    COL.SAMPLE, COL.TARGET, COL.RUN, COL.WELL, COL.REPLICATES,
    COL.AMP_STATUS, COL.FLAGS, COL.EXCLUDE_REASON, COL.REPEAT_NEEDED,
)

# Column order for the final "clean summary" export
CLEAN_SUMMARY_COLUMNS = (
    COL.SAMPLE, COL.TARGET,
    COL.AMP_STATUS,
    COL.MEAN_CQ, COL.CQ_SD, COL.SQ_MEAN,
    COL.CORRECTION_FACTOR, COL.CORRECTED_SQ_MEAN,
    COL.REPLICATES, COL.RUNS, COL.WELL,
    COL.FLAGS,
)


@dataclass
class PipelineConfig:
    """Configuration for the qPCR pipeline."""

    # File paths
    workdir: Path
    plate_map_path: Path | None = None
    raw_dir: Path | None = None
    target_mapping_path: Path | None = None
    log_dir: Path | None = None
    results_dir: Path | None = None

    # Quality control thresholds
    pos_cq_strong: float = 35.0  # Cq < 35 considered strong
    low_cq_min: float = 0.1  # Minimum Cq for low Cq flag
    low_cq_max: float = 10.0  # Maximum Cq for low Cq flag
    high_sd: float = 2.0  # High replicate variability threshold

    # Weak single positive exclusion
    exclude_weak_single_positives: bool = True
    weak_single_cq_threshold: float = 35.0  # Cq > this threshold = weak positive

    # Internal control correction (Step 3C). Disabled unless a target name is
    # given (e.g. via the --internal-control CLI flag) — there's no default
    # target, so there's nothing to correct against until the user names one.
    internal_control_target: str | None = None
    # Regex with one capturing group used to pull a "Batch" tag out of a run
    # name, for pairing assay plates with their internal-control plate. Only
    # needed if your run names don't follow "ASSAY_BATCH_DATE" (underscore-
    # delimited, batch = the middle token) — e.g. r"(BATCH\d+)".
    internal_control_batch_regex: str | None = None

    # Bad flags for exclusion
    bad_flags: set[str] | None = None

    # Control detection keywords
    negative_control_keywords: set[str] | None = None
    positive_control_keywords: set[str] | None = None

    # Header detection
    max_scan_lines: int = 200

    # Plate dimensions
    plate_rows: int = 8  # A-H
    plate_cols: int = 12  # 1-12

    def __post_init__(self):
        """Set default paths and mutable defaults based on workdir."""
        if self.plate_map_path is None:
            self.plate_map_path = self.workdir / "plate_maps.csv"
        if self.raw_dir is None:
            self.raw_dir = self.workdir / "raw"
        if self.target_mapping_path is None:
            self.target_mapping_path = self.workdir / "target_mapping.csv"
        if self.log_dir is None:
            self.log_dir = self.workdir / LOGS_DIR_NAME
        if self.results_dir is None:
            self.results_dir = self.workdir / RESULTS_DIR_NAME
        if self.bad_flags is None:
            self.bad_flags = {FLAG.LOW_CQ, FLAG.HIGH_SD}
        if self.negative_control_keywords is None:
            self.negative_control_keywords = set(DEFAULT_NEGATIVE_CONTROL_KEYWORDS)
        if self.positive_control_keywords is None:
            self.positive_control_keywords = set(DEFAULT_POSITIVE_CONTROL_KEYWORDS)


# ---------------------------------------------------------------------------
# Instrument column aliases
#
# Maps each canonical column name to the raw header spellings different
# instrument exports use for it. Drives both header detection (which required
# fields must be present) and the rename step in load_qpcr_file().
# ---------------------------------------------------------------------------
COLUMN_ALIASES = {
    COL.WELL: ['Well', 'Well Position'],
    COL.FLUOR: ['Fluor', 'Reporter'],
    COL.CQ: ['Cq'],
    COL.CQ_MEAN_RAW: ['Cq Mean'],
    COL.CQ_SD_RAW: ['Cq Std. Dev', 'Cq SD'],
    COL.SQ_RAW: ['Starting Quantity (SQ)', 'SQ', 'Quantity', 'Quantity Mean'],
    COL.SAMPLE: ['Sample'],
}

# Column groups that must each have at least one alias present in a candidate
# header row for find_qpcr_header_row() to recognize it. A "Cq" requirement
# is satisfied by either a raw Cq column or a Cq Mean column, since some
# instrument exports only report the latter.
HEADER_DETECTION_GROUPS = {
    'well': COLUMN_ALIASES[COL.WELL],
    'cq': COLUMN_ALIASES[COL.CQ] + COLUMN_ALIASES[COL.CQ_MEAN_RAW],
    'quantity': COLUMN_ALIASES[COL.SQ_RAW],
    'fluor': COLUMN_ALIASES[COL.FLUOR],
}

# Column aliases (lowercased) that indicate an Applied Biosystems export
# rather than a Bio-Rad one.
ABI_INDICATORS = {'well position', 'reporter'}
