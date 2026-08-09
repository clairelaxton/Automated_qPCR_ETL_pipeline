"""
Top-level orchestration for the qPCR pipeline: wires together plate-map
loading, per-file processing, summarization, target mapping, QC flagging,
inter-run collapsing, and export into the sequence QPCRPipeline.run() runs.
"""
import sys
import pandas as pd
from pathlib import Path
import logging

from .config import (
    PipelineConfig, EXCLUSIONS_COLUMNS, CLEAN_SUMMARY_COLUMNS,
    WEAK_SINGLE_POSITIVE_REASON, COL, STATUS, REPEAT, SUFFIX, SHEET,
    LOGS_DIR_NAME, RESULTS_DIR_NAME,
)
from .logger_setup import setup_logger
from .file_io import load_plate_map
from .processor import process_all_qpcr_files
from .statistics import create_summary_table
from .target_mapping import apply_target_mapping, closest_match
from .quality_control import flag_samples, collapse_interrun
from .internal_control import apply_internal_control_correction


logger = logging.getLogger(__name__)


class QPCRPipeline:
    """Runs the full qPCR analysis pipeline for one working directory."""

    def __init__(self, workdir: Path | None = None, config: PipelineConfig | None = None):
        """
        Set up logging, validate inputs, and compute output file paths.

        Args:
            workdir: Working directory (defaults to current directory)
            config: Pipeline configuration (optional; a default is built from workdir)
        """
        used_default_workdir = workdir is None
        if workdir is None:
            workdir = Path.cwd()
        else:
            workdir = Path(workdir).resolve()

        if config is None:
            config = PipelineConfig(workdir=workdir)
        elif config.workdir != workdir:
            self._retarget_workdir(config, workdir)

        self.config = config
        self.logger = setup_logger(workdir, config.log_dir)
        if used_default_workdir:
            self.logger.info(f"No path provided. Using current working directory: {workdir}")
        else:
            self.logger.info(f"Using provided working directory: {workdir}")

        self._validate_inputs()

        # Output paths
        self.folder_name = workdir.name
        self.results_dir = self.config.results_dir
        self.results_dir.mkdir(exist_ok=True)
        self.logger.info(f"Folder name identified as: {self.folder_name}")
        self.combined_csv = workdir / f"{self.folder_name}{SUFFIX.COMBINED_RESULTS}"
        self.clean_summary_xlsx = self.results_dir / f"{self.folder_name}{SUFFIX.CLEAN_SUMMARY}"
        self.exclusions_xlsx = self.results_dir / f"{self.folder_name}{SUFFIX.EXCLUSIONS}"
        self.logger.info("Output files will be saved to:")
        self.logger.info(f"  Combined results CSV: {self.combined_csv}")
        self.logger.info(f"  Results folder:       {self.results_dir}")
        self.logger.info(f"  Clean summary Excel:  {self.clean_summary_xlsx}")
        self.logger.info(f"  Exclusions Excel:     {self.exclusions_xlsx}")

    @staticmethod
    def _retarget_workdir(config: PipelineConfig, new_workdir: Path) -> None:
        """
        Point an already-built config at new_workdir. Any path field that was
        auto-derived from the config's old workdir (i.e. never explicitly
        customized to somewhere unrelated) is re-derived against the new one
        too — otherwise only reassigning config.workdir would silently leave
        raw_dir/plate_map_path/etc. pointed at the old directory.
        """
        old_workdir = config.workdir
        derived = {
            "plate_map_path": old_workdir / "plate_maps.csv",
            "raw_dir": old_workdir / "raw",
            "target_mapping_path": old_workdir / "target_mapping.csv",
            "log_dir": old_workdir / LOGS_DIR_NAME,
            "results_dir": old_workdir / RESULTS_DIR_NAME,
        }
        for field, old_default in derived.items():
            if getattr(config, field) == old_default:
                setattr(config, field, new_workdir / old_default.relative_to(old_workdir))
        config.workdir = new_workdir

    def _validate_inputs(self) -> None:
        """Exit with an error if plate_maps.csv or the raw/ folder (with CSVs) is missing."""
        if not self.config.plate_map_path.exists() or not self.config.plate_map_path.is_file():
            self.logger.error(f"Missing required file: {self.config.plate_map_path}")
            self.logger.error("Please ensure 'plate_maps.csv' is in your working directory.")
            sys.exit(1)

        if not self.config.raw_dir.exists() or not self.config.raw_dir.is_dir():
            self.logger.error(f"Missing required folder: {self.config.raw_dir}")
            self.logger.error("Please ensure there is a folder named 'raw' containing .csv files.")
            sys.exit(1)

        raw_csv_files = list(self.config.raw_dir.glob("*.csv"))
        if not raw_csv_files:
            self.logger.error("⚠️ No CSV files found in the 'raw' folder.")
            sys.exit(1)

        self.logger.info("✅ Found plate_maps.csv file.")
        self.logger.info("✅ Found 'raw' data folder.")
        self.logger.info(f"Found {len(raw_csv_files)} CSV file(s) in 'raw' folder.")

    def run(self) -> None:
        """Run every pipeline step in sequence and write all output files."""
        # Step 1: Load plate map
        plate_map = load_plate_map(self.config.plate_map_path, self.config)

        # Step 2: Process qPCR files
        qpcr_files = sorted(self.config.raw_dir.glob("*.csv"))

        # Check for missing runs
        qpcr_file_basenames = [p.stem for p in qpcr_files]
        unique_runs = pd.unique(plate_map[COL.RUN])
        missing_runs = [run for run in unique_runs if run not in qpcr_file_basenames]
        if missing_runs:
            self.logger.warning(
                "⚠️ Missing qPCR data files for runs referenced in plate_maps.csv "
                "— that plate map's samples will be skipped entirely. "
                f"Available raw/ files: {', '.join(qpcr_file_basenames) or '(none)'}"
            )
            for run in missing_runs:
                suggestion = closest_match(run, qpcr_file_basenames)
                hint = f" — did you mean '{suggestion}'?" if suggestion else ""
                self.logger.warning(f"  - '{run}'{hint}")

        combined, skipped_files = process_all_qpcr_files(
            qpcr_files,
            plate_map,
            self.config
        )

        if combined.empty:
            self.logger.error("⚠️ No valid qPCR data was processed. Exiting.")
            sys.exit(1)

        # Save combined results
        combined.to_csv(self.combined_csv, index=False, encoding='utf-8')
        self.logger.info(f"✅ Combined results saved to: {self.combined_csv}")

        # Step 3: Create summary
        summary = create_summary_table(combined)

        # Step 3B: Apply target mapping
        summary = apply_target_mapping(summary, self.config)

        # Step 3C: Optional internal control correction
        try:
            summary = apply_internal_control_correction(summary, self.config)
        except ValueError as e:
            self.logger.error(f"❌ {e}")
            sys.exit(1)

        # Step 3D: Flag samples
        summary = flag_samples(summary, self.config)

        # Step 4: Collapse inter-run duplicates
        final_summary, excluded_df, controls_excluded = collapse_interrun(
            summary,
            self.config
        )

        # Step 5: Clean and export
        self._export_final_outputs(final_summary, excluded_df, controls_excluded)

    def _exclude_weak_single_positives(
        self,
        final_summary: pd.DataFrame,
        excluded_df: pd.DataFrame,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Pull weak single-replicate positives (Amp, 1 replicate, high Mean Cq) out
        of the final summary and record them as exclusions needing a repeat.
        """
        if not self.config.exclude_weak_single_positives:
            return final_summary, excluded_df

        weak_mask = (
            (final_summary[COL.AMP_STATUS].astype(str) == STATUS.AMP) &
            (pd.to_numeric(final_summary[COL.REPLICATES], errors="coerce").fillna(0).astype(int) == 1) &
            (pd.to_numeric(final_summary[COL.MEAN_CQ], errors="coerce") > self.config.weak_single_cq_threshold)
        )
        weak_singles = final_summary.loc[weak_mask].copy()

        if weak_singles.empty:
            self.logger.info("🧹 No weak single positives detected to exclude.")
            return final_summary, excluded_df

        self.logger.info(
            f"🧹 Excluding {len(weak_singles)} weak single positives "
            f"(Amp, 1 replicate, Mean Cq > {self.config.weak_single_cq_threshold})."
        )

        weak_log = weak_singles.copy()
        weak_log[COL.EXCLUDE_REASON] = WEAK_SINGLE_POSITIVE_REASON
        weak_log[COL.REPEAT_NEEDED] = REPEAT.YES
        if COL.RUN not in weak_log.columns and COL.RUNS in weak_log.columns:
            weak_log = weak_log.rename(columns={COL.RUNS: COL.RUN})
        weak_log = weak_log.reindex(columns=list(EXCLUSIONS_COLUMNS))

        excluded_df = pd.concat([excluded_df, weak_log], ignore_index=True)
        final_summary = final_summary.loc[~weak_mask].copy()
        return final_summary, excluded_df

    def _build_clean_summary(self, final_summary: pd.DataFrame) -> pd.DataFrame:
        """Rename Run to Runs, restrict/order columns to CLEAN_SUMMARY_COLUMNS, and sort by Sample/Target."""
        final_summary = final_summary.rename(columns={COL.RUN: COL.RUNS})

        present_cols = [c for c in CLEAN_SUMMARY_COLUMNS if c in final_summary.columns]
        clean_summary = final_summary[present_cols].copy()

        sort_cols = [c for c in [COL.SAMPLE, COL.TARGET] if c in clean_summary.columns]
        if sort_cols:
            clean_summary = clean_summary.sort_values(by=sort_cols).reset_index(drop=True)
        return clean_summary

    def _build_repeat_summary(self, excluded_df: pd.DataFrame) -> pd.DataFrame:
        """Group excluded rows marked Repeat_Needed=YES by Sample into a Sample/Targets-Needed table."""
        columns = [COL.SAMPLE, "Targets Needed"]
        if excluded_df is None or excluded_df.empty:
            return pd.DataFrame(columns=columns)

        needs_repeat = excluded_df[COL.REPEAT_NEEDED].astype(str).str.strip().str.upper().eq(REPEAT.YES)
        repeat_rows = excluded_df.loc[needs_repeat, [COL.SAMPLE, COL.TARGET]].dropna()
        if repeat_rows.empty:
            return pd.DataFrame(columns=columns)

        return (
            repeat_rows
            .groupby(COL.SAMPLE, as_index=False)[COL.TARGET]
            .agg(lambda s: ", ".join(sorted({t.strip() for t in map(str, s) if t and t.strip()})))
            .rename(columns={COL.TARGET: "Targets Needed"})
            .sort_values(by=COL.SAMPLE)
            .reset_index(drop=True)
        )

    def _normalize_exclusions_table(self, excluded_df: pd.DataFrame) -> pd.DataFrame:
        """Normalize Repeat_Needed casing, restrict/order columns to EXCLUSIONS_COLUMNS, and sort."""
        if excluded_df is None or excluded_df.empty:
            return pd.DataFrame(columns=list(EXCLUSIONS_COLUMNS))

        excl_all = excluded_df.copy()
        excl_all[COL.REPEAT_NEEDED] = excl_all[COL.REPEAT_NEEDED].astype(str).str.strip().str.upper()
        excl_all = excl_all.reindex(columns=[c for c in EXCLUSIONS_COLUMNS if c in excl_all.columns])

        sort_cols = [c for c in [COL.REPEAT_NEEDED, COL.SAMPLE, COL.TARGET] if c in excl_all.columns]
        if sort_cols:
            excl_all = excl_all.sort_values(by=sort_cols).reset_index(drop=True)
        return excl_all

    def _export_final_outputs(
        self,
        final_summary: pd.DataFrame,
        excluded_df: pd.DataFrame,
        controls_excluded: pd.DataFrame
    ) -> None:
        """Exclude weak single positives, then write the clean summary and the exclusions workbook."""
        if final_summary is None or final_summary.empty:
            self.logger.error("STEP 5: final_summary is empty or missing.")
            sys.exit(1)

        if excluded_df is None:
            excluded_df = pd.DataFrame(columns=list(EXCLUSIONS_COLUMNS))

        final_summary, excluded_df = self._exclude_weak_single_positives(final_summary, excluded_df)

        final_clean_summary = self._build_clean_summary(final_summary)
        final_clean_summary.to_excel(self.clean_summary_xlsx, index=False)
        self.logger.info(f"✅ Clean summary saved to: {self.clean_summary_xlsx}")

        repeat_summary_by_sample = self._build_repeat_summary(excluded_df)
        excl_all = self._normalize_exclusions_table(excluded_df)

        with pd.ExcelWriter(self.exclusions_xlsx, engine="openpyxl") as xw:
            excl_all.to_excel(xw, sheet_name=SHEET.EXCLUDED_SAMPLES, index=False)
            repeat_summary_by_sample.to_excel(xw, sheet_name=SHEET.SAMPLES_TO_REPEAT, index=False)
            if not controls_excluded.empty:
                controls_excluded.to_excel(xw, sheet_name=SHEET.CONTROLS_EXCLUDED, index=False)

        self.logger.info(f"📄 Exclusions workbook saved to: {self.exclusions_xlsx}")
