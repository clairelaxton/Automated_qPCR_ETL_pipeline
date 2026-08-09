"""
Logging setup for the qPCR pipeline: a console handler at INFO level plus a
timestamped file handler at DEBUG level, both attached to the root logger so
every module's logger.info()/warning()/error() calls are captured.
"""
import logging
from pathlib import Path
from datetime import datetime

from .config import LOGS_DIR_NAME


def setup_logger(workdir: Path, log_dir: Path | None = None) -> logging.Logger:
    """
    Configure the root logger with console + file handlers and return it.

    Args:
        workdir: Working directory
        log_dir: Directory for the log file (defaults to workdir/logs)

    Returns:
        The configured root logger instance
    """
    if log_dir is None:
        log_dir = workdir / LOGS_DIR_NAME

    log_dir.mkdir(parents=True, exist_ok=True)

    # Create timestamped log file
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M")
    log_file = log_dir / f"qpcr_pipeline_{timestamp}.log"

    # Configure root logger so all module logs are captured.
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Clear any existing handlers
    logger.handlers.clear()

    # Console handler (simple output)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_format = logging.Formatter("%(levelname)s: %(message)s")
    console_handler.setFormatter(console_format)

    # File handler (detailed log with timestamps)
    file_handler = logging.FileHandler(log_file, encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_format = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    file_handler.setFormatter(file_format)

    # Attach handlers
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    logger.info(f"🪵 Logging initialized. Detailed log will be saved to: {log_file}")

    return logger
