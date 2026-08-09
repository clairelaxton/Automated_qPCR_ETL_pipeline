"""
qPCR Analysis Pipeline

A comprehensive pipeline for processing and analyzing qPCR data from
multiple instrument platforms (Bio-Rad CFX and Applied Biosystems).
"""

__version__ = "1.0.0"

from .pipeline import QPCRPipeline

__all__ = ["QPCRPipeline"]
