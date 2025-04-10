# humann3_tools/__init__.py
"""
HUMAnN3 Tools - A comprehensive toolkit for metagenomic analysis with HUMAnN3.

This package provides a modular workflow for processing and analyzing metagenomic data:
1. Quality control and host depletion with KneadData
2. Functional profiling with HUMAnN3
3. Joining, normalizing, and unstratifying HUMAnN3 output files
4. Statistical testing
5. Differential abundance analysis
6. Visualization
"""

__version__ = "0.1.0"

# Import key functions for easy access
from humann3_tools.logger import setup_logger, log_print

# Import main functions
from humann3_tools.core.humann3 import (
    run_full_pipeline,
    process_humann3_files_only,
    analyze_existing_humann3_files,
    run_pathway_differential_abundance,
    run_gene_differential_abundance,
    run_preprocessing_and_analysis  
)
