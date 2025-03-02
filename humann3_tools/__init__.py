# humann3_tools/__init__.py
"""
HUMAnN3 Tools - A package for processing and analyzing HUMAnN3 output data.

This package provides functions for:
1. Processing HUMAnN3 output files (normalization, joining, splitting)
2. Downstream analysis (statistical tests, visualization)
3. Differential abundance analysis (ALDEx2, ANCOM, ANCOM-BC)
"""

__version__ = "0.1.0"

# Import main functions for easy access
from humann3_tools.main import (
    run_full_pipeline,
    process_humann3_files_only,
    analyze_existing_humann3_files,
    run_pathway_differential_abundance,
    run_gene_differential_abundance,
    run_preprocessing_and_analysis  
)

# Make logger functions available at the package level
from humann3_tools.logger import setup_logger, log_print

# Import key utility functions
from humann3_tools.utils.file_utils import check_file_exists, sanitize_filename

# Import differential abundance functions
from humann3_tools.analysis.differential_abundance import (
    aldex2_like,
    ancom,
    ancom_bc,
    run_differential_abundance_analysis
)

# Import preprocessing functions for easier access
from humann3_tools.preprocessing.kneaddata import (
    run_kneaddata,
    run_kneaddata_parallel
)

from humann3_tools.preprocessing.pipeline import (
    run_preprocessing_pipeline,
    run_preprocessing_pipeline_parallel
)