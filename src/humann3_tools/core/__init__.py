# humann3_tools/core/__init__.py
"""
Core functionality for HUMAnN3 Tools.

This package contains the core functionality for the different steps of the workflow:
- kneaddata.py: KneadData processing
- humann3.py: HUMAnN3 processing  
- join_unstratify.py: Join and unstratify operations
"""

from src.humann3_tools.core.kneaddata import (
    check_kneaddata_installation,
    run_kneaddata,
    run_kneaddata_parallel
)

from src.humann3_tools.core.humann3 import (
    run_full_pipeline,
    run_preprocessing_and_analysis,
    process_humann3_files_only,
    analyze_existing_humann3_files,
    run_pathway_differential_abundance,
    run_gene_differential_abundance
)

from src.humann3_tools.core.join_unstratify import (
    process_join_unstratify
)
