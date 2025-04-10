# src/humann3_tools/preprocessing/__init__.py
"""
Preprocessing module for HUMAnN3 Tools.

This module provides functions to preprocess raw sequence files using KneadData and run HUMAnN3 on the results.
"""

# Import from core to maintain backward compatibility
from humann3_tools.core.kneaddata import (
    run_kneaddata, 
    run_kneaddata_parallel, 
    check_kneaddata_installation
)

from humann3_tools.core.humann3 import (
    run_humann3, 
    run_humann3_parallel, 
    check_humann3_installation
)

# Import pipeline functions if available
try:
    from humann3_tools.preprocessing.pipeline import (
        run_preprocessing_pipeline, 
        run_preprocessing_pipeline_parallel
    )
except ImportError:
    pass
