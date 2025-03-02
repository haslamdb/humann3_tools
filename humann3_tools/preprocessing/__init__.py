# humann3_tools/preprocessing/__init__.py
"""
Preprocessing module for HUMAnN3 Tools.

This module provides functions to preprocess raw sequence files using KneadData and run HUMAnN3 on the results.
"""

from humann3_tools.preprocessing.kneaddata import run_kneaddata, check_kneaddata_installation
from humann3_tools.preprocessing.humann3_run import run_humann3, check_humann3_installation
from humann3_tools.preprocessing.pipeline import run_preprocessing_pipeline