# humann3_tools/core/__init__.py
"""
Core functionality for HUMAnN3 Tools.

This package contains the core functionality for the different steps of the workflow:
- kneaddata.py: KneadData processing
- humann3.py: HUMAnN3 processing  
- join_unstratify.py: Join and unstratify operations
"""

from humann3_tools.core.kneaddata import (
    check_kneaddata_installation,
    run_kneaddata,
    run_kneaddata_parallel
)

from humann3_tools.preprocessing.humann3_run import (
    check_humann3_installation,
    run_humann3,
    run_humann3_parallel
)

from humann3_tools.core.join_unstratify import (
    process_join_unstratify
)
