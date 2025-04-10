# humann3_tools/utils/__init__.py
"""
Utility modules for humann3_tools.

This package contains various utility modules:
- cmd_utils.py: Utilities for running shell commands
- file_utils.py: Utilities for file handling
- input_handler.py: Utilities for handling different input methods
- resource_utils.py: Utilities for monitoring and managing system resources
"""

# Import key functions
from humann3_tools.utils.file_utils import (
    check_file_exists,
    strip_suffixes_from_file_headers,
    sanitize_filename
)

from humann3_tools.utils.cmd_utils import run_cmd

from humann3_tools.utils.input_handler import (
    get_input_files,
    find_sample_files,
    collect_files_from_metadata,
    read_samples_file
)

from humann3_tools.utils.resource_utils import (
    track_peak_memory,
    limit_memory_usage
)
