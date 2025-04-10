# humann3_tools/humann3_tools/utils/__init__.py
"""Utility functions for humann3_tools."""

from humann3_tools.utils.file_utils import (
    check_file_exists,
    check_file_exists_with_logger,
    sanitize_filename,
    strip_suffix
)

from humann3_tools.utils.cmd_utils import run_cmd

from humann3_tools.utils.sample_utils import (
    validate_sample_key,
    validate_sample_key_noninteractive,
    check_input_files_exist
)

from humann3_tools.utils.resource_utils import (
    calculate_optimal_resources,
    get_memory_usage,
    log_resource_usage,
    estimate_memory_requirements,
    check_resource_availability,
    monitor_memory_usage,
    stop_memory_monitoring,
    track_peak_memory,
    limit_memory_usage
)

from humann3_tools.utils.metadata_utils import (
    collect_samples_from_metadata,
    find_sample_files,
    prompt_for_sequence_file_patterns,
    read_samples_file
)
