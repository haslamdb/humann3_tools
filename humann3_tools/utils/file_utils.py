# humann3_tools/humann3_tools/utils/file_utils.py
import os
import logging

def check_file_exists(filepath, description):
    """Check if a file exists and is readable."""
    logger = logging.getLogger('humann3_analysis')
    
    if not os.path.isfile(filepath):
        logger.error(f"ERROR: {description} file does not exist: {filepath}")
        return False
    
    if not os.access(filepath, os.R_OK):
        logger.error(f"ERROR: {description} file exists but is not readable: {filepath}")
        return False
    
    return True

def check_file_exists_with_logger(filepath, description, logger):
    """
    Check if file exists & is readable (logger-based).
    """
    if not os.path.isfile(filepath):
        logger.error(f"{description} file does not exist: {filepath}")
        return False
    if not os.access(filepath, os.R_OK):
        logger.error(f"{description} file is not readable: {filepath}")
        return False
    return True

def sanitize_filename(filename):
    """Replace invalid filename characters with underscores."""
    import re
    return re.sub(r'[<>:"/\\|?*]', '_', filename)

def strip_suffix(col):
    """Remove HUMAnN3 abundance suffix (if present) from column names."""
    return col.replace(".paired_Abundance-CPM", "")