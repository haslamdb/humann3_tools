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
    suffixes = [
        ".paired_Abundance-CPM", 
        "_Abundance-CPM",
        ".paired_Abundance-RELAB", 
        "_Abundance-RELAB",
        "-cpm",
        "-relab"
    ]
    for suffix in suffixes:
        if col.lower().endswith(suffix.lower()):
            return col[:-len(suffix)]
    return col

def strip_suffixes_from_file_headers(file_path, logger=None):
    """
    Remove HUMAnN3 abundance suffixes from column headers in a file.
    
    Args:
        file_path: Path to the file with headers to strip
        logger: Optional logger for messages
        
    Returns:
        True if successful, False otherwise
    """
    if logger is None:
        logger = logging.getLogger('humann3_analysis')
    
    try:
        # Read the file
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        if not lines:
            logger.warning(f"Empty file: {file_path}")
            return False
        
        # Find the first non-comment line (header line)
        header_index = 0
        for i, line in enumerate(lines):
            if not line.startswith('#'):
                header_index = i
                break
        
        # Get the header line and split it
        header = lines[header_index].strip()
        cols = header.split('\t')
        
        # Apply strip_suffix to each column except the first one (which is usually the feature ID)
        cols = [cols[0]] + [strip_suffix(col) for col in cols[1:]]
        
        # Update the header line
        lines[header_index] = '\t'.join(cols) + '\n'
        
        # Write the file back
        with open(file_path, 'w') as f:
            f.writelines(lines)
        
        logger.info(f"Stripped suffixes from headers in: {file_path}")
        return True
    except Exception as e:
        if logger:
            logger.error(f"Error stripping suffixes from headers in {file_path}: {str(e)}")
        return False