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
    """
    Enhanced function to remove HUMAnN3 abundance suffix from column names.
    Handles various suffix patterns more robustly.
    """
    # First, try exact matches with known suffixes
    suffixes = [
        ".paired_Abundance-CPM", 
        "_Abundance-CPM",
        ".paired_Abundance-RELAB", 
        "_Abundance-RELAB",
        "-cpm",
        "-relab",
        ".cpm",
        ".relab"
    ]
    
    for suffix in suffixes:
        if col.lower().endswith(suffix.lower()):
            return col[:-len(suffix)]
    
    # If no exact match, try pattern-based detection
    
    # Pattern 1: Sample.UNIT or Sample.something.UNIT
    parts = col.split('.')
    if len(parts) > 1:
        last_part = parts[-1].lower()
        if any(unit in last_part for unit in ['cpm', 'relab', 'abundance']):
            return '.'.join(parts[:-1])
    
    # Pattern 2: Sample_UNIT or Sample-UNIT or variations
    for marker in ['_abundance', '-abundance', '_cpm', '-cpm', '_relab', '-relab']:
        pos = col.lower().find(marker)
        if pos > 0:
            return col[:pos]
    
    # No match found
    return col

def strip_suffixes_from_file_headers(file_path, logger=None):
    """
    Remove HUMAnN3 abundance suffixes from column headers in a file.
    Uses enhanced suffix detection with multiple patterns.
    
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
        new_cols = [cols[0]]
        change_count = 0
        
        for col in cols[1:]:
            new_col = strip_suffix(col)
            new_cols.append(new_col)
            if new_col != col:
                change_count += 1
                logger.debug(f"Stripped suffix: '{col}' -> '{new_col}'")
        
        # Check if any columns were actually modified
        if change_count == 0:
            logger.info(f"No suffix patterns found in headers of: {file_path}")
            return True
        
        # Update the header line
        lines[header_index] = '\t'.join(new_cols) + '\n'
        
        # Write the file back
        with open(file_path, 'w') as f:
            f.writelines(lines)
        
        logger.info(f"Stripped {change_count} suffixes from headers in: {file_path}")
        return True
    except Exception as e:
        if logger:
            logger.error(f"Error stripping suffixes from headers in {file_path}: {str(e)}")
        return False