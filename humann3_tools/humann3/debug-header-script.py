#!/usr/bin/env python3
"""
Debug script to examine column headers in HUMAnN3 output files and test suffix removal.

Usage:
    python debug_headers.py <file_path>

Example:
    python debug_headers.py /path/to/pathway_abundance-cpm_unstratified.tsv
"""

import sys
import os

def examine_file_headers(file_path):
    """
    Examine the headers in a file to help debug suffix stripping issues.
    """
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        if not lines:
            print(f"Empty file: {file_path}")
            return
        
        # Find the first non-comment line (header line)
        header_index = 0
        for i, line in enumerate(lines):
            if not line.startswith('#'):
                header_index = i
                break
        
        # Get the header line and split it
        header = lines[header_index].strip()
        cols = header.split('\t')
        
        print(f"File: {file_path}")
        print(f"Header line index: {header_index}")
        print(f"Number of columns: {len(cols)}")
        print(f"First column: '{cols[0]}'")
        
        # Print a sample of column headers
        sample_size = min(5, len(cols) - 1)
        print(f"\nSample of {sample_size} column headers:")
        for i in range(1, 1 + sample_size):
            print(f"  Column {i}: '{cols[i]}'")
        
        # Test different suffix patterns
        potential_suffixes = [
            ".paired_Abundance-CPM", 
            "_Abundance-CPM",
            ".paired_Abundance-RELAB", 
            "_Abundance-RELAB",
            "-cpm",
            "-relab",
            ".cpm",
            ".relab"
        ]
        
        print("\nTesting suffix patterns on the first few columns:")
        for i in range(1, 1 + sample_size):
            col = cols[i]
            print(f"\n  Testing column {i}: '{col}'")
            
            for suffix in potential_suffixes:
                lower_col = col.lower()
                lower_suffix = suffix.lower()
                ends_with_suffix = lower_col.endswith(lower_suffix)
                
                print(f"    Suffix '{suffix}': {'MATCH' if ends_with_suffix else 'no match'}")
                
                if ends_with_suffix:
                    without_suffix = col[:-len(suffix)]
                    print(f"      Would become: '{without_suffix}'")

        # Test enhanced suffix removal
        print("\nTesting enhanced suffix removal logic:")
        for i in range(1, 1 + sample_size):
            col = cols[i]
            print(f"\n  Column {i}: '{col}'")
            
            # Look for known sample/group pattern with unit suffix
            parts = col.split('.')
            if len(parts) > 1:
                print(f"    Contains period(s). Parts: {parts}")
                
                # Check if the last part is a unit suffix
                last_part = parts[-1].lower()
                if any(unit in last_part for unit in ['cpm', 'relab', 'abundance']):
                    print(f"    Last part '{last_part}' appears to be a unit suffix")
                    new_col = '.'.join(parts[:-1])
                    print(f"    After cleanup: '{new_col}'")
                else:
                    print(f"    Last part '{last_part}' doesn't appear to be a unit suffix")
            
            # If no period, check for underscore or dash patterns
            elif '_' in col or '-' in col:
                print(f"    Contains underscore(s) or dash(es)")
                
                # Attempt to identify the position where the suffix starts
                # Common patterns: SampleName_cpm, SampleName-Abundance-CPM, etc.
                for marker in ['_abundance', '-abundance', '_cpm', '-cpm', '_relab', '-relab']:
                    pos = col.lower().find(marker)
                    if pos > 0:
                        print(f"    Found marker '{marker}' at position {pos}")
                        new_col = col[:pos]
                        print(f"    After cleanup: '{new_col}'")
                        break
                else:
                    print(f"    No common unit marker found")
            
            else:
                print(f"    No periods, underscores, or dashes found - likely no suffix")
                
        return cols
    except Exception as e:
        print(f"Error examining file: {str(e)}")
        return None

def strip_suffix_improved(col):
    """
    Enhanced version of strip_suffix with better pattern matching.
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
            new_col = col[:-len(suffix)]
            print(f"Exact match! '{col}' -> '{new_col}' (removed '{suffix}')")
            return new_col
    
    # If no exact match, try pattern-based detection
    
    # Pattern 1: Sample.UNIT or Sample.something.UNIT
    parts = col.split('.')
    if len(parts) > 1:
        last_part = parts[-1].lower()
        if any(unit in last_part for unit in ['cpm', 'relab', 'abundance']):
            new_col = '.'.join(parts[:-1])
            print(f"Pattern match! '{col}' -> '{new_col}' (removed '.{parts[-1]}')")
            return new_col
    
    # Pattern 2: Sample_UNIT or Sample-UNIT or variations
    for marker in ['_abundance', '-abundance', '_cpm', '-cpm', '_relab', '-relab']:
        pos = col.lower().find(marker)
        if pos > 0:
            new_col = col[:pos]
            print(f"Pattern match! '{col}' -> '{new_col}' (removed '{col[pos:]}')")
            return new_col
    
    # No match found
    print(f"No suffix pattern found in '{col}'")
    return col

def test_strip_suffix_on_file(file_path):
    """
    Test the improved strip_suffix function on a real file.
    """
    print(f"\nTesting improved strip_suffix on file: {file_path}")
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Find the header line
        header_index = 0
        for i, line in enumerate(lines):
            if not line.startswith('#'):
                header_index = i
                break
        
        # Get and split the header
        header = lines[header_index].strip()
        cols = header.split('\t')
        
        # Apply the improved strip_suffix
        new_cols = [cols[0]]  # Keep first column unchanged
        change_count = 0
        
        for col in cols[1:]:
            new_col = strip_suffix_improved(col)
            new_cols.append(new_col)
            if new_col != col:
                change_count += 1
        
        print(f"\nTotal columns processed: {len(cols) - 1}")
        print(f"Columns changed: {change_count}")
        
        # Ask if user wants to save the changes
        if change_count > 0:
            choice = input("\nDo you want to save these changes to the file? (y/n): ")
            if choice.lower() == 'y':
                # Update the header line
                lines[header_index] = '\t'.join(new_cols) + '\n'
                
                # Write back to the file
                with open(file_path, 'w') as f:
                    f.writelines(lines)
                print(f"Changes saved to {file_path}")
            else:
                print("Changes not saved")
    
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <file_path>")
        sys.exit(1)
    
    file_path = sys.argv[1]
    if not os.path.exists(file_path):
        print(f"Error: File does not exist: {file_path}")
        sys.exit(1)
    
    examine_file_headers(file_path)
    test_strip_suffix_on_file(file_path)
