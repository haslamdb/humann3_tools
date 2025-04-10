#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import glob
import logging
from typing import List, Dict, Tuple, Optional, Union

# Set up logging
logger = logging.getLogger('humann3_workflow')
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)

def find_sample_files(sample_id: str, 
                      search_dir: str, 
                      file_pattern: Optional[str] = None,
                      r1_suffix: Optional[str] = None, 
                      r2_suffix: Optional[str] = None,
                      paired: bool = False) -> List[str]:
    """
    Find sequence files for a sample based on patterns.
    
    Args:
        sample_id: Sample identifier
        search_dir: Directory to search for sequence files
        file_pattern: Optional pattern to use for file search
        r1_suffix: Suffix for R1 files in paired-end mode 
        r2_suffix: Suffix for R2 files in paired-end mode
        paired: Whether to look for paired-end files
        
    Returns:
        List of file paths matching the sample
    """
    files = []
    
    # If a specific file pattern is provided
    if file_pattern:
        pattern = file_pattern.replace('{sample}', sample_id)
        search_pattern = os.path.join(search_dir, pattern)
        files = sorted(glob.glob(search_pattern))
        
    # For paired-end reads
    elif paired and r1_suffix and r2_suffix:
        r1_pattern = f"{sample_id}{r1_suffix}"
        r2_pattern = f"{sample_id}{r2_suffix}"
        r1_search = os.path.join(search_dir, r1_pattern)
        r2_search = os.path.join(search_dir, r2_pattern)
        
        r1_files = sorted(glob.glob(r1_search))
        r2_files = sorted(glob.glob(r2_search))
        
        # Ensure we have matching pairs
        if len(r1_files) == len(r2_files) and len(r1_files) > 0:
            # Return interleaved list [r1_1, r2_1, r1_2, r2_2, ...]
            for r1, r2 in zip(r1_files, r2_files):
                files.extend([r1, r2])
        else:
            logger.warning(f"Couldn't find matching paired files for sample {sample_id}")
    
    # For single-end reads
    else:
        # Try common patterns
        common_patterns = [
            f"{sample_id}.fastq", 
            f"{sample_id}.fq", 
            f"{sample_id}.fastq.gz", 
            f"{sample_id}.fq.gz",
            f"{sample_id}_*.fastq.gz",
            f"{sample_id}*.fastq.gz"
        ]
        
        for pattern in common_patterns:
            search_pattern = os.path.join(search_dir, pattern)
            matches = glob.glob(search_pattern)
            if matches:
                files.extend(matches)
                break
    
    return files

def collect_samples_from_metadata(metadata_file: str, 
                                 seq_dir: str, 
                                 sample_col: Optional[str] = None,
                                 r1_col: Optional[str] = None, 
                                 r2_col: Optional[str] = None,
                                 file_pattern: Optional[str] = None,
                                 r1_suffix: Optional[str] = None, 
                                 r2_suffix: Optional[str] = None,
                                 paired: bool = False) -> Dict[str, List[str]]:
    """
    Collect sample information from metadata file and find associated sequence files.
    
    Args:
        metadata_file: Path to metadata CSV file
        seq_dir: Directory containing sequence files
        sample_col: Column name for sample identifiers
        r1_col: Column name for R1 file paths
        r2_col: Column name for R2 file paths
        file_pattern: Pattern for finding files (e.g., "{sample}_S*_R*.fastq.gz")
        r1_suffix: Suffix for R1 files (e.g., "_R1.fastq.gz")
        r2_suffix: Suffix for R2 files (e.g., "_R2.fastq.gz")
        paired: Whether samples are paired-end
        
    Returns:
        Dictionary mapping sample IDs to sequence file paths
    """
    # Read metadata file
    try:
        metadata = pd.read_csv(metadata_file)
    except Exception as e:
        logger.error(f"Error reading metadata file: {e}")
        return {}
    
    # Identify sample ID column
    if not sample_col:
        # Try to auto-detect
        common_sample_cols = ['SampleID', 'Sample_ID', 'SampleName', 'sample_id', 'sample_name', 'Sample', 'ID']
        for col in common_sample_cols:
            if col in metadata.columns:
                sample_col = col
                logger.info(f"Auto-detected sample ID column: {sample_col}")
                break
        
        if not sample_col:
            sample_col = metadata.columns[0]
            logger.warning(f"Could not auto-detect sample ID column, using the first column: {sample_col}")
    
    # Check if there are explicit file path columns
    has_file_cols = r1_col in metadata.columns
    if paired:
        has_file_cols = has_file_cols and r2_col in metadata.columns
    
    samples = {}
    
    # Process each row in the metadata
    for _, row in metadata.iterrows():
        sample_id = str(row[sample_col])
        
        # Skip rows with missing sample IDs
        if not sample_id or pd.isna(sample_id):
            continue
        
        # Case 1: File paths are in the metadata
        if has_file_cols:
            if paired:
                r1_path = row[r1_col]
                r2_path = row[r2_col]
                
                # Check if paths are relative or absolute
                if not os.path.isabs(r1_path):
                    r1_path = os.path.join(seq_dir, r1_path)
                if not os.path.isabs(r2_path):
                    r2_path = os.path.join(seq_dir, r2_path)
                
                samples[sample_id] = [r1_path, r2_path]
            else:
                file_path = row[r1_col]
                if not os.path.isabs(file_path):
                    file_path = os.path.join(seq_dir, file_path)
                samples[sample_id] = [file_path]
        
        # Case 2: Need to find files based on patterns
        else:
            files = find_sample_files(
                sample_id=sample_id,
                search_dir=seq_dir,
                file_pattern=file_pattern,
                r1_suffix=r1_suffix,
                r2_suffix=r2_suffix,
                paired=paired
            )
            
            if files:
                samples[sample_id] = files
            else:
                logger.warning(f"No sequence files found for sample {sample_id}")
    
    # Log the results
    logger.info(f"Found sequence files for {len(samples)} out of {len(metadata)} samples")
    
    return samples

def interactive_file_pattern_collection(paired: bool = False) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Interactively collect file pattern information from the user.
    
    Args:
        paired: Whether to collect paired-end patterns
        
    Returns:
        Tuple of (file_pattern, r1_suffix, r2_suffix)
    """
    print("\nHow would you like to identify sequence files for each sample?")
    print("1. Using a pattern with {sample} placeholder")
    print("2. Using standard suffixes")
    
    choice = input("Enter choice (1 or 2): ").strip()
    
    if choice == "1":
        print("\nEnter a pattern using {sample} as a placeholder for the sample ID.")
        print("Examples:")
        print("  {sample}.fastq.gz")
        print("  {sample}_S*_R*.fastq.gz")
        
        file_pattern = input("Pattern: ").strip()
        return file_pattern, None, None
    
    elif choice == "2":
        if paired:
            print("\nEnter suffixes for paired-end files.")
            print("Examples:")
            print("  _R1.fastq.gz and _R2.fastq.gz")
            print("  _1.fq.gz and _2.fq.gz")
            
            r1_suffix = input("R1 suffix: ").strip()
            r2_suffix = input("R2 suffix: ").strip()
            
            return None, r1_suffix, r2_suffix
        else:
            print("\nEnter suffix for single-end files.")
            print("Examples:")
            print("  .fastq.gz")
            print("  _trimmed.fq.gz")
            
            r1_suffix = input("File suffix: ").strip()
            
            return None, r1_suffix, None
    
    else:
        print("Invalid choice. Using auto-detection.")
        return None, None, None

def main():
    parser = argparse.ArgumentParser(description="Metadata-driven sample collection for HUMAnN3 workflows")
    
    # Required arguments
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV file")
    parser.add_argument("--seq-dir", required=True, help="Directory containing sequence files")
    
    # Optional arguments
    parser.add_argument("--output", default="sample_files.txt", help="Output file to save sample file paths")
    parser.add_argument("--sample-col", help="Column name for sample identifiers")
    parser.add_argument("--r1-col", help="Column name for R1 file paths")
    parser.add_argument("--r2-col", help="Column name for R2 file paths")
    parser.add_argument("--file-pattern", help="Pattern for finding files (e.g., '{sample}_S*_R*.fastq.gz')")
    parser.add_argument("--r1-suffix", help="Suffix for R1 files (e.g., '_R1.fastq.gz')")
    parser.add_argument("--r2-suffix", help="Suffix for R2 files (e.g., '_R2.fastq.gz')")
    parser.add_argument("--paired", action="store_true", help="Whether samples are paired-end")
    parser.add_argument("--interactive", action="store_true", help="Interactively collect file pattern information")
    
    args = parser.parse_args()
    
    # Set up file logging
    file_handler = logging.FileHandler("sample_collection.log")
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(file_handler)
    
    # Interactively collect file pattern information if needed
    if args.interactive and not (args.r1_col or args.file_pattern or args.r1_suffix):
        file_pattern, r1_suffix, r2_suffix = interactive_file_pattern_collection(args.paired)
        
        if file_pattern:
            args.file_pattern = file_pattern
        if r1_suffix:
            args.r1_suffix = r1_suffix
        if r2_suffix:
            args.r2_suffix = r2_suffix
    
    # Collect samples from metadata
    samples = collect_samples_from_metadata(
        metadata_file=args.metadata,
        seq_dir=args.seq_dir,
        sample_col=args.sample_col,
        r1_col=args.r1_col,
        r2_col=args.r2_col,
        file_pattern=args.file_pattern,
        r1_suffix=args.r1_suffix,
        r2_suffix=args.r2_suffix,
        paired=args.paired
    )
    
    if not samples:
        logger.error("No sample files found. Exiting.")
        sys.exit(1)
    
    # Write results to output file
    with open(args.output, 'w') as f:
        for sample_id, files in samples.items():
            files_str = " ".join(files)
            f.write(f"{sample_id}\t{files_str}\n")
    
    logger.info(f"Sample file mapping saved to {args.output}")
    
    # Also print a summary
    print(f"\nSummary of collected samples:")
    print(f"Found sequence files for {len(samples)} samples")
    
    # Print first few samples as examples
    sample_preview = list(samples.items())[:5]
    print("\nExamples:")
    for sample_id, files in sample_preview:
        print(f"  Sample {sample_id}: {len(files)} file(s)")
        for file in files[:2]:  # Show max 2 files per sample
            print(f"    - {os.path.basename(file)}")
        if len(files) > 2:
            print(f"    - ... ({len(files)-2} more)")
    
    # Print command for using these files with the HUMAnN3 workflow
    print("\nTo use these files with the HUMAnN3 workflow:")
    file_list = []
    for files in samples.values():
        file_list.extend(files)
    
    print(f"\nhumann3-tools --run-preprocessing --input-fastq {' '.join(file_list[:3])} ... (total {len(file_list)} files)")
    print("\nOr to use the sample file with a future enhanced workflow:")
    print(f"humann3-tools --run-preprocessing --samples-file {args.output} ...\n")

if __name__ == "__main__":
    main()
