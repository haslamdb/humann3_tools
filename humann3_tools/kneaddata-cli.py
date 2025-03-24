#!/usr/bin/env python3
"""
HUMAnN3 Tools KneadData CLI

This standalone script handles the trimming and filtering of reads using KneadData.
It can be used independently from the main humann3_tools pipeline.

Examples:
  # Basic usage with paired-end reads:
  humann3-kneaddata --input-fastq sample1_R1.fastq.gz sample1_R2.fastq.gz --paired --threads 8 --output-dir ./kneaddata_output
  
  # Running with multiple samples:
  humann3-kneaddata --input-fastq sample1_R1.fastq.gz sample1_R2.fastq.gz sample2_R1.fastq.gz sample2_R2.fastq.gz --paired --threads 8
  
  # Specifying database path:
  humann3-kneaddata --input-fastq sample1.fastq.gz --reference-dbs /path/to/kneaddata_db --output-dir ./kneaddata_output

Dependencies:
  - KneadData (v0.7.0+)
  - Python 3.6+
"""

import os
import sys
import time
import argparse
import logging
from typing import List, Dict, Optional

from humann3_tools.logger import setup_logger, log_print
from humann3_tools.preprocessing.kneaddata import (
    check_kneaddata_installation,
    run_kneaddata,
    run_kneaddata_parallel
)
from humann3_tools.utils.metadata_utils import (
    collect_samples_from_metadata,
    find_sample_files
)

def parse_args():
    """Parse command line arguments for the KneadData CLI."""
    parser = argparse.ArgumentParser(
        description="Run KneadData for quality control and host depletion"
    )
    
    # Required arguments
    parser.add_argument("--input-fastq", nargs="+", required=True,
                      help="Input FASTQ file(s) for preprocessing")
    
    # Output options
    parser.add_argument("--output-dir", default="./kneaddata_output",
                      help="Directory for output files")
    
    # KneadData options
    parser.add_argument("--reference-dbs", nargs="+", required=True,
                      help="Path(s) to KneadData reference database(s)")
    parser.add_argument("--paired", action="store_true",
                      help="Input files are paired-end reads")
    parser.add_argument("--decontaminate-pairs", default="strict", 
                      choices=["strict", "lenient", "unpaired"],
                      help="Method for decontaminating paired-end reads (default: strict)")
    
    # Performance options
    parser.add_argument("--threads", type=int, default=1,
                      help="Number of threads to use")
    parser.add_argument("--use-parallel", action="store_true",
                      help="Use parallel processing (process multiple samples simultaneously)")
    parser.add_argument("--max-parallel", type=int, default=None,
                      help="Maximum number of samples to process in parallel")
    
    # Metadata-driven options
    parser.add_argument("--use-metadata", action="store_true",
                      help="Use metadata file to locate sequence files")
    parser.add_argument("--metadata-file",
                      help="Path to metadata CSV file (when using --use-metadata)")
    parser.add_argument("--seq-dir",
                      help="Directory containing sequence files (when using --use-metadata)")
    parser.add_argument("--sample-col",
                      help="Column name for sample IDs in metadata (when using --use-metadata)")
    parser.add_argument("--r1-suffix",
                      help="Suffix for R1 sequence file paths (when using --use-metadata with paired data)")
    parser.add_argument("--r2-suffix",
                      help="Suffix for R2 sequence file paths (when using --use-metadata with paired data)")
    parser.add_argument("--file-pattern",
                      help="File pattern to match for samples (can use {sample})")
    
    # Logging options
    parser.add_argument("--log-file",
                      help="Path to log file")
    parser.add_argument("--log-level", default="INFO",
                      choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                      help="Logging level")
    
    # Additional KneadData options
    parser.add_argument("--additional-options", nargs="+",
                      help="Additional options to pass to KneadData (format: key=value)")
    
    return parser.parse_args()

def main():
    """Main function to run KneadData."""
    # Parse arguments
    args = parse_args()
    
    # Setup logging
    log_level = getattr(logging, args.log_level.upper())
    logger = setup_logger(log_file=args.log_file, log_level=log_level)
    
    start_time = time.time()
    log_print("Starting HUMAnN3 Tools KneadData Pipeline", level="info")
    
    # Check KneadData installation
    kneaddata_ok, kneaddata_version = check_kneaddata_installation()
    if not kneaddata_ok:
        log_print(f"ERROR: KneadData not properly installed: {kneaddata_version}", level="error")
        sys.exit(1)
    log_print(f"Using KneadData version: {kneaddata_version}", level="info")
    
    # Process input files
    input_files = args.input_fastq
    
    # If using metadata-driven approach
    if args.use_metadata:
        if not args.metadata_file or not args.seq_dir:
            log_print("ERROR: --metadata-file and --seq-dir are required when using --use-metadata", level="error")
            sys.exit(1)
            
        samples_dict = collect_samples_from_metadata(
            metadata_file=args.metadata_file,
            seq_dir=args.seq_dir,
            sample_col=args.sample_col,
            r1_suffix=args.r1_suffix,
            r2_suffix=args.r2_suffix,
            file_pattern=args.file_pattern,
            paired=args.paired
        )
        
        # Reset input_files
        input_files = []
        
        for sample_id, files in samples_dict.items():
            input_files.extend(files)
            log_print(f"Added files for sample {sample_id}: {[os.path.basename(f) for f in files]}", level="info")
        
        log_print(f"Collected {len(input_files)} sequence files from metadata", level="info")
        
        if not input_files:
            log_print("ERROR: No input files collected from metadata. Cannot proceed.", level="error")
            sys.exit(1)
    
    # Prepare KneadData options
    kneaddata_options = {}
    if args.paired:
        kneaddata_options["decontaminate-pairs"] = args.decontaminate_pairs
    
    # Add any additional options
    if args.additional_options:
        for option in args.additional_options:
            if '=' in option:
                key, value = option.split('=', 1)
                kneaddata_options[key] = value
            else:
                kneaddata_options[option] = True
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Run KneadData
    if args.use_parallel:
        # Parallel processing
        log_print("Running KneadData in parallel mode", level="info")
        results = run_kneaddata_parallel(
            input_files=input_files,
            output_dir=args.output_dir,
            threads=args.threads,
            max_parallel=args.max_parallel,
            reference_dbs=args.reference_dbs,
            paired=args.paired,
            additional_options=kneaddata_options,
            logger=logger
        )
    else:
        # Standard processing
        log_print("Running KneadData in standard mode", level="info")
        results = run_kneaddata(
            input_files=input_files,
            output_dir=args.output_dir,
            threads=args.threads,
            reference_dbs=args.reference_dbs, 
            paired=args.paired,
            additional_options=kneaddata_options,
            logger=logger
        )
    
    # Generate a summary of results
    log_print("\nKneadData Pipeline Summary:", level="info")
    
    if results:
        sample_count = len(results)
        output_count = sum(len(files) for files in results.values())
        log_print(f"Successfully processed {sample_count} samples", level="info")
        log_print(f"Generated {output_count} output files", level="info")
        
        # Print output locations
        log_print("\nOutput Directory:", level="info")
        log_print(f"  {os.path.abspath(args.output_dir)}", level="info")
        
        # Print next steps
        log_print("\nNext Steps:", level="info")
        log_print("  To run HUMAnN3 on these files, use:", level="info")
        log_print("  humann3-run --input-files [kneaddata_output_files] --output-dir [humann3_output_dir]", level="info")
    else:
        log_print("No samples were successfully processed", level="warning")
    
    # Print elapsed time
    elapsed_time = time.time() - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    hours, minutes = divmod(minutes, 60)
    
    log_print(f"\nTotal processing time: {int(hours)}h {int(minutes)}m {int(seconds)}s", level="info")
    
    return 0 if results else 1

if __name__ == "__main__":
    sys.exit(main())
