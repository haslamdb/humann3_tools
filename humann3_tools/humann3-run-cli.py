#!/usr/bin/env python3
"""
HUMAnN3 Tools HUMAnN3 Runner CLI

This standalone script runs HUMAnN3 on processed sequence files, typically
coming from KneadData. It can be used independently from the main pipeline.

Examples:
  # Basic usage with single sequence file:
  humann3-run --input-files cleaned_sample.fastq --threads 8 --output-dir ./humann3_output
  
  # Running with multiple samples:
  humann3-run --input-files sample1.fastq sample2.fastq --threads 8 --output-dir ./humann3_output
  
  # Specifying database paths:
  humann3-run --input-files sample.fastq --nucleotide-db /path/to/chocophlan --protein-db /path/to/uniref

Dependencies:
  - HUMAnN3 (v3.0.0+)
  - Python 3.6+
"""

import os
import sys
import time
import argparse
import logging
from typing import List, Dict, Optional

from humann3_tools.logger import setup_logger, log_print
from humann3_tools.preprocessing.humann3_run import (
    check_humann3_installation,
    run_humann3,
    run_humann3_parallel
)

def parse_args():
    """Parse command line arguments for the HUMAnN3 CLI."""
    parser = argparse.ArgumentParser(
        description="Run HUMAnN3 on processed sequence files"
    )
    
    # Required arguments
    parser.add_argument("--input-files", nargs="+", required=True,
                      help="Input sequence file(s) (typically from KneadData)")
    
    # Output options
    parser.add_argument("--output-dir", default="./humann3_output",
                      help="Directory for output files")
    parser.add_argument("--pathabundance-dir",
                      help="Directory for pathway abundance files (default: {output-dir}/PathwayAbundance)")
    parser.add_argument("--pathcoverage-dir",
                      help="Directory for pathway coverage files (default: {output-dir}/PathwayCoverage)")
    parser.add_argument("--genefamilies-dir",
                      help="Directory for gene families files (default: {output-dir}/GeneFamilies)")  
    parser.add_argument("--metaphlan-dir",
                      help="Directory for metaphlan bugs list files (default: {output-dir}/MetaphlanFiles)")
    
    # HUMAnN3 database options
    parser.add_argument("--nucleotide-db",
                      help="Path to HUMAnN3 nucleotide database (ChocoPhlAn)")
    parser.add_argument("--protein-db",
                      help="Path to HUMAnN3 protein database (UniRef)")
    
    # Performance options
    parser.add_argument("--threads", type=int, default=1,
                      help="Number of threads to use")
    parser.add_argument("--use-parallel", action="store_true",
                      help="Use parallel processing (process multiple samples simultaneously)")
    parser.add_argument("--max-parallel", type=int, default=None,
                      help="Maximum number of samples to process in parallel")
    
    # Logging options
    parser.add_argument("--log-file",
                      help="Path to log file")
    parser.add_argument("--log-level", default="INFO",
                      choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                      help="Logging level")
    
    # Additional HUMAnN3 options
    parser.add_argument("--bypass-prescreen", action="store_true",
                      help="Bypass the MetaPhlAn taxonomic prescreen")
    parser.add_argument("--bypass-nucleotide-index", action="store_true",
                      help="Bypass the nucleotide index database")
    parser.add_argument("--bypass-translated-search", action="store_true",
                      help="Bypass the translated search")
    parser.add_argument("--additional-options", nargs="+",
                      help="Additional options to pass to HUMAnN3 (format: key=value)")
    
    return parser.parse_args()

def main():
    """Main function to run HUMAnN3."""
    # Parse arguments
    args = parse_args()
    
    # Setup logging
    log_level = getattr(logging, args.log_level.upper())
    logger = setup_logger(log_file=args.log_file, log_level=log_level)
    
    start_time = time.time()
    log_print("Starting HUMAnN3 Tools HUMAnN3 Runner", level="info")
    
    # Check HUMAnN3 installation
    humann3_ok, humann3_version = check_humann3_installation()
    if not humann3_ok:
        log_print(f"ERROR: HUMAnN3 not properly installed: {humann3_version}", level="error")
        sys.exit(1)
    log_print(f"Using HUMAnN3 version: {humann3_version}", level="info")
    
    # Process input files
    input_files = args.input_files
    if not input_files:
        log_print("ERROR: No input files specified", level="error")
        sys.exit(1)
    
    # Verify input files exist
    for input_file in input_files:
        if not os.path.exists(input_file):
            log_print(f"ERROR: Input file not found: {input_file}", level="error")
            sys.exit(1)
    
    log_print(f"Processing {len(input_files)} input files", level="info")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Setup output directories
    pathabundance_dir = args.pathabundance_dir or os.path.join(args.output_dir, "PathwayAbundance")
    pathcoverage_dir = args.pathcoverage_dir or os.path.join(args.output_dir, "PathwayCoverage")
    genefamilies_dir = args.genefamilies_dir or os.path.join(args.output_dir, "GeneFamilies")
    metaphlan_dir = args.metaphlan_dir or os.path.join(args.output_dir, "MetaphlanFiles")
    
    os.makedirs(pathabundance_dir, exist_ok=True)
    os.makedirs(pathcoverage_dir, exist_ok=True)
    os.makedirs(genefamilies_dir, exist_ok=True)
    os.makedirs(metaphlan_dir, exist_ok=True)
    
    # Prepare HUMAnN3 options
    humann3_options = {}
    
    if args.bypass_prescreen:
        humann3_options["bypass-prescreen"] = True
    if args.bypass_nucleotide_index:
        humann3_options["bypass-nucleotide-index"] = True
    if args.bypass_translated_search:
        humann3_options["bypass-translated-search"] = True
    
    # Add any additional options
    if args.additional_options:
        for option in args.additional_options:
            if '=' in option:
                key, value = option.split('=', 1)
                humann3_options[key] = value
            else:
                humann3_options[option] = True
    
    # Run HUMAnN3
    if args.use_parallel:
        # Parallel processing
        log_print("Running HUMAnN3 in parallel mode", level="info")
        results = run_humann3_parallel(
            input_files=input_files,
            output_dir=args.output_dir,
            threads=args.threads,
            max_parallel=args.max_parallel,
            nucleotide_db=args.nucleotide_db,
            protein_db=args.protein_db,
            additional_options=humann3_options,
            logger=logger
        )
    else:
        # Standard processing
        log_print("Running HUMAnN3 in standard mode", level="info")
        results = run_humann3(
            input_files=input_files,
            output_dir=args.output_dir,
            threads=args.threads,
            nucleotide_db=args.nucleotide_db,
            protein_db=args.protein_db,
            additional_options=humann3_options,
            logger=logger,
            pathabdirectory=pathabundance_dir,
            genedirectory=genefamilies_dir,
            pathcovdirectory=pathcoverage_dir,
            metadirectory=metaphlan_dir
        )
    
    # Generate a summary of results
    log_print("\nHUMAnN3 Runner Summary:", level="info")
    
    if results:
        sample_count = len(results)
        
        # Count output files by type
        path_count = sum(1 for outputs in results.values() if isinstance(outputs, dict) and outputs.get('pathabundance'))
        gene_count = sum(1 for outputs in results.values() if isinstance(outputs, dict) and outputs.get('genefamilies'))
        pathcov_count = sum(1 for outputs in results.values() if isinstance(outputs, dict) and outputs.get('pathcoverage'))
        metaphlan_count = sum(1 for outputs in results.values() if isinstance(outputs, dict) and outputs.get('metaphlan'))
        
        log_print(f"Successfully processed {sample_count} samples", level="info")
        log_print(f"Generated:", level="info")
        log_print(f"  - {path_count} pathway abundance files", level="info")
        log_print(f"  - {gene_count} gene families files", level="info")
        log_print(f"  - {pathcov_count} pathway coverage files", level="info")
        log_print(f"  - {metaphlan_count} MetaPhlAn bugs list files", level="info")
        
        # Print output locations
        log_print("\nOutput Directories:", level="info")
        log_print(f"  Pathway Abundance: {os.path.abspath(pathabundance_dir)}", level="info")
        log_print(f"  Gene Families: {os.path.abspath(genefamilies_dir)}", level="info")
        log_print(f"  Pathway Coverage: {os.path.abspath(pathcoverage_dir)}", level="info")
        log_print(f"  MetaPhlAn Bugs List: {os.path.abspath(metaphlan_dir)}", level="info")
        
        # Print next steps
        log_print("\nNext Steps:", level="info")
        log_print("  To join and process these files, use:", level="info")
        log_print("  humann3-join --input-dir PathwayAbundance --pathabundance --output-dir joined_output", level="info")
        log_print("  humann3-join --input-dir GeneFamilies --genefamilies --output-dir joined_output", level="info")
        log_print("  or use the combined command:", level="info")
        log_print("  join_unstratify_humann_output --pathway-dir PathwayAbundance --gene-dir GeneFamilies --output-dir joined_output", level="info")
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
