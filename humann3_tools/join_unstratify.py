# humann3_tools/humann3_tools/join_unstratify.py
import os
import sys
import argparse
import logging
import time

from humann3_tools.logger import setup_logger, log_print
from humann3_tools.utils.sample_utils import validate_sample_key, check_input_files_exist
from humann3_tools.humann3.pathway_processing import process_pathway_abundance
from humann3_tools.humann3.gene_processing import process_gene_families

def parse_args():
    """Parse command line arguments for join_unstratify_humann_output."""
    parser = argparse.ArgumentParser(
        description="Join and unstratify HUMAnN3 gene family and pathway abundance files without running downstream analysis"
    )
    
    # Required arguments
    parser.add_argument("--sample-key", required=True, help="CSV file with columns for sample names and metadata")
    parser.add_argument("--pathway-dir", required=True, help="Directory containing raw pathway abundance files")
    parser.add_argument("--gene-dir", required=True, help="Directory containing raw gene family files")
    
    # Output options
    parser.add_argument("--output-dir", default="./HUMAnN3_Processed", help="Directory for output files")
    parser.add_argument("--output-prefix", default="ProcessedFiles", help="Prefix for output files")
    
    # Processing options
    parser.add_argument("--skip-pathway", action="store_true", help="Skip pathway processing")
    parser.add_argument("--skip-gene", action="store_true", help="Skip gene family processing")
    parser.add_argument("--units", default="cpm", choices=["cpm", "relab"], help="Units for normalization (default: cpm)")
    
    # Additional options
    parser.add_argument("--no-interactive", action="store_true", help="Non-interactive mode for sample key column selection")
    parser.add_argument("--log-file", default=None, help="Path to log file")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], 
                      help="Logging level (default: INFO)")
    
    return parser.parse_args()

def join_unstratify_humann_output():
    """
    Main function to join and unstratify HUMAnN3 gene family and pathway files.
    This function only processes the HUMAnN3 output files without running 
    preprocessing or downstream analysis.
    """
    args = parse_args()
    
    # Setup logging
    logger = setup_logger(log_file=args.log_file, log_level=getattr(logging, args.log_level))
    log_print("Starting HUMAnN3 Join and Unstratify Pipeline", level="info")
    start_time = time.time()
    
    # Process sample metadata
    samples, selected_columns = validate_sample_key(args.sample_key, no_interactive=args.no_interactive)
    
    # Check input files
    valid_path_samples, valid_gene_samples = check_input_files_exist(samples, args.pathway_dir, args.gene_dir)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Process pathway files
    pathway_unstrat_file = None
    if not args.skip_pathway and valid_path_samples:
        log_print("Processing pathway abundance files...", level="info")
        pathway_unstrat_file = process_pathway_abundance(
            valid_path_samples,
            args.pathway_dir,
            args.output_dir,
            args.output_prefix,
            selected_columns=selected_columns,
            units=args.units
        )
        
        if pathway_unstrat_file:
            log_print(f"Pathway processing completed. Unstratified file: {pathway_unstrat_file}", level="info")
        else:
            log_print("Pathway processing failed to produce an unstratified file", level="warning")
    else:
        if args.skip_pathway:
            log_print("Skipping pathway processing as requested", level="info")
        else:
            log_print("No valid pathway files found; skipping pathway processing", level="warning")
    
    # Process gene family files
    gene_unstrat_file = None
    if not args.skip_gene and valid_gene_samples:
        log_print("Processing gene family files...", level="info")
        gene_unstrat_file = process_gene_families(
            valid_gene_samples,
            args.gene_dir,
            args.output_dir,
            args.output_prefix,
            selected_columns=selected_columns,
            units=args.units
        )
        
        if gene_unstrat_file:
            log_print(f"Gene family processing completed. Unstratified file: {gene_unstrat_file}", level="info")
        else:
            log_print("Gene family processing failed to produce an unstratified file", level="warning")
    else:
        if args.skip_gene:
            log_print("Skipping gene family processing as requested", level="info")
        else:
            log_print("No valid gene family files found; skipping gene family processing", level="warning")
    
    # Report summary
    elapsed = time.time() - start_time
    mm, ss = divmod(elapsed, 60)
    
    log_print(f"Pipeline completed in {int(mm)}m {int(ss)}s", level="info")
    log_print("Summary:", level="info")
    
    if pathway_unstrat_file:
        log_print(f"  - Pathway unstratified file: {pathway_unstrat_file}", level="info")
    if gene_unstrat_file:
        log_print(f"  - Gene family unstratified file: {gene_unstrat_file}", level="info")
    
    if not pathway_unstrat_file and not gene_unstrat_file:
        log_print("  - No output files were generated", level="warning")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(join_unstratify_humann_output())
