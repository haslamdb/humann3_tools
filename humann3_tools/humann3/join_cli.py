# humann3_tools/humann3/join_cli.py

import os
import argparse
import logging
import sys
import subprocess
import glob

from humann3_tools.logger import setup_logger, log_print
from humann3_tools.utils.cmd_utils import run_cmd

def parse_args():
    """Parse command line arguments for the join/normalize CLI."""
    parser = argparse.ArgumentParser(
        description="Join, normalize, and unstratify HUMAnN3 output files"
    )
    
    # Required arguments
    parser.add_argument("--input-dir", required=True, 
                      help="Directory containing HUMAnN3 output files")
    
    # File type options
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pathabundance", action="store_true", 
                     help="Process pathabundance files")
    group.add_argument("--pathcoverage", action="store_true", 
                     help="Process pathcoverage files")
    group.add_argument("--genefamilies", action="store_true", 
                     help="Process genefamilies files")
    
    # Output options
    parser.add_argument("--output-dir", default="./joined_output",
                      help="Directory for output files")
    parser.add_argument("--output-basename", 
                      help="Base filename for output (default derived from file type)")
    parser.add_argument("--units", default="cpm", choices=["cpm", "relab"],
                      help="Units for normalization (default: cpm)")
    
    # Additional options
    parser.add_argument("--update-snames", action="store_true",
                      help="Update sample names during normalization")
    parser.add_argument("--file-pattern", 
                      help="Glob pattern for input files (default based on file type)")
    parser.add_argument("--log-file", default=None,
                      help="Path to log file")
    parser.add_argument("--log-level", default="INFO", 
                      choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                      help="Logging level (default: INFO)")
    
    return parser.parse_args()

def main():
    """Main function to join, normalize, and unstratify HUMAnN3 output files."""
    args = parse_args()
    
    # Setup logging
    logger = setup_logger(log_file=args.log_file, 
                         log_level=getattr(logging, args.log_level.upper()))
    log_print("Starting join/normalize/unstratify process", level="info")
    
    # Check if input directory exists
    if not os.path.isdir(args.input_dir):
        log_print(f"ERROR: Input directory not found: {args.input_dir}", level="error")
        sys.exit(1)
    
    # Determine file type and pattern
    if args.pathabundance:
        file_type = "pathabundance"
        default_pattern = "*pathabundance.tsv"
    elif args.pathcoverage:
        file_type = "pathcoverage"
        default_pattern = "*pathcoverage.tsv"
    elif args.genefamilies:
        file_type = "genefamilies"
        default_pattern = "*genefamilies.tsv"
    
    # Use provided pattern or default
    pattern = args.file_pattern if args.file_pattern else default_pattern
    
    # Find input files
    search_pattern = os.path.join(args.input_dir, pattern)
    input_files = glob.glob(search_pattern)
    
    if not input_files:
        log_print(f"ERROR: No files found matching pattern: {search_pattern}", level="error")
        sys.exit(1)
    
    log_print(f"Found {len(input_files)} {file_type} files", level="info")
    
    # Create output directories
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Set output basename if not provided
    if not args.output_basename:
        args.output_basename = f"{file_type}_{args.units}"
    
    # Create normalized files directory
    norm_dir = os.path.join(args.output_dir, "normalized")
    os.makedirs(norm_dir, exist_ok=True)
    
    # 1. Normalize each file
    log_print(f"Normalizing files to {args.units} units...", level="info")
    normalized_files = []
    
    for input_file in input_files:
        basename = os.path.basename(input_file)
        sample_name = basename.replace(f"_{file_type}.tsv", "").replace(f".{file_type}.tsv", "")
        output_norm = os.path.join(norm_dir, f"{sample_name}_{file_type}_{args.units}.tsv")
        
        cmd = [
            "humann_renorm_table",
            "--input", input_file,
            "--output", output_norm,
            "--units", args.units
        ]
        
        if args.update_snames:
            cmd.append("--update-snames")
        
        success = run_cmd(cmd, exit_on_error=False)
        if success:
            normalized_files.append(output_norm)
            log_print(f"Normalized: {os.path.basename(input_file)} → {os.path.basename(output_norm)}", level="debug")
        else:
            log_print(f"Failed to normalize: {os.path.basename(input_file)}", level="warning")
    
    if not normalized_files:
        log_print("ERROR: No files were successfully normalized", level="error")
        sys.exit(1)
    
    log_print(f"Successfully normalized {len(normalized_files)} files", level="info")
    
    # 2. Join normalized files
    log_print("Joining normalized files...", level="info")
    joined_output = os.path.join(args.output_dir, f"{args.output_basename}.tsv")
    
    join_cmd = [
        "humann_join_tables",
        "-i", norm_dir,
        "-o", joined_output
    ]
    
    success = run_cmd(join_cmd, exit_on_error=False)
    if not success:
        log_print("ERROR: Failed to join normalized files", level="error")
        sys.exit(1)
    
    log_print(f"Successfully joined files: {joined_output}", level="info")
    
    # 3. Split stratified table
    log_print("Splitting stratified table...", level="info")
    
    split_cmd = [
        "humann_split_stratified_table",
        "-i", joined_output,
        "-o", args.output_dir
    ]
    
    success = run_cmd(split_cmd, exit_on_error=False)
    if not success:
        log_print("ERROR: Failed to split stratified table", level="error")
        sys.exit(1)
    
    # Find unstratified file
    unstrat_pattern = os.path.join(args.output_dir, f"*{args.output_basename}*_unstratified.tsv")
    unstrat_files = glob.glob(unstrat_pattern)
    
    if not unstrat_files:
        log_print("WARNING: Could not find unstratified output file", level="warning")
    else:
        unstrat_file = unstrat_files[0]
        log_print(f"Successfully created unstratified file: {os.path.basename(unstrat_file)}", level="info")
        
        # Rename the unstratified file to a standard name
        final_name = os.path.join(args.output_dir, f"{args.output_basename}_unstratified.tsv")
        try:
            if os.path.abspath(unstrat_file) != os.path.abspath(final_name):
                os.rename(unstrat_file, final_name)
                log_print(f"Renamed to: {os.path.basename(final_name)}", level="info")
        except Exception as e:
            log_print(f"WARNING: Could not rename file: {str(e)}", level="warning")
    
    log_print("Join/normalize/unstratify process complete", level="info")
    return 0

if __name__ == "__main__":
    sys.exit(main())
