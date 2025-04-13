#!/usr/bin/env python3
"""
HUMAnN3 Tools - Main CLI Interface

A comprehensive toolkit for metagenomic analysis with HUMAnN3, from quality control to visualization.

Available Commands:
  kneaddata   - Quality control and host depletion using KneadData
  humann3     - Run HUMAnN3 on preprocessed sequence files
  join        - Join, normalize, and unstratify HUMAnN3 output files
  stats       - Run statistical tests on HUMAnN3 output data
  diff        - Run differential abundance analysis
  viz         - Create visualizations from HUMAnN3 output data

Example Usage:
  humann3-tools kneaddata --input-files sample_R1.fastq.gz sample_R2.fastq.gz --paired --reference-dbs human
  humann3-tools humann3 --input-dir kneaddata_output --output-dir humann3_output
  humann3-tools join --input-dir PathwayAbundance --pathabundance --output-dir joined_output
  humann3-tools stats --abundance-file joined_output/pathway_abundance_cpm_unstratified.tsv --metadata-file metadata.csv
  humann3-tools diff --abundance-file joined_output/pathway_abundance_cpm_unstratified.tsv --metadata-file metadata.csv
  humann3-tools viz --abundance-file joined_output/pathway_abundance_cpm_unstratified.tsv --metadata-file metadata.csv

For more information on any command, use:
  humann3-tools [command] --help
"""

import os
import sys
import argparse
import importlib
import logging
import warnings
import pkg_resources

# Import CLI modules
from src.humann3_tools.cli import kneaddata_cli
from src.humann3_tools.cli import humann3_cli
from src.humann3_tools.cli import join_cli
from src.humann3_tools.cli import stats_cli
from src.humann3_tools.cli import diff_cli
from src.humann3_tools.cli import viz_cli

# Set up logging
logger = logging.getLogger('humann3_tools')
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
logger.addHandler(handler)

# Package version - avoid pkg_resources for compatibility
#__version__ = "0.2.0"  # Hardcode version to avoid dependency issues

def setup_subparsers(parser):
    """
    Set up the subparsers for each command module.
    
    Args:
        parser: The main argument parser
        
    Returns:
        The parser with subparsers added
    """
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Add basic subparsers directly
    
    # Humann3 subparser
    humann3_parser = subparsers.add_parser('humann3', help='HUMAnN3 functional profiling')
    humann3_parser.add_argument('--input-dir', help="Directory containing KneadData output files")
    humann3_parser.add_argument('--output-dir', help="Output directory")
    humann3_parser.add_argument('--threads', type=int, default=1, help="Number of threads")
    humann3_parser.add_argument('--input-files', nargs="+", help="Input FASTQ files")
    humann3_parser.add_argument('--paired', action='store_true', help="Input is paired-end reads")
    
    # Join subparser
    join_parser = subparsers.add_parser('join', help='Join and normalize HUMAnN3 output')
    join_parser.add_argument('--input-dir', required=True, help="Directory containing HUMAnN3 output files")
    join_parser.add_argument('--output-dir', default="./joined_output", help="Output directory")
    join_parser.add_argument('--pathabundance', action='store_true', help="Process pathway abundance files")
    join_parser.add_argument('--pathcoverage', action='store_true', help="Process pathway coverage files")
    join_parser.add_argument('--genefamilies', action='store_true', help="Process gene family files")
    
    # Stats subparser
    stats_parser = subparsers.add_parser('stats', help='Statistical testing')
    stats_parser.add_argument('--abundance-file', help="Abundance file (from join step)")
    stats_parser.add_argument('--metadata-file', help="Metadata file (CSV)")
    
    # Diff subparser
    diff_parser = subparsers.add_parser('diff', help='Differential abundance analysis')
    diff_parser.add_argument('--abundance-file', help="Abundance file (from join step)")
    diff_parser.add_argument('--metadata-file', help="Metadata file (CSV)")
    
    # Viz subparser
    viz_parser = subparsers.add_parser('viz', help='Visualization of results')
    viz_parser.add_argument('--abundance-file', help="Abundance file (from join step)")
    viz_parser.add_argument('--metadata-file', help="Metadata file (CSV)")
    
    # Kneaddata subparser
    kneaddata_parser = subparsers.add_parser('kneaddata', help='KneadData preprocessing')
    kneaddata_parser.add_argument('--input-files', nargs="+", help="Input FASTQ files")
    kneaddata_parser.add_argument('--output-dir', help="Output directory")
    kneaddata_parser.add_argument('--paired', action='store_true', help="Input is paired-end reads")
    
    return parser

def main():
    """
    Main entry point for HUMAnN3 Tools CLI.
    
    Parses the command line arguments and dispatches to the appropriate module.
    """
    # Handle warnings
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    
    # First just parse the command name
    if len(sys.argv) < 2:
        # No command provided, show help
        parser = argparse.ArgumentParser(
            description="HUMAnN3 Tools - A comprehensive toolkit for metagenomic analysis",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=__doc__
        )
        parser.add_argument('--version', action='version', version=f'HUMAnN3 Tools v{__version__}')
        parser.print_help()
        return 0
    
    command = sys.argv[1]
    
    # Remove the command from argv for further processing
    command_args = sys.argv[2:]
    
    # Execute appropriate command
    if command == 'humann3':
        logger.info(f"Running HUMAnN3...")
        return humann3_cli.main()
        
    elif command == 'join':
        logger.info(f"Running Join...")
        return join_cli.main()
        
    elif command == 'kneaddata':
        logger.info("Running KneadData...")
        return kneaddata_cli.main()
        
    elif command == 'stats':
        logger.info("Running Statistics...")
        return stats_cli.main()
        
    elif command == 'diff':
        logger.info("Running Differential Analysis...")
        return diff_cli.main()
        
    elif command == 'viz':
        logger.info("Running Visualization...")
        return viz_cli.main()
    
    elif command == '--help' or command == '-h':
        # Show help
        parser = argparse.ArgumentParser(
            description="HUMAnN3 Tools - A comprehensive toolkit for metagenomic analysis",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=__doc__
        )
        parser.add_argument('--version', action='version', version=f'HUMAnN3 Tools v{__version__}')
        parser.print_help()
        return 0
        
    else:
        logger.error(f"Unknown command: {command}")
        # Show available commands
        print("\nAvailable commands:")
        print("  humann3    - Run HUMAnN3 on preprocessed sequence files")
        print("  join       - Join and normalize HUMAnN3 output files")
        print("  kneaddata  - Quality control and host depletion using KneadData")
        print("  stats      - Run statistical tests on HUMAnN3 output data")
        print("  diff       - Run differential abundance analysis")
        print("  viz        - Create visualizations from HUMAnN3 output data")
        print("\nFor more information on any command, use:")
        print("  humann3-tools [command] --help")
        return 1

if __name__ == "__main__":
    sys.exit(main())