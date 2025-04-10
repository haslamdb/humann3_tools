#!/usr/bin/env python3
"""
HUMAnN3 Tools - Simplified HUMAnN3 CLI Module
"""

import argparse
import sys
import os

def parse_args(args=None, parent_parser=None):
    """
    Parse command line arguments for the HUMAnN3 module.
    
    Args:
        args: List of arguments to parse (default: sys.argv[1:])
        parent_parser: Parent parser to add arguments to (default: create new parser)
        
    Returns:
        Parsed arguments if args is provided, otherwise the configured parser
    """
    # Use existing parser or create new one
    parser = parent_parser or argparse.ArgumentParser(
        description="Run HUMAnN3 on preprocessed sequence files"
    )
    
    # Add arguments
    parser.add_argument("--input-dir", help="Directory containing KneadData output files")
    parser.add_argument("--output-dir", default="./humann3_output", help="Output directory")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads")
    parser.add_argument("--paired", action="store_true", help="Input files are paired-end reads")
    
    # Return the parser or parsed args
    if args is not None:
        return parser.parse_args(args)
    return parser

def main(args=None):
    """
    Main function to run HUMAnN3 processing.
    
    Args:
        args: Command line arguments (optional)
    """
    # Parse arguments
    if isinstance(args, list):
        args = parse_args(args)
    else:
        args = parse_args()
    
    # Print the arguments
    print("Running HUMAnN3 processing...")
    for arg, value in vars(args).items():
        print(f"  {arg}: {value}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())