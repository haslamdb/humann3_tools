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

# Set up logging
logger = logging.getLogger('humann3_tools')
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
logger.addHandler(handler)

# Package version
__version__ = "0.1.0"

def setup_subparsers(parser):
    """
    Set up the subparsers for each command module.
    
    Args:
        parser: The main argument parser
        
    Returns:
        The parser with subparsers added
    """
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Define commands and their modules
    commands = {
        'kneaddata': 'KneadData preprocessing',
        'humann3': 'HUMAnN3 functional profiling',
        'join': 'Join and normalize HUMAnN3 output',
        'stats': 'Statistical testing',
        'diff': 'Differential abundance analysis',
        'viz': 'Visualization of results'
    }
    
    # Add subparsers for each command
    for command, description in commands.items():
        module_name = f'{command}_cli'
        try:
            # Import the module to get its parser
            module = importlib.import_module(f'humann3_tools.cli.{module_name}')
            subparser = subparsers.add_parser(command, help=description)
            
            # Get the module's argument parser
            if hasattr(module, 'parse_args'):
                # Copy arguments from the module's parser
                module_parser = argparse.ArgumentParser()
                module.parse_args([], parser=module_parser)
                for action in module_parser._actions:
                    if action.dest != 'help':  # Skip help action
                        subparser._add_action(action)
        except ImportError:
            # Just add a basic subparser if module can't be imported
            subparser = subparsers.add_parser(command, help=description)
    
    return parser

def main():
    """
    Main entry point for HUMAnN3 Tools CLI.
    
    Parses the command line arguments and dispatches to the appropriate module.
    """
    # Handle warnings
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    
    # Create the main parser
    parser = argparse.ArgumentParser(
        description="HUMAnN3 Tools - A comprehensive toolkit for metagenomic analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    # Add version argument
    parser.add_argument('--version', action='version', version=f'HUMAnN3 Tools v{__version__}')
    
    # Add subparsers for each command
    parser = setup_subparsers(parser)
    
    # Parse just enough to get the command
    args, remaining = parser.parse_known_args()
    
    # If no command was provided, show help
    if not args.command:
        parser.print_help()
        return 0
    
    # Import the module for the specified command
    module_name = f"humann3_tools.cli.{args.command}_cli"
    try:
        module = importlib.import_module(module_name)
    except ImportError as e:
        # If we're in development mode, try to import from local directory
        try:
            sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
            module = importlib.import_module(f"humann3_tools.cli.{args.command}_cli")
        except ImportError:
            logger.error(f"Module {module_name} not found: {e}")
            return 1
    
    # If --help was provided as the first argument, show command help
    if args.help or '--help' in remaining:
        if hasattr(module, 'parse_args'):
            module.parse_args(['--help'])
        else:
            parser.parse_args([args.command, '--help'])
        return 0
    
    # Run the module's main function
    if hasattr(module, 'main'):
        # If we have remaining arguments, parse them with the module's parser
        if args.command and remaining:
            module_parser = argparse.ArgumentParser()
            if hasattr(module, 'parse_args'):
                module_parser = module.parse_args([], parser=module_parser)
                module_args = module_parser.parse_args(remaining)
                return module.main(module_args)
            else:
                # Just pass the remaining arguments as they are
                return module.main()
        else:
            # No additional arguments provided
            return module.main()
    else:
        logger.error(f"Module {module_name} does not have a main function")
        return 1

if __name__ == "__main__":
    sys.exit(main())
