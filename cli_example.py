#!/usr/bin/env python3
"""
Simple CLI example
"""

import argparse

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="Simple CLI example")
    parser.add_argument("--input", help="Input file or directory")
    parser.add_argument("--output", help="Output file or directory")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads")
    
    # Add subcommands
    subparsers = parser.add_subparsers(dest="command", help="Command to run")
    
    # Add humann3 subcommand
    humann3_parser = subparsers.add_parser("humann3", help="Run HUMAnN3")
    humann3_parser.add_argument("--input-dir", help="Input directory")
    humann3_parser.add_argument("--output-dir", help="Output directory")
    
    # Add join subcommand
    join_parser = subparsers.add_parser("join", help="Join files")
    join_parser.add_argument("--input-dir", help="Input directory")
    join_parser.add_argument("--output-dir", help="Output directory")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Print arguments
    print("Arguments:")
    for arg in vars(args):
        print(f"  {arg}: {getattr(args, arg)}")
    
    # Run command
    if args.command == "humann3":
        print("Running HUMAnN3...")
    elif args.command == "join":
        print("Joining files...")
    else:
        parser.print_help()
    
    return 0

if __name__ == "__main__":
    main()