# humann3_tools/humann3_tools/cli.py
import os
import sys
import argparse
import time
import logging
import traceback

from humann3_tools.logger import setup_logger, log_print
from humann3_tools.utils.sample_utils import validate_sample_key, check_input_files_exist
from humann3_tools.humann3.pathway_processing import process_pathway_abundance
from humann3_tools.humann3.gene_processing import process_gene_families
from humann3_tools.analysis.metadata import read_and_process_metadata
from humann3_tools.analysis.visualizations import read_and_process_gene_families, read_and_process_pathways
from humann3_tools.analysis.statistical import run_statistical_tests
from humann3_tools.utils.file_utils import check_file_exists_with_logger
from humann3_tools.preprocessing.pipeline import (
    run_preprocessing_pipeline,
    run_preprocessing_pipeline_parallel,
)
from humann3_tools.preprocessing.kneaddata import check_kneaddata_installation
from humann3_tools.preprocessing.humann3_run import check_humann3_installation
from humann3_tools.utils.resource_utils import limit_memory_usage
from humann3_tools.utils.metadata_utils import collect_samples_from_metadata, find_sample_files


def main():
    """Main entry point for the humann3_tools CLI."""
    parser = argparse.ArgumentParser(description="HUMAnN3 Tools: Process and analyze HUMAnN3 output")

    # --- 1. Global Options ---
    global_group = parser.add_argument_group("Global Options")
    global_group.add_argument("--log-file", default=None, help="Path to combined log file")
    global_group.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default=INFO)",
    )
    global_group.add_argument(
        "--max-memory", type=int, default=None, help="Maximum memory usage in MB (default: unlimited)"
    )
    global_group.add_argument(
        "--list-files", action="store_true", help="Just list input files in --pathway-dir and --gene-dir, then exit"
    )
    global_group.add_argument(
        "--no-interactive", action="store_true", help="Non-interactive mode for sample key column selection"
    )

    # --- 2. Input/Output Options ---
    io_group = parser.add_argument_group("Input/Output Options")
    io_group.add_argument("--sample-key", required=True, help="CSV file with columns for sample names and metadata")
    io_group.add_argument(
        "--output-dir",
        default="./Humann3Output",
        help="Directory where HUMAnN3-processed files and downstream results will go",
    )
    io_group.add_argument(
        "--output-prefix",
        default="ProcessedFiles",
        help="Prefix for intermediate HUMAnN3 output directories/files",
    )

    # --- 3. Preprocessing (KneadData/HUMAnN3) Options ---
    preprocessing_group = parser.add_argument_group("Preprocessing Options")
    preprocessing_group.add_argument(
        "--run-preprocessing", action="store_true", help="Run preprocessing (KneadData and HUMAnN3) on raw sequence files"
    )
    preprocessing_group.add_argument("--input-fastq", nargs="+", help="Input FASTQ file(s) for preprocessing")
    preprocessing_group.add_argument("--paired", action="store_true", help="Input files are paired-end reads")
    preprocessing_group.add_argument(
        "--kneaddata-dbs", nargs="+", help="Path(s) to KneadData reference database(s). Can specify multiple databases."
    )
    preprocessing_group.add_argument("--humann3-nucleotide-db", help="Path to HUMAnN3 nucleotide database (ChocoPhlAn)")
    preprocessing_group.add_argument("--humann3-protein-db", help="Path to HUMAnN3 protein database (UniRef)")
    preprocessing_group.add_argument(
        "--threads", type=int, default=1, help="Number of threads to use for preprocessing"
    )
    preprocessing_group.add_argument(
        "--kneaddata-output-dir",
        help="Directory for KneadData output files (default: {output-dir}/PreprocessedData/kneaddata_output)",
    )
    preprocessing_group.add_argument(
        "--humann3-output-dir",
        help="Directory for HUMAnN3 output files (default: {output-dir}/PreprocessedData/humann3_output)",
    )

    # --- 3.1 Metadata driven workflow ---
    metadata_group = parser.add_argument_group("Metadata-driven workflow options")
    metadata_group.add_argument("--use-metadata", action="store_true", help="Read samples and file paths from metadata")
    metadata_group.add_argument(
        "--seq-dir", help="Directory containing sequence files (when using --use-metadata)"
    )
    metadata_group.add_argument(
        "--sample-col", help="Column name for sample IDs in metadata (when using --use-metadata)"
    )
    metadata_group.add_argument("--r1-col", help="Column name for R1 sequence file paths")
    metadata_group.add_argument("--r2-col", help="Column name for R2 sequence file paths")
    metadata_group.add_argument(
        "--file-pattern", help="File pattern to match for samples (can use {sample})"
    )
    metadata_group.add_argument("--r1-suffix", help="Suffix to append for R1 sequence file paths")
    metadata_group.add_argument("--r2-suffix", help="Suffix to append for R2 sequence file paths")
    metadata_group.add_argument(
        "--samples-file", help="Tab-delimited file with sample IDs and sequence file paths"
    )

    # --- 4. HUMAnN3 Processing Options ---
    humann3_group = parser.add_argument_group("HUMAnN3 Processing Options")
    humann3_group.add_argument(
        "--pathway-dir", required=True, help="Directory containing raw pathway abundance files"
    )
    humann3_group.add_argument("--gene-dir", required=True, help="Directory containing raw gene family files")
    humann3_group.add_argument("--skip-pathway", action="store_true", help="Skip HUMAnN3 pathway processing")
    humann3_group.add_argument("--skip-gene", action="store_true", help="Skip HUMAnN3 gene family processing")
    humann3_group.add_argument(
        "--annotations-dir", default="annotated", help="Directory name for additional HUMAnN3 annotations (if used)"
    )

    # --- 5. Downstream Analysis Options ---
    downstream_group = parser.add_argument_group("Downstream Analysis Options")
    downstream_group.add_argument("--skip-downstream", action="store_true", help="Skip downstream analysis")
    downstream_group.add_argument(
        "--group-col", default="Group", help="The column name to use for grouping in stats (default: 'Group')"
    )

    # --- 6. Differential Abundance Analysis Options ---
    diff_group = parser.add_argument_group("Differential Abundance Options")
    diff_group.add_argument(
        "--run-diff-abundance",
        action="store_true",
        help="Run differential abundance analysis using ANCOM, ALDEx2, and/or ANCOM-BC",
    )
    diff_group.add_argument(
        "--diff-methods",
        default="aldex2,ancom,ancom-bc",
        help="Comma-separated list of methods to use (default: aldex2,ancom,ancom-bc)",
    )
    diff_group.add_argument(
        "--exclude-unmapped", action="store_true", help="Exclude unmapped features from differential abundance analysis"
    )

    # --- 7. Parallel Processing Options ---
    parallel_group = parser.add_argument_group("Parallel Processing Options")
    parallel_group.add_argument(
        "--threads-per-sample", type=int, default=1, help="Number of threads to use per sample"
    )
    parallel_group.add_argument(
        "--max-parallel", type=int, default=None, help="Maximum number of samples to process in parallel"
    )
    parallel_group.add_argument(
        "--use-parallel", action="store_true", help="Use parallel processing for preprocessing steps"
    )

    args = parser.parse_args()

    # Setup logging
    logger = setup_logger(log_file=args.log_file, log_level=getattr(logging, args.log_level))
    log_print("Starting HUMAnN3 Tools Pipeline", level="info")

    start_time = time.time()

    # Handle metadata-driven workflow - collect input files
    input_files = []
    if args.run_preprocessing:
        if args.samples_file:
            # Read from samples file
            from humann3_tools.utils.metadata_utils import read_samples_file
            samples_dict = read_samples_file(args.samples_file)
            for sample_id, files in samples_dict.items():
                input_files.extend(files)
            log_print(f"Loaded {len(input_files)} sequence files from samples file", level='info')
        
        elif args.use_metadata and args.seq_dir:
            # Collect samples from metadata
            from humann3_tools.utils.metadata_utils import collect_samples_from_metadata
            samples_dict = collect_samples_from_metadata(
                metadata_file=args.sample_key,
                seq_dir=args.seq_dir,
                sample_col=args.sample_col,
                r1_col=args.r1_col,
                r2_col=args.r2_col,
                file_pattern=args.file_pattern,
                r1_suffix=args.r1_suffix,
                r2_suffix=args.r2_suffix,
                paired=args.paired
            )
            
            for sample_id, files in samples_dict.items():
                input_files.extend(files)
            
            log_print(f"Collected {len(input_files)} sequence files from {len(samples_dict)} samples", level='info')
        
        # If input files were collected, override args.input_fastq
        if input_files:
            args.input_fastq = input_files
            
        # Now proceed with the original input_fastq check
        if not args.input_fastq:
            log_print("ERROR: No input files specified. Use --input-fastq, --samples-file, or --use-metadata", level='error')
            sys.exit(1)

    # If only listing files, do that and exit
    if args.list_files:
        log_print(f"Files in pathway dir: {args.pathway_dir}", level="info")
        if os.path.isdir(args.pathway_dir):
            for f in sorted(os.listdir(args.pathway_dir)):
                log_print("  " + f, level="info")
        else:
            log_print("  Pathway dir not found.", level="warning")

        log_print(f"Files in gene dir: {args.gene_dir}", level="info")
        if os.path.isdir(args.gene_dir):
            for f in sorted(os.listdir(args.gene_dir)):
                log_print("  " + f, level="info")
        else:
            log_print("  Gene dir not found.", level="warning")

        sys.exit(0)

    # Validate sample key first - we'll need it for both paths
    samples, selected_columns = validate_sample_key(args.sample_key, no_interactive=args.no_interactive)

    ## Preprocessing and HUMAnN3 alignment ###
    preprocessing_results = None
    if args.run_preprocessing:
        # If memory limit is specified, try to apply it
        if args.max_memory:
            success = limit_memory_usage(args.max_memory)
            if success:
                log_print(f"Set memory limit to {args.max_memory} MB", level="info")
            else:
                log_print("Failed to set memory limit, proceeding with unlimited memory", level="warning")

        if not args.input_fastq:
            log_print("ERROR: --input-fastq is required when using --run-preprocessing", level="error")
            sys.exit(1)

        # Check installations
        kneaddata_ok, kneaddata_msg = check_kneaddata_installation()
        if not kneaddata_ok:
            log_print(f"ERROR: KneadData not properly installed: {kneaddata_msg}", level="error")
            sys.exit(1)

        humann3_ok, humann3_msg = check_humann3_installation()
        if not humann3_ok:
            log_print(f"ERROR: HUMAnN3 not properly installed: {humann3_msg}", level="error")
            sys.exit(1)

        # Create preprocessing output directory
        preproc_dir = os.path.join(args.output_dir, "processed_files")
        os.makedirs(preproc_dir, exist_ok=True)

        # Use specified output dirs or create defaults
        kneaddata_output_dir = args.kneaddata_output_dir
        if kneaddata_output_dir is None:
            kneaddata_output_dir = os.path.join(preproc_dir, "kneaddata_output")

        humann3_output_dir = args.humann3_output_dir
        if humann3_output_dir is None:
            humann3_output_dir = os.path.join(preproc_dir, "humann3_output")

        # Choose between regular or parallel processing
        if args.use_parallel:
            log_print("Using parallel preprocessing pipeline", level="info")
            preprocessing_results = run_preprocessing_pipeline_parallel(
                input_files=args.input_fastq,
                output_dir=preproc_dir,
                threads_per_sample=args.threads_per_sample,
                max_parallel=args.max_parallel,
                kneaddata_dbs=args.kneaddata_dbs,
                nucleotide_db=args.humann3_nucleotide_db,
                protein_db=args.humann3_protein_db,
                paired=args.paired,
                kneaddata_output_dir=kneaddata_output_dir,
                humann3_output_dir=humann3_output_dir,
                logger=logger,
            )
        else:
            log_print("Using standard preprocessing pipeline", level="info")
            preprocessing_results = run_preprocessing_pipeline(
                input_files=args.input_fastq,
                output_dir=preproc_dir,
                threads=args.threads,
                kneaddata_dbs=args.kneaddata_dbs,
                nucleotide_db=args.humann3_nucleotide_db,
                protein_db=args.humann3_protein_db,
                paired=args.paired,
                kneaddata_output_dir=kneaddata_output_dir,
                humann3_output_dir=humann3_output_dir,
                logger=logger,
            )



