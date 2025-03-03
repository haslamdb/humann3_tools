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
    preprocessing_group.add_argument(
    "--paired", action="store_true", help="Input files are paired-end reads (default: False)"
)

    preprocessing_group.add_argument(
        "--decontaminate-pairs", default="strict", choices=["strict", "lenient", "unpaired"], 
        help="Method for decontaminating paired-end reads (default: strict)"
    )
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

    if args.use_metadata and args.seq_dir:  
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
        
        # Reset input_files to avoid duplication
        input_files = []
        
        log_print(f"Found {len(samples_dict)} samples with the following files:", level='info')
        for sample_id, files in samples_dict.items():
            file_info = ", ".join([os.path.basename(f) for f in files])
            log_print(f"  Sample {sample_id}: {file_info}", level='debug')
        
        # Process each sample according to paired/unpaired mode
        for sample_id, files in samples_dict.items():
            if args.paired:
                # For paired reads, we need exactly 2 files
                if len(files) == 2:
                    # Add both files for this sample
                    input_files.extend(files)
                    log_print(f"Added paired files for sample {sample_id}", level='debug')
                else:
                    log_print(f"Skipping sample {sample_id}: found {len(files)} files, need exactly 2 for paired mode", level='warning')
            else:
                # For single-end reads, just add the first file
                if files:
                    input_files.append(files[0])
                    log_print(f"Added single file for sample {sample_id}", level='debug')
                else:
                    log_print(f"Skipping sample {sample_id}: no files found", level='warning')
        
        log_print(f"Collected {len(input_files)} sequence files from {len(samples_dict)} samples", level='info')
        
        # If input files were collected, override args.input_fastq
        if input_files:
            args.input_fastq = input_files
            # Log the collected input files
            log_print("Input files for processing:", level='debug')
            for i, file in enumerate(args.input_fastq):
                log_print(f"  {i+1}: {os.path.basename(file)}", level='debug')
        else:
            log_print("ERROR: No input files collected from metadata. Cannot proceed.", level='error')
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

    # Validate sample key first 
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
        # Determine if files are paired or unpaired
        is_paired = args.paired

        kneaddata_options = {}
        if is_paired:
            kneaddata_options["decontaminate-pairs"] = args.decontaminate_pairs

        # For parallel processing, we use threads_per_sample instead of threads
        threads_per_sample = args.threads
        if hasattr(args, 'threads_per_sample') and args.threads_per_sample:
            threads_per_sample = args.threads_per_sample
        
        # Log paired status and file count
        log_print(f"Running in {'paired' if is_paired else 'single-end'} mode with {len(args.input_fastq)} input files", level="debug")
        for i, file in enumerate(args.input_fastq[:min(10, len(args.input_fastq))]):
            log_print(f"  Input file {i+1}: {os.path.basename(file)}", level="debug")
        
        # Run the preprocessing pipeline with updated parameters
        preprocessing_results = run_preprocessing_pipeline_parallel(
            input_files=args.input_fastq,
            output_dir=preproc_dir,
            threads_per_sample=threads_per_sample,  # Use threads_per_sample instead of threads
            max_parallel=args.max_parallel,  # Add this parameter
            kneaddata_dbs=args.kneaddata_dbs,
            nucleotide_db=args.humann3_nucleotide_db,
            protein_db=args.humann3_protein_db,
            paired=is_paired,  # Consistently use this value
            kneaddata_options=kneaddata_options,  
            kneaddata_output_dir=kneaddata_output_dir,
            humann3_output_dir=humann3_output_dir,
            logger=logger,
        )

    if args.run_preprocessing and preprocessing_results:
        log_print("Preprocessing completed successfully. Continuing to HUMAnN3 file processing...", level='info')
        
        # Extract path information from preprocessing results
        humann3_results = preprocessing_results.get('humann3_results', {})
        
        if not humann3_results:
            log_print("ERROR: No HUMAnN3 results found after preprocessing", level='error')
            sys.exit(1)
        
        # Convert HUMAnN3 results to the format expected by the processing functions
        valid_path_samples = []
        valid_gene_samples = []
        
        for sample_id, files in humann3_results.items():
            if files.get('pathabundance'):
                path_file = files['pathabundance']
                valid_path_samples.append((sample_id, path_file))
                log_print(f"Found pathway file for {sample_id}: {os.path.basename(path_file)}", level='debug')
                
            if files.get('genefamilies'):
                gene_file = files['genefamilies']
                valid_gene_samples.append((sample_id, gene_file))
                log_print(f"Found gene file for {sample_id}: {os.path.basename(gene_file)}", level='debug')
        
        log_print(f"Found {len(valid_path_samples)} pathway files and {len(valid_gene_samples)} gene family files", level='info')
        
    # Now we can run the processing steps using these files
        
        # Process pathways
        pathway_unstrat_file = None
        if not args.skip_pathway and valid_path_samples:
            from humann3_tools.humann3.pathway_processing import process_pathway_abundance
            pathway_unstrat_file = process_pathway_abundance(
                valid_path_samples,
                args.humann3_output_dir,  # Use the HUMAnN3 output directory as source
                args.output_dir,
                args.output_prefix,
                selected_columns=selected_columns
            )
        else:
            if args.skip_pathway:
                log_print("Skipping HUMAnN3 pathway processing", level="info")
            else:
                log_print("No valid pathway files; skipping pathway stage", level="warning")

        # Process gene families
        gene_unstrat_file = None
        if not args.skip_gene and valid_gene_samples:
            from humann3_tools.humann3.gene_processing import process_gene_families
            gene_unstrat_file = process_gene_families(
                valid_gene_samples,
                args.humann3_output_dir,  # Use the HUMAnN3 output directory as source
                args.output_dir,
                args.output_prefix,
                selected_columns=selected_columns
            )
        else:
            if args.skip_gene:
                log_print("Skipping HUMAnN3 gene processing", level="info")
            else:
                log_print("No valid gene files; skipping gene stage", level="warning")
        
        # Continue with downstream analysis if not skipped
        if args.skip_downstream:
            log_print("Skipping downstream analysis stage", level="info")
        else:
            if not pathway_unstrat_file and not gene_unstrat_file:
                log_print("No unstratified HUMAnN3 outputs found; cannot run downstream analysis", level="warning")
            else:
                try:
                    # Create analysis output directory
                    downstream_out = os.path.join(args.output_dir, "DownstreamAnalysis")
                    os.makedirs(downstream_out, exist_ok=True)
                    logger.info(f"Downstream analysis output will be in: {downstream_out}")

                    # Read sample metadata
                    from humann3_tools.analysis.metadata import read_and_process_metadata
                    from humann3_tools.utils.file_utils import check_file_exists_with_logger
                    
                    if not check_file_exists_with_logger(args.sample_key, "Sample key", logger):
                        log_print("Cannot proceed with downstream analysis; missing sample key", level="error")
                    else:
                        sample_key_df = read_and_process_metadata(args.sample_key, logger)

                        # Process gene families if available
                        from humann3_tools.analysis.visualizations import read_and_process_gene_families
                        if gene_unstrat_file and check_file_exists_with_logger(gene_unstrat_file, "Gene families", logger):
                            read_and_process_gene_families(gene_unstrat_file, sample_key_df, downstream_out, logger)

                        # Process pathways if available
                        from humann3_tools.analysis.visualizations import read_and_process_pathways
                        from humann3_tools.analysis.statistical import run_statistical_tests
                        
                        if pathway_unstrat_file and check_file_exists_with_logger(pathway_unstrat_file, "Pathways", logger):
                            pathways_merged = read_and_process_pathways(
                                pathway_unstrat_file,
                                sample_key_df,
                                downstream_out,
                                logger
                            )
                            # Run statistical tests
                            run_statistical_tests(pathways_merged, downstream_out, logger, group_col=args.group_col)
                except Exception as e:
                    logger.error(f"Downstream analysis failed: {e}")
                    logger.error(traceback.format_exc())
        
        # Skip the normal processing path since we've already done it
        log_print("Full pipeline completed successfully", level='info')
        elapsed = time.time() - start_time
        hh, rr = divmod(elapsed, 3600)
        mm, ss = divmod(rr, 60)
        log_print(f"Pipeline finished in {int(hh)}h {int(mm)}m {int(ss)}s", level="info")
        sys.exit(0)
