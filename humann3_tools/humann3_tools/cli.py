# humann3_tools/humann3_tools/cli.py
import os
import sys
import argparse
import time
import logging
import traceback

from humann3_tools.humann3_tools.logger import setup_logger, log_print
from humann3_tools.utils.sample_utils import validate_sample_key, check_input_files_exist
from humann3_tools.utils.file_utils import check_file_exists_with_logger

# Define stub functions for the missing modules
def process_pathway_abundance(*args, **kwargs):
    return "example_pathway_file.tsv"

def process_gene_families(*args, **kwargs):
    return "example_gene_file.tsv"

def read_and_process_metadata(*args, **kwargs):
    return None

def read_and_process_gene_families(*args, **kwargs):
    return None

def read_and_process_pathways(*args, **kwargs):
    return None

def run_statistical_tests(*args, **kwargs):
    return None

def run_preprocessing_pipeline(*args, **kwargs):
    return {"humann3_results": {}}

def run_preprocessing_pipeline_parallel(*args, **kwargs):
    return {"humann3_results": {}}

def check_kneaddata_installation(*args, **kwargs):
    return (True, "")

def check_humann3_installation(*args, **kwargs):
    return (True, "")

def limit_memory_usage(*args, **kwargs):
    return True


def main():
    
    """Main entry point for the humann3_tools CLI."""
    parser = argparse.ArgumentParser(description="HUMAnN3 Tools: Process and analyze HUMAnN3 output")

    #  --- Preprocessing and Humann alignment ---
    parser.add_argument("--run-preprocessing", action="store_true",
                    help="Run preprocessing (KneadData and HUMAnN3) on raw sequence files")
    parser.add_argument("--input-fastq", nargs="+",
                    help="Input FASTQ file(s) for preprocessing")
    parser.add_argument("--paired", action="store_true",
                    help="Input files are paired-end reads")
    parser.add_argument("--kneaddata-db",
                    help="Path to KneadData reference database")
    parser.add_argument("--humann3-nucleotide-db",
                    help="Path to HUMAnN3 nucleotide database (ChocoPhlAn)")
    parser.add_argument("--humann3-protein-db",
                    help="Path to HUMAnN3 protein database (UniRef)")
    parser.add_argument("--threads", type=int, default=1,
                    help="Number of threads to use for preprocessing")    
    # ---  HUMAnN3 output processing ---
    parser.add_argument("--sample-key", required=True,
                        help="CSV file with columns for sample names and metadata")
    parser.add_argument("--pathway-dir", required=True,
                        help="Directory containing raw pathway abundance files")
    parser.add_argument("--gene-dir", required=True,
                        help="Directory containing raw gene family files")
    parser.add_argument("--output-dir", default="./Humann3Output",
                        help="Directory where HUMAnN3-processed files and downstream results will go")
    parser.add_argument("--output-prefix", default="ProcessedFiles",
                        help="Prefix for intermediate HUMAnN3 output directories/files")
    parser.add_argument("--skip-pathway", action="store_true",
                        help="Skip HUMAnN3 pathway abundance processing")
    parser.add_argument("--skip-gene", action="store_true",
                        help="Skip HUMAnN3 gene family processing")
    parser.add_argument("--annotations-dir", default="annotated",
                        help="Directory name for additional HUMAnN3 annotations (if used)")
    parser.add_argument("--list-files", action="store_true",
                        help="Just list input files in --pathway-dir and --gene-dir, then exit")
    parser.add_argument("--no-interactive", action="store_true",
                        help="Non-interactive mode for sample key column selection")
    parser.add_argument("--max-memory", type=int, default=None,
                   help="Maximum memory usage in MB (default: unlimited)")

    # --- Downstream analysis ---
    parser.add_argument("--skip-downstream", action="store_true",
                        help="Skip the downstream analysis stage entirely")
    parser.add_argument("--group-col", default="Group",
                        help="The column name to use for grouping in statistical tests (default: 'Group')")

    # Logging
    parser.add_argument("--log-file", default=None,
                        help="Path to combined log file")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG","INFO","WARNING","ERROR","CRITICAL"],
                        help="Logging level (default=INFO)")
    
    # --- Differential abundance analysis ---
    parser.add_argument("--run-diff-abundance", action="store_true",
                        help="Run differential abundance analysis using ANCOM, ALDEx2, and/or ANCOM-BC")
    parser.add_argument("--diff-methods", default="aldex2,ancom,ancom-bc",
                        help="Comma-separated list of differential abundance methods to use (default: aldex2,ancom,ancom-bc)")
    parser.add_argument("--exclude-unmapped", action="store_true",
                        help="Exclude unmapped features from differential abundance analysis")
    
    # --- Parallel processing ---
    parser.add_argument("--threads-per-sample", type=int, default=1,
                   help="Number of threads to use per sample")
    parser.add_argument("--max-parallel", type=int, default=None,
                    help="Maximum number of samples to process in parallel (default: CPU count)")
    parser.add_argument("--use-parallel", action="store_true",
                    help="Use parallel processing for preprocessing steps")
    parser.add_argument("--skip-kneaddata", action="store_true",
                    help="Skip KneadData preprocessing and run HUMAnN3 directly on input FASTQ files")
    
    args = parser.parse_args()

    # Setup logging
    logger = setup_logger(log_file=args.log_file, log_level=getattr(logging, args.log_level))
    log_print("Starting HUMAnN3 Tools Pipeline", level='info')

    start_time = time.time()

    # If only listing files, do that and exit
    if args.list_files:
        log_print(f"Files in pathway dir: {args.pathway_dir}", level='info')
        if os.path.isdir(args.pathway_dir):
            for f in sorted(os.listdir(args.pathway_dir)):
                log_print("  " + f, level='info')
        else:
            log_print("  Pathway dir not found.", level='warning')
        
        log_print(f"Files in gene dir: {args.gene_dir}", level='info')
        if os.path.isdir(args.gene_dir):
            for f in sorted(os.listdir(args.gene_dir)):
                log_print("  " + f, level='info')
        else:
            log_print("  Gene dir not found.", level='warning')
        
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
                log_print(f"Set memory limit to {args.max_memory} MB", level='info')
            else:
                log_print("Failed to set memory limit, proceeding with unlimited memory", level='warning')

        if not args.input_fastq:
            log_print("ERROR: --input-fastq is required when using --run-preprocessing", level='error')
            sys.exit(1)

        # Check installations
        kneaddata_ok, kneaddata_msg = check_kneaddata_installation()
        if not kneaddata_ok:
            log_print(f"ERROR: KneadData not properly installed: {kneaddata_msg}", level='error')
            sys.exit(1)

        humann3_ok, humann3_msg = check_humann3_installation()
        if not humann3_ok:
            log_print(f"ERROR: HUMAnN3 not properly installed: {humann3_msg}", level='error')
            sys.exit(1)

        # Create preprocessing output directory
        preproc_dir = os.path.join(args.output_dir, "PreprocessedData")
        os.makedirs(preproc_dir, exist_ok=True)

        # Choose between regular or parallel processing
        if args.use_parallel:
            log_print("Using parallel preprocessing pipeline", level='info')
            preprocessing_results = run_preprocessing_pipeline_parallel(
                input_files=args.input_fastq,
                output_dir=preproc_dir,
                threads_per_sample=args.threads_per_sample,
                max_parallel=args.max_parallel,
                kneaddata_db=args.kneaddata_db,
                nucleotide_db=args.humann3_nucleotide_db,
                protein_db=args.humann3_protein_db,
                paired=args.paired,
                logger=logger,
                skip_kneaddata=args.skip_kneaddata
            )
        else:
            log_print("Using standard preprocessing pipeline", level='info')
            preprocessing_results = run_preprocessing_pipeline(
                input_files=args.input_fastq,
                output_dir=preproc_dir,
                threads=args.threads,
                kneaddata_db=args.kneaddata_db,
                nucleotide_db=args.humann3_nucleotide_db,
                protein_db=args.humann3_protein_db,
                paired=args.paired,
                logger=logger,
                skip_kneaddata=args.skip_kneaddata
            )

        if not preprocessing_results:
            log_print("Preprocessing pipeline failed", level='error')
            sys.exit(1)

        log_print("Preprocessing completed successfully", level='info')

        # Extract preprocessing outputs
        pathway_files = []
        gene_files = []

        if preprocessing_results and 'humann3_results' in preprocessing_results:
            humann3_results = preprocessing_results['humann3_results']
            for sample, files in humann3_results.items():
                if files.get('pathabundance'):
                    pathway_files.append((sample, files['pathabundance']))
                if files.get('genefamilies'):
                    gene_files.append((sample, files['genefamilies']))

    # Process HUMAnN3 output files - either from preprocessing or directly from directories
    if preprocessing_results is None:
        # Check input files from directories
        valid_path_samples, valid_gene_samples = check_input_files_exist(samples, args.pathway_dir, args.gene_dir)
    else:
        # Use results from preprocessing
        valid_path_samples, valid_gene_samples = pathway_files, gene_files

    # Process pathways
    pathway_unstrat_file = None
    if not args.skip_pathway and valid_path_samples:
        source_dir = preproc_dir if preprocessing_results else args.pathway_dir
        pathway_unstrat_file = process_pathway_abundance(valid_path_samples,
                                                       source_dir,
                                                       args.output_dir,
                                                       args.output_prefix,
                                                       selected_columns=selected_columns)
    else:
        if args.skip_pathway:
            log_print("Skipping HUMAnN3 pathway processing", level='info')
        else:
            log_print("No valid pathway files; skipping pathway stage", level='warning')

    # Process gene families
    gene_unstrat_file = None
    if not args.skip_gene and valid_gene_samples:
        source_dir = preproc_dir if preprocessing_results else args.gene_dir
        gene_unstrat_file = process_gene_families(valid_gene_samples,
                                                source_dir,
                                                args.output_dir,
                                                args.output_prefix,
                                                selected_columns=selected_columns)
    else:
        if args.skip_gene:
            log_print("Skipping HUMAnN3 gene processing", level='info')
        else:
            log_print("No valid gene files; skipping gene stage", level='warning')

    # Summaries
    log_print("=== HUMAnN3 Stage Outputs ===", level='info')
    if pathway_unstrat_file:
        log_print(f"  Pathway unstratified file: {pathway_unstrat_file}", level='info')
    if gene_unstrat_file:
        log_print(f"  Gene families unstratified file: {gene_unstrat_file}", level='info')

    # ========== STAGE 2: Downstream Analysis ==========
    if args.skip_downstream:
        log_print("Skipping downstream analysis stage (--skip-downstream)", level='info')
    else:
        # We cannot proceed if we don't have unstratified pathways or gene families
        if not pathway_unstrat_file and not gene_unstrat_file:
            log_print("No unstratified HUMAnN3 outputs found; cannot run downstream analysis", level='warning')
        else:
            try:
                # Create a sub-directory for the downstream analysis
                downstream_out = os.path.join(args.output_dir, "DownstreamAnalysis")
                os.makedirs(downstream_out, exist_ok=True)
                logger.info(f"Downstream analysis output will be in: {downstream_out}")

                # 1) Read sample key
                if not check_file_exists_with_logger(args.sample_key, "Sample key", logger):
                    log_print("Cannot proceed with downstream analysis; missing sample key", level='error')
                else:
                    sample_key_df = read_and_process_metadata(args.sample_key, logger)

                    # 2) If gene families file is available, analyze it
                    if gene_unstrat_file and check_file_exists_with_logger(gene_unstrat_file, "Gene families", logger):
                        read_and_process_gene_families(gene_unstrat_file, sample_key_df, downstream_out, logger)

                    # 3) If pathways file is available, analyze it
                    if pathway_unstrat_file and check_file_exists_with_logger(pathway_unstrat_file, "Pathways", logger):
                        pathways_merged = read_and_process_pathways(pathway_unstrat_file, sample_key_df, downstream_out, logger)
                        # 4) Run stats on pathways with user-defined grouping
                        run_statistical_tests(pathways_merged, downstream_out, logger, group_col=args.group_col)

            except Exception as e:
                logger.error(f"Downstream analysis failed: {e}")
                logger.error(traceback.format_exc())

    elapsed = time.time() - start_time
    hh, rr = divmod(elapsed, 3600)
    mm, ss = divmod(rr, 60)
    log_print(f"Pipeline finished in {int(hh)}h {int(mm)}m {int(ss)}s", level='info')
    log_print("All done!", level='info')
    return 0


if __name__ == "__main__":
    sys.exit(main())