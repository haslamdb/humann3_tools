# humann3_tools/preprocessing/pipeline.py
import os
import logging
from humann3_tools.preprocessing.kneaddata import run_kneaddata, check_kneaddata_installation, run_kneaddata_parallel
from humann3_tools.preprocessing.humann3_run import run_humann3, check_humann3_installation, run_humann3_parallel
from humann3_tools.logger import log_print
from humann3_tools.utils.resource_utils import (
    track_peak_memory, 
    monitor_memory_usage, 
    stop_memory_monitoring
)

def run_preprocessing_pipeline(input_files, output_dir, threads=1, 
                              kneaddata_dbs=None, nucleotide_db=None, protein_db=None,
                              kneaddata_options=None, humann3_options=None, 
                              paired=False, kneaddata_output_dir=None, humann3_output_dir=None,
                              logger=None):
    """
    Run the full preprocessing pipeline: KneadData → HUMAnN3.
    
    Args:
        input_files: List of input FASTQ files
        output_dir: Base directory for outputs
        threads: Number of threads to use
        kneaddata_dbs: Path to KneadData reference database(s)
        nucleotide_db: Path to HUMAnN3 nucleotide database
        protein_db: Path to HUMAnN3 protein database
        kneaddata_options: Dict of additional KneadData options
        humann3_options: Dict of additional HUMAnN3 options
        paired: Whether input files are paired
        kneaddata_output_dir: Custom directory for KneadData outputs
        humann3_output_dir: Custom directory for HUMAnN3 outputs
        logger: Logger instance
    """
    # Use specified output directories or create defaults
    if kneaddata_output_dir is None:
        kneaddata_output = os.path.join(output_dir, "kneaddata_output")
    else:
        kneaddata_output = kneaddata_output_dir
        
    if humann3_output_dir is None:
        humann3_output = os.path.join(output_dir, "humann3_output")
    else:
        humann3_output = humann3_output_dir
        
    os.makedirs(kneaddata_output, exist_ok=True)
    os.makedirs(humann3_output, exist_ok=True)
    
    if logger is None:
        logger = logging.getLogger('humann3_analysis')
    
    # Check installations
    kneaddata_ok, kneaddata_version = check_kneaddata_installation()
    if not kneaddata_ok:
        logger.error(f"KneadData not properly installed: {kneaddata_version}")
        return None
    
    humann3_ok, humann3_version = check_humann3_installation()
    if not humann3_ok:
        logger.error(f"HUMAnN3 not properly installed: {humann3_version}")
        return None
    
    logger.info(f"Starting preprocessing pipeline with {len(input_files)} input files")
    logger.info(f"KneadData version: {kneaddata_version}")
    logger.info(f"HUMAnN3 version: {humann3_version}")
    logger.info(f"KneadData output dir: {kneaddata_output}")
    logger.info(f"HUMAnN3 output dir: {humann3_output}")
    
    # Step 1: Run KneadData
    logger.info("Starting KneadData step...")
    kneaddata_files = run_kneaddata(
        input_files=input_files,
        output_dir=kneaddata_output,
        threads=threads,
        reference_dbs=kneaddata_dbs, 
        paired=paired,
        additional_options=kneaddata_options,
        logger=logger
    )
    
    if not kneaddata_files:
        logger.error("KneadData step failed, stopping pipeline")
        return None
    
    logger.info(f"KneadData completed successfully with {len(kneaddata_files)} output files")
    
    # Filter out contaminant files
    kneaddata_files = [f for f in kneaddata_files if "contam" not in os.path.basename(f).lower()]
    logger.info(f"After filtering, {len(kneaddata_files)} files will be used for HUMAnN3")

        
    # Concatenate paired files for HUMAnN3
    logger.info("Preparing KneadData outputs for HUMAnN3...")
    paired_files = {}

    # First, check that we have valid files to work with
    if not kneaddata_files:
        logger.error("No valid KneadData output files to prepare for HUMAnN3")
        return None

    for file in kneaddata_files:
        # Log file information for debugging
        file_size = os.path.getsize(file)
        logger.info(f"Processing KneadData output file: {os.path.basename(file)} (Size: {file_size} bytes)")
        
        if file_size == 0:
            logger.warning(f"Skipping empty file: {file}")
            continue
            
        # Extract the base sample name - handle different KneadData output naming patterns
        filename = os.path.basename(file)
        
        # Try to extract the sample name from different possible patterns
        base_name = None
        if "_paired_1.fastq" in filename:
            base_name = filename.split("_paired_1.fastq")[0]
        elif "_paired_2.fastq" in filename:
            base_name = filename.split("_paired_2.fastq")[0]
        elif "_clean.fastq" in filename:
            base_name = filename.split("_clean.fastq")[0]
        else:
            # Just use the filename without extension as a fallback
            base_name = os.path.splitext(filename)[0]
        
        # Group files by sample name
        if base_name not in paired_files:
            paired_files[base_name] = []
        paired_files[base_name].append(file)

    # Prepare files for HUMAnN3
    humann3_input_files = []

    # Process each sample's files
    for sample, files in paired_files.items():
        logger.info(f"Sample {sample} has {len(files)} KneadData output files")
        
        # Check if we need to concatenate files
        if len(files) > 1 and paired:
            # Sort files to ensure consistent order (R1 before R2)
            files.sort()
            
            # Create concatenated file in the same directory
            concatenated_file = os.path.join(os.path.dirname(files[0]), f"{sample}_paired_concat.fastq")
            
            logger.info(f"Concatenating {len(files)} files for sample {sample} to {os.path.basename(concatenated_file)}")
            
            # Concatenate the files
            with open(concatenated_file, 'w') as outfile:
                total_lines = 0
                for file in files:
                    logger.debug(f"  Adding file: {os.path.basename(file)} (Size: {os.path.getsize(file)} bytes)")
                    
                    try:
                        with open(file, 'r') as infile:
                            file_content = infile.read()
                            outfile.write(file_content)
                            line_count = len(file_content.splitlines())
                            total_lines += line_count
                            logger.debug(f"  Added {line_count} lines from {os.path.basename(file)}")
                    except Exception as e:
                        logger.error(f"Error reading/writing file {file}: {str(e)}")
            
            # Verify the concatenated file was created properly
            if os.path.exists(concatenated_file) and os.path.getsize(concatenated_file) > 0:
                logger.info(f"Successfully created concatenated file: {os.path.basename(concatenated_file)} "
                        f"(Size: {os.path.getsize(concatenated_file)} bytes, Lines: {total_lines})")
                humann3_input_files.append(concatenated_file)
            else:
                logger.error(f"Failed to create valid concatenated file for {sample}. "
                        f"File exists: {os.path.exists(concatenated_file)}, "
                        f"Size: {os.path.getsize(concatenated_file) if os.path.exists(concatenated_file) else 'N/A'}")
                # Fall back to using individual files
                logger.info(f"Falling back to using individual files for sample {sample}")
                humann3_input_files.extend(files)
        else:
            # For single files or non-paired mode, use them directly
            logger.info(f"Using {len(files)} file(s) directly for sample {sample}")
            humann3_input_files.extend(files)

    logger.info(f"Prepared {len(humann3_input_files)} input files for HUMAnN3")

    # Now check that we have valid files to run HUMAnN3 on
    if not humann3_input_files:
        logger.error("No valid input files prepared for HUMAnN3")
        return None

    # Verify all input files exist and have content
    for file in humann3_input_files:
        if not os.path.exists(file):
            logger.error(f"HUMAnN3 input file does not exist: {file}")
            return None
        
        if os.path.getsize(file) == 0:
            logger.error(f"HUMAnN3 input file is empty: {file}")
            return None
        
        logger.info(f"Validated HUMAnN3 input file: {os.path.basename(file)} (Size: {os.path.getsize(file)} bytes)")
    
    # Step 2: Run HUMAnN3
    logger.info("Starting HUMAnN3 step...")
    humann3_results = run_humann3(
        input_files=humann3_input_files, 
        output_dir=humann3_output,
        threads=threads,
        nucleotide_db=nucleotide_db,
        protein_db=protein_db,
        additional_options=humann3_options,
        logger=logger
    )
    if not humann3_results:
        logger.error("HUMAnN3 step failed")
        return None
    
    logger.info(f"HUMAnN3 completed successfully for {len(humann3_results)} samples")
    
    # Return combined results
    return {
        'kneaddata_files': kneaddata_files,
        'humann3_results': humann3_results
    }

def run_preprocessing_pipeline_parallel(input_files, output_dir, threads_per_sample=1, 
                                       max_parallel=None, kneaddata_dbs=None, 
                                       nucleotide_db=None, protein_db=None,
                                       kneaddata_options=None, humann3_options=None, 
                                       paired=False, kneaddata_output_dir=None, humann3_output_dir=None,
                                       logger=None):
    """
    Run the full preprocessing pipeline in parallel: KneadData → HUMAnN3.
    
    Args:
        input_files: List of input FASTQ files
        output_dir: Base directory for outputs
        threads_per_sample: Number of threads per sample
        max_parallel: Maximum number of parallel samples (None = CPU count)
        kneaddata_dbs: Path(s) to KneadData reference database(s). Can be a string or a list of paths.
        nucleotide_db: Path to HUMAnN3 nucleotide database
        protein_db: Path to HUMAnN3 protein database
        kneaddata_options: Dict of additional KneadData options
        humann3_options: Dict of additional HUMAnN3 options
        paired: Whether input files are paired
        kneaddata_output_dir: Custom directory for KneadData outputs
        humann3_output_dir: Custom directory for HUMAnN3 outputs
        logger: Logger instance
        
    Returns:
        Dict of final HUMAnN3 output file paths by sample and type
    """
    if logger is None:
        logger = logging.getLogger('humann3_analysis')
    
    # Check installations
    kneaddata_ok, kneaddata_version = check_kneaddata_installation()
    if not kneaddata_ok:
        logger.error(f"KneadData not properly installed: {kneaddata_version}")
        return None
    
    humann3_ok, humann3_version = check_humann3_installation()
    if not humann3_ok:
        logger.error(f"HUMAnN3 not properly installed: {humann3_version}")
        return None
    
    logger.info(f"Starting parallel preprocessing pipeline with {len(input_files)} input files")
    logger.info(f"KneadData version: {kneaddata_version}")
    logger.info(f"HUMAnN3 version: {humann3_version}")
    
    # Use specified output directories or create defaults
    if kneaddata_output_dir is None:
        kneaddata_output = os.path.join(output_dir, "kneaddata_output")
    else:
        kneaddata_output = kneaddata_output_dir
        
    if humann3_output_dir is None:
        humann3_output = os.path.join(output_dir, "humann3_output")
    else:
        humann3_output = humann3_output_dir
    
    logger.info(f"KneadData output dir: {kneaddata_output}")
    logger.info(f"HUMAnN3 output dir: {humann3_output}")
    
    os.makedirs(kneaddata_output, exist_ok=True)
    os.makedirs(humann3_output, exist_ok=True)
    
    # Step 1: Run KneadData in parallel
    logger.info("Starting KneadData step in parallel...")
    kneaddata_results = run_kneaddata_parallel(
        input_files=input_files,
        output_dir=kneaddata_output,
        threads=threads_per_sample,
        max_parallel=max_parallel,
        reference_dbs=kneaddata_dbs,
        paired=paired,
        additional_options=kneaddata_options,
        logger=logger
    )
    
    if not kneaddata_results:
        logger.error("KneadData step failed, stopping pipeline")
        return None
    
    # Flatten the list of KneadData output files
    kneaddata_files = []
    for sample_id, files in kneaddata_results.items():
        if isinstance(files, list):
            kneaddata_files.extend(files)
            
    # Filter out contaminant files
    kneaddata_files = [f for f in kneaddata_files if "contam" not in os.path.basename(f).lower()]
    logger.info(f"After filtering, {len(kneaddata_files)} files will be used for HUMAnN3")
    
    logger.info(f"KneadData completed with {len(kneaddata_files)} output files")
    
    # Step 2: Run HUMAnN3 in parallel
    logger.info("Starting HUMAnN3 step in parallel...")
    humann3_results = run_humann3_parallel(
        input_files=kneaddata_files,
        output_dir=humann3_output,
        threads=threads_per_sample,
        max_parallel=max_parallel,
        nucleotide_db=nucleotide_db,
        protein_db=protein_db,
        additional_options=humann3_options,
        logger=logger
    )
    
    if not humann3_results:
        logger.error("HUMAnN3 step failed")
        return None
    
    logger.info(f"HUMAnN3 completed for {len(humann3_results)} samples")
    
    # Return combined results
    return {
        'kneaddata_files': kneaddata_files,
        'humann3_results': humann3_results
    }

# @track_peak_memory
# def run_preprocessing_pipeline_parallel(input_files, output_dir, threads_per_sample=1, 
#                                        max_parallel=None, kneaddata_dbs=None, 
#                                        nucleotide_db=None, protein_db=None,
#                                        kneaddata_options=None, humann3_options=None, 
#                                        paired=False, kneaddata_output_dir=None, humann3_output_dir=None,
#                                        logger=None):
#     """
#     Run the full preprocessing pipeline in parallel: KneadData → HUMAnN3.
    
#     Args:
#         input_files: List of input FASTQ files
#         output_dir: Base directory for outputs
#         threads_per_sample: Number of threads per sample
#         max_parallel: Maximum number of parallel samples (None = CPU count)
#         kneaddata_dbs: Path(s) to KneadData reference database(s). Can be a string or a list of paths.
#         nucleotide_db: Path to HUMAnN3 nucleotide database
#         protein_db: Path to HUMAnN3 protein database
#         kneaddata_options: Dict of additional KneadData options
#         humann3_options: Dict of additional HUMAnN3 options
#         paired: Whether input files are paired
#         kneaddata_output_dir: Custom directory for KneadData outputs
#         humann3_output_dir: Custom directory for HUMAnN3 outputs
#         logger: Logger instance
        
#     Returns:
#         Dict of final HUMAnN3 output file paths by sample and type
#     """
#     if logger is None:
#         logger = logging.getLogger('humann3_analysis')
    
#     # Check installations
#     kneaddata_ok, kneaddata_version = check_kneaddata_installation()
#     if not kneaddata_ok:
#         logger.error(f"KneadData not properly installed: {kneaddata_version}")
#         return None
    
#     humann3_ok, humann3_version = check_humann3_installation()
#     if not humann3_ok:
#         logger.error(f"HUMAnN3 not properly installed: {humann3_version}")
#         return None
    
#     logger.info(f"Starting parallel preprocessing pipeline with {len(input_files)} input files")
#     logger.info(f"KneadData version: {kneaddata_version}")
#     logger.info(f"HUMAnN3 version: {humann3_version}")
    
#     # Use specified output directories or create defaults
#     if kneaddata_output_dir is None:
#         kneaddata_output = os.path.join(output_dir, "kneaddata_output")
#     else:
#         kneaddata_output = kneaddata_output_dir
        
#     if humann3_output_dir is None:
#         humann3_output = os.path.join(output_dir, "humann3_output")
#     else:
#         humann3_output = humann3_output_dir
    
#     logger.info(f"KneadData output dir: {kneaddata_output}")
#     logger.info(f"HUMAnN3 output dir: {humann3_output}")
    
#     os.makedirs(kneaddata_output, exist_ok=True)
#     os.makedirs(humann3_output, exist_ok=True)
    
#     # Step 1: Run KneadData in parallel
#     logger.info("Starting KneadData step in parallel...")
#     kneaddata_results = run_kneaddata_parallel(
#         input_files=input_files,
#         output_dir=kneaddata_output,
#         threads=threads_per_sample,
#         max_parallel=max_parallel,
#         reference_dbs=kneaddata_dbs,
#         paired=paired,
#         additional_options=kneaddata_options,
#         logger=logger
#     )
    
#     if not kneaddata_results:
#         logger.error("KneadData step failed, stopping pipeline")
#         return None
    
#     # Flatten the list of KneadData output files
#     kneaddata_files = []
#     for sample_id, files in kneaddata_results.items():
#         if isinstance(files, list):
#             kneaddata_files.extend(files)
            
#     # Filter out contaminant files
#     kneaddata_files = [f for f in kneaddata_files if "contam" not in os.path.basename(f).lower()]
#     logger.info(f"After filtering, {len(kneaddata_files)} files will be used for HUMAnN3")
    
#     logger.info(f"KneadData completed with {len(kneaddata_files)} output files")
    
#     # Step 2: Run HUMAnN3 in parallel
#     logger.info("Starting HUMAnN3 step in parallel...")
#     humann3_results = run_humann3_parallel(
#         input_files=kneaddata_files,
#         output_dir=humann3_output,
#         threads=threads_per_sample,
#         max_parallel=max_parallel,
#         nucleotide_db=nucleotide_db,
#         protein_db=protein_db,
#         additional_options=humann3_options,
#         logger=logger
#     )
    
#     if not humann3_results:
#         logger.error("HUMAnN3 step failed")
#         return None
    
#     logger.info(f"HUMAnN3 completed for {len(humann3_results)} samples")
    
#     # Return combined results
#     return {
#         'kneaddata_files': kneaddata_files,
#         'humann3_results': humann3_results
#     }