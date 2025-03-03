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
        kneaddata_dbs: Path to KneadData reference database(s). Can be a string or a list of paths.
        nucleotide_db: Path to HUMAnN3 nucleotide database
        protein_db: Path to HUMAnN3 protein database
        kneaddata_options: Dict of additional KneadData options
        humann3_options: Dict of additional HUMAnN3 options
        paired: Whether input files are paired
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
    
    logger.info(f"Starting preprocessing pipeline with {len(input_files)} input files")
    logger.info(f"KneadData version: {kneaddata_version}")
    logger.info(f"HUMAnN3 version: {humann3_version}")
    
    # Create output directories - use provided paths if available, else use defaults
    if kneaddata_output_dir is None:
        kneaddata_output = os.path.join(output_dir, "kneaddata_output")
    else:
        kneaddata_output = kneaddata_output_dir
        
    if humann3_output_dir is None:
        humann3_output = os.path.join(output_dir, "humann3_output")
    else:
        humann3_output = humann3_output_dir

    # Ensure the base output directory exists
    os.makedirs(output_dir, exist_ok=True)
    # Create the specific output directories
    os.makedirs(kneaddata_output, exist_ok=True)
    os.makedirs(humann3_output, exist_ok=True)

    # Log the directories being used
    logger.info(f"Using KneadData output directory: {kneaddata_output}")
    logger.info(f"Using HUMAnN3 output directory: {humann3_output}")
    
# Step 1: Run KneadData
    logger.info(f"Starting KneadData step (paired={paired})...")
    
    # For paired reads, we need to group input files properly
    if paired:
        # Ensure we have an even number of files for paired-end mode
        if len(input_files) % 2 != 0:
            logger.error("Paired mode requires an even number of input files")
            return None
            
        # Process files in pairs
        kneaddata_files = []
        for i in range(0, len(input_files), 2):
            if i+1 < len(input_files):
                # Run KneadData on this pair
                pair_files = [input_files[i], input_files[i+1]]
                sample_name = os.path.basename(pair_files[0]).split('_')[0]  # Get basename for logging
                
                logger.info(f"Processing paired files for sample {sample_name}")
                
                # Create sample-specific output directory
                sample_output = os.path.join(kneaddata_output, sample_name)
                os.makedirs(sample_output, exist_ok=True)
                
                pair_results = run_kneaddata(
                    input_files=pair_files,
                    output_dir=sample_output,
                    threads=threads,
                    reference_dbs=kneaddata_dbs,
                    paired=True,  # Always True in this block
                    additional_options=kneaddata_options,
                    logger=logger
                )
                
                if pair_results:
                    kneaddata_files.extend(pair_results)
                else:
                    logger.warning(f"No KneadData output files for sample {sample_name}")
    else:
        # Single-end mode: process each file individually
        kneaddata_files = []
        for file in input_files:
            sample_name = os.path.basename(file).split('_')[0]  # Get basename for logging
            
            logger.info(f"Processing single-end file for sample {sample_name}")
            
            # Create sample-specific output directory
            sample_output = os.path.join(kneaddata_output, sample_name)
            os.makedirs(sample_output, exist_ok=True)
            
            file_results = run_kneaddata(
                input_files=[file],
                output_dir=sample_output,
                threads=threads,
                reference_dbs=kneaddata_dbs,
                paired=False,  # Always False in this block
                additional_options=kneaddata_options,
                logger=logger
            )
            
            if file_results:
                kneaddata_files.extend(file_results)
            else:
                logger.warning(f"No KneadData output files for sample {sample_name}")

    
    if not kneaddata_files:
        logger.error("KneadData step failed, stopping pipeline")
        return None
    
    logger.info(f"KneadData completed successfully with {len(kneaddata_files)} output files")

    # After KneadData has run and produced files:
    logger.info("Preparing KneadData outputs for HUMAnN3...")
    sample_files = {}
    humann3_input_files = []

    # First, check that we have valid files to work with
    if not kneaddata_files:
        logger.error("No valid KneadData output files to prepare for HUMAnN3")
        return None

    # Log all found files for debugging
    logger.info(f"Found {len(kneaddata_files)} KneadData output files:")
    for i, file in enumerate(kneaddata_files[:10]):  # Show first 10 to avoid log spam
        logger.info(f"  {i+1}. {os.path.basename(file)}")
    if len(kneaddata_files) > 10:
        logger.info(f"  ... and {len(kneaddata_files) - 10} more files")

    # Extract sample names more carefully
    for file in kneaddata_files:
        filename = os.path.basename(file)
        
        # Define reliable patterns for extracting sample names
        sample_name = None
        
        # Pattern 1: Extract from paired files
        if "_paired_1.fastq" in filename:
            parts = filename.split("_paired_1.fastq")[0]
            sample_name = parts.split("_kneaddata")[0] if "_kneaddata" in parts else parts
        elif "_paired_2.fastq" in filename:
            parts = filename.split("_paired_2.fastq")[0]
            sample_name = parts.split("_kneaddata")[0] if "_kneaddata" in parts else parts
        
        # Skip files we can't categorize
        if not sample_name:
            logger.debug(f"Skipping file with unclear sample name: {filename}")
            continue
            
        # Initialize list for this sample if not exists
        if sample_name not in sample_files:
            sample_files[sample_name] = []
        
        # Add file to sample's list (only if not already added)
        if file not in sample_files[sample_name]:
            sample_files[sample_name].append(file)
            logger.debug(f"Added {filename} to sample {sample_name}")

    # Process each sample's paired files
    logger.info(f"Found {len(sample_files)} samples after organizing files")
    
    for sample, files in sample_files.items():
        logger.info(f"Processing sample {sample} with {len(files)} KneadData output files")
        
        # Group files by paired number
        paired_1_files = [f for f in files if "_paired_1.fastq" in os.path.basename(f)]
        paired_2_files = [f for f in files if "_paired_2.fastq" in os.path.basename(f)]
        
        # Handle the case where we have duplicate files
        if len(paired_1_files) > 1:
            logger.warning(f"Found {len(paired_1_files)} R1 files for sample {sample}, using the first one")
            paired_1_files = [paired_1_files[0]]
            
        if len(paired_2_files) > 1:
            logger.warning(f"Found {len(paired_2_files)} R2 files for sample {sample}, using the first one")
            paired_2_files = [paired_2_files[0]]
        
        # Check if we have a proper pair
        if len(paired_1_files) == 1 and len(paired_2_files) == 1:
            # Create concatenated file
            concatenated_file = os.path.join(os.path.dirname(paired_1_files[0]), f"{sample}_paired_concat.fastq")
            logger.info(f"Concatenating paired files for sample {sample} to {os.path.basename(concatenated_file)}")
            
            try:
                # Create concatenated file
                with open(concatenated_file, 'w') as outfile:
                    for file in [paired_1_files[0], paired_2_files[0]]:
                        logger.debug(f"  Adding file: {os.path.basename(file)} (Size: {os.path.getsize(file)} bytes)")
                        with open(file, 'r') as infile:
                            outfile.write(infile.read())
                
                # Verify the concatenated file
                if os.path.exists(concatenated_file) and os.path.getsize(concatenated_file) > 0:
                    logger.info(f"Successfully created concatenated file: {os.path.basename(concatenated_file)} "
                            f"(Size: {os.path.getsize(concatenated_file)} bytes)")
                    humann3_input_files.append(concatenated_file)
                else:
                    logger.error(f"Failed to create valid concatenated file for {sample}.")
            except Exception as e:
                logger.error(f"Error concatenating files for sample {sample}: {str(e)}")
        else:
            logger.warning(f"Sample {sample} doesn't have exactly one R1 and one R2 file after deduplication")

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

    logger.info(f"KneadData completed with {len(kneaddata_files)} total output files")

    # Prepare files for HUMAnN3
    logger.info("Preparing KneadData outputs for HUMAnN3...")
    sample_files = {}
    humann3_input_files = []

    # Organize files by sample
    for file in kneaddata_files:
        # Extract the base sample name - focus on our specific pattern
        filename = os.path.basename(file)
        
        # Check if file matches our target pattern
        if "_paired_1.fastq" in filename or "_paired_2.fastq" in filename:
            # Extract sample name by removing the suffix
            if "_paired_1.fastq" in filename:
                sample_name = filename.replace("_paired_1.fastq", "")
            else:
                sample_name = filename.replace("_paired_2.fastq", "")
                
            # Further cleanup if needed (e.g., remove any kneaddata prefix)
            if sample_name.endswith("_kneaddata"):
                sample_name = sample_name[:-10]  # Remove "_kneaddata" suffix
            
            # Initialize list for this sample if not exists
            if sample_name not in sample_files:
                sample_files[sample_name] = []
            
            # Add file to sample's list
            sample_files[sample_name].append(file)
            logger.debug(f"Added {filename} to sample {sample_name}")

    # Process each sample's paired files
    for sample, files in sample_files.items():
        logger.info(f"Processing sample {sample} with {len(files)} KneadData output files")
        
        # Skip if we don't have exactly 2 files (paired-end)
        if len(files) != 2:
            logger.warning(f"Sample {sample} has {len(files)} files instead of expected 2 paired files. Skipping.")
            continue
        
        # Sort files to ensure R1 comes before R2
        files.sort()  # This should put _paired_1 before _paired_2
        
        # Verify we have the correct paired files
        file1 = os.path.basename(files[0])
        file2 = os.path.basename(files[1])
        
        if not ("_paired_1.fastq" in file1 and "_paired_2.fastq" in file2):
            logger.warning(f"Files for sample {sample} don't match expected pattern: {file1}, {file2}. Skipping.")
            continue
        
        # Create concatenated file
        concatenated_file = os.path.join(os.path.dirname(files[0]), f"{sample}_paired_concat.fastq")
        logger.info(f"Concatenating paired files for sample {sample} to {os.path.basename(concatenated_file)}")
        
        try:
            # Create concatenated file
            with open(concatenated_file, 'w') as outfile:
                for file in files:
                    logger.debug(f"  Adding file: {os.path.basename(file)} (Size: {os.path.getsize(file)} bytes)")
                    with open(file, 'r') as infile:
                        outfile.write(infile.read())
            
            # Verify the concatenated file
            if os.path.exists(concatenated_file) and os.path.getsize(concatenated_file) > 0:
                logger.info(f"Successfully created concatenated file: {os.path.basename(concatenated_file)} "
                        f"(Size: {os.path.getsize(concatenated_file)} bytes)")
                humann3_input_files.append(concatenated_file)
            else:
                logger.error(f"Failed to create valid concatenated file for {sample}.")
        except Exception as e:
            logger.error(f"Error concatenating files for sample {sample}: {str(e)}")

    logger.info(f"Prepared {len(humann3_input_files)} input files for HUMAnN3")

    # Now check that we have valid files to run HUMAnN3 on
    if not humann3_input_files:
        logger.error("No valid input files prepared for HUMAnN3")
        return None

    # Continue with HUMAnN3 execution...
    # Replace the original kneaddata_files with our new concatenated files
    kneaddata_files = humann3_input_files

    # Step 2: Run HUMAnN3 in parallel
    logger.info("Starting HUMAnN3 step in parallel...")
    humann3_results = run_humann3_parallel(
        input_files=humann3_input_files,
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