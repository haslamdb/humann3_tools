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

    # Group files by sample name (remove _paired_1, _paired_2 suffixes)
    for file in kneaddata_files:
        # Extract the base sample name
        filename = os.path.basename(file)
        if "_paired_1.fastq" in filename or "_paired_2.fastq" in filename:
            base_name = filename.replace("_paired_1.fastq", "").replace("_paired_2.fastq", "")
            if base_name not in paired_files:
                paired_files[base_name] = []
            paired_files[base_name].append(file)

    # Concatenate paired files for each sample
    concatenated_files = []
    for sample, files in paired_files.items():
        if len(files) == 2:
            # Sort to ensure _paired_1 comes before _paired_2
            files.sort()
            concatenated_file = os.path.join(os.path.dirname(files[0]), f"{sample}_paired_concat.fastq")
            
            # Concatenate the files
            logger.info(f"Concatenating paired files for sample {sample}")
            with open(concatenated_file, 'w') as outfile:
                for file in files:
                    with open(file, 'r') as infile:
                        for line in infile:
                            outfile.write(line)
            
            concatenated_files.append(concatenated_file)
        else:
            # If not a proper pair, use the files individually
            logger.info(f"Found {len(files)} files for sample {sample}, using individually")
            concatenated_files.extend(files)

    # Use concatenated files for HUMAnN3
    humann3_input_files = concatenated_files if concatenated_files else kneaddata_files
    logger.info(f"Prepared {len(humann3_input_files)} input files for HUMAnN3")
    
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