# humann3_tools/preprocessing/pipeline.py
import os
import logging
from humann3_tools.preprocessing.kneaddata import run_kneaddata, check_kneaddata_installation
from humann3_tools.preprocessing.humann3_run import run_humann3, check_humann3_installation
from humann3_tools.logger import log_print
from humann3_tools.utils.resource_utils import (
    track_peak_memory, 
    monitor_memory_usage, 
    stop_memory_monitoring
)

def run_preprocessing_pipeline(input_files, output_dir, threads=1, 
                              kneaddata_db=None, nucleotide_db=None, protein_db=None,
                              kneaddata_options=None, humann3_options=None, 
                              paired=False, logger=None):
    """
    Run the full preprocessing pipeline: KneadData → HUMAnN3.
    
    Args:
        input_files: List of input FASTQ files
        output_dir: Base directory for outputs
        threads: Number of threads to use
        kneaddata_db: Path to KneadData reference database
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
    
    # Create output directories
    kneaddata_output = os.path.join(output_dir, "kneaddata_output")
    humann3_output = os.path.join(output_dir, "humann3_output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Step 1: Run KneadData
    logger.info("Starting KneadData step...")
    kneaddata_files = run_kneaddata(
        input_files=input_files,
        output_dir=kneaddata_output,
        threads=threads,
        reference_db=kneaddata_db,
        paired=paired,
        additional_options=kneaddata_options,
        logger=logger
    )
    
    if not kneaddata_files:
        logger.error("KneadData step failed, stopping pipeline")
        return None
    
    logger.info(f"KneadData completed successfully with {len(kneaddata_files)} output files")
    
    # Step 2: Run HUMAnN3
    logger.info("Starting HUMAnN3 step...")
    humann3_results = run_humann3(
        input_files=kneaddata_files,
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


# Add support for parallels processing

def run_preprocessing_pipeline_parallel(input_files, output_dir, threads_per_sample=1, 
                                       max_parallel=None, kneaddata_db=None, 
                                       nucleotide_db=None, protein_db=None,
                                       kneaddata_options=None, humann3_options=None, 
                                       paired=False, logger=None):
    """
    Run the full preprocessing pipeline in parallel: KneadData → HUMAnN3.
    
    Args:
        input_files: List of input FASTQ files
        output_dir: Base directory for outputs
        threads_per_sample: Number of threads per sample
        max_parallel: Maximum number of parallel samples (None = CPU count)
        kneaddata_db: Path to KneadData reference database
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
    
    logger.info(f"Starting parallel preprocessing pipeline with {len(input_files)} input files")
    logger.info(f"KneadData version: {kneaddata_version}")
    logger.info(f"HUMAnN3 version: {humann3_version}")
    
    # Create output directories
    kneaddata_output = os.path.join(output_dir, "kneaddata_output")
    humann3_output = os.path.join(output_dir, "humann3_output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Step 1: Run KneadData in parallel
    logger.info("Starting KneadData step in parallel...")
    kneaddata_results = run_kneaddata_parallel(
        input_files=input_files,
        output_dir=kneaddata_output,
        threads=threads_per_sample,
        max_parallel=max_parallel,
        reference_db=kneaddata_db,
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

@track_peak_memory
def run_preprocessing_pipeline_parallel(input_files, output_dir, threads_per_sample=1, 
                                      max_parallel=None, kneaddata_db=None, 
                                      nucleotide_db=None, protein_db=None,
                                      kneaddata_options=None, humann3_options=None, 
                                      paired=False, logger=None):
    """Run the full preprocessing pipeline in parallel."""
    if logger is None:
        logger = logging.getLogger('humann3_analysis')
    
    # Start memory monitoring
    monitor = monitor_memory_usage(logger, threshold_mb=8000, interval=30)
    
    try:
        # Your existing code here...
        
        # Step 1: Run KneadData in parallel
        logger.info("Starting KneadData step in parallel...")
        kneaddata_results = run_kneaddata_parallel(...)
        
        # Step 2: Run HUMAnN3 in parallel
        logger.info("Starting HUMAnN3 step in parallel...")
        humann3_results = run_humann3_parallel(...)
        
        # Return results
        return {...}
    
    finally:
        # Stop memory monitoring
        stop_memory_monitoring(monitor)