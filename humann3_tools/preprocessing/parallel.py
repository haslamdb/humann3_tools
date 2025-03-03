# humann3_tools/preprocessing/parallel.py
import os
import time
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

def process_sample_parallel(sample_tuple, function, **kwargs):
    """
    Process a single sample with the provided function.
    
    Args:
        sample_tuple: Tuple of (sample_id, file_path)
        function: Function to run on the sample
        **kwargs: Additional arguments to pass to the function
        
    Returns:
        Tuple of (sample_id, result)
    """
    # Extract sample_id and file_path correctly, handling different formats
    # Print for debugging
    logger = logging.getLogger('humann3_analysis')
    logger.debug(f"Sample tuple received: {sample_tuple}")
    
    # Handle different tuple formats
    if isinstance(sample_tuple, tuple) and len(sample_tuple) == 2:
        sample_id, file_path = sample_tuple
        
        # If file_path is itself a tuple (for paired files), handle it
        if isinstance(file_path, tuple):
            # For kneaddata, the second argument should be paired_file
            # Modify the function call to extract the paired files
            r1_file, r2_file = file_path
            logger.info(f"Processing paired files for sample {sample_id}: {os.path.basename(r1_file)}, {os.path.basename(r2_file)}")
            
            # For paired files, r1_file is the main input, r2_file is the paired input
            start_time = time.time()
            result = function(r1_file, sample_id=sample_id, paired_file=r2_file, **kwargs)
            
            elapsed = time.time() - start_time
            logger.info(f"Finished processing sample {sample_id} in {elapsed:.2f} seconds")
            
            return sample_id, result
    
    # Original code for non-paired files
    sample_id, file_path = sample_tuple
    
    start_time = time.time()
    logger.info(f"Started processing sample {sample_id}")
    
    result = function(file_path, sample_id=sample_id, **kwargs)
    
    elapsed = time.time() - start_time
    logger.info(f"Finished processing sample {sample_id} in {elapsed:.2f} seconds")
    
    return sample_id, result

def run_parallel(sample_list, function, max_workers=None, **kwargs):
    """
    Run a function on multiple samples in parallel.
    
    Args:
        sample_list: List of tuples (sample_id, file_path)
        function: Function to run on each sample
        max_workers: Maximum number of parallel processes (None = CPU count)
        **kwargs: Additional arguments to pass to the function
        
    Returns:
        Dictionary mapping sample_ids to results
    """
    logger = logging.getLogger('humann3_analysis')
    logger.info(f"Starting parallel processing of {len(sample_list)} samples with {max_workers} workers")
    
    results = {}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit tasks
        future_to_sample = {
            executor.submit(process_sample_parallel, sample_tuple, function, **kwargs): sample_tuple[0]
            for sample_tuple in sample_list
        }
        
        # Process results as they complete
        for future in as_completed(future_to_sample):
            sample_id = future_to_sample[future]
            try:
                sample_id, result = future.result()
                results[sample_id] = result
                logger.info(f"Successfully processed sample {sample_id}")
            except Exception as e:
                logger.error(f"Error processing sample {sample_id}: {str(e)}")
    
    logger.info(f"Completed parallel processing. Successfully processed {len(results)} of {len(sample_list)} samples")
    return results

def run_parallel(sample_list, function, max_workers=None, **kwargs):
    """Run a function on multiple samples in parallel with progress bar."""
    logger = logging.getLogger('humann3_analysis')
    logger.info(f"Starting parallel processing of {len(sample_list)} samples with {max_workers} workers")
    
    results = {}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit tasks
        future_to_sample = {
            executor.submit(process_sample_parallel, sample_tuple, function, **kwargs): sample_tuple[0]
            for sample_tuple in sample_list
        }
        
        # Process results with progress bar
        for future in tqdm(as_completed(future_to_sample), total=len(future_to_sample), 
                          desc="Processing samples", unit="sample"):
            sample_id = future_to_sample[future]
            try:
                sample_id, result = future.result()
                results[sample_id] = result
            except Exception as e:
                logger.error(f"Error processing sample {sample_id}: {str(e)}")
    
    return results

