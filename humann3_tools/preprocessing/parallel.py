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
    sample_id, file_path = sample_tuple
    logger = logging.getLogger('humann3_analysis')
    
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

