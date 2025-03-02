# humann3_tools/utils/resource_utils.py
import os
import multiprocessing
import logging
import psutil

def calculate_optimal_resources(available_threads, num_samples, min_threads_per_sample=1):
    """
    Calculate optimal thread allocation between samples and processes.
    
    Args:
        available_threads: Total available CPU threads
        num_samples: Number of samples to process
        min_threads_per_sample: Minimum threads to allocate per sample
        
    Returns:
        Tuple of (threads_per_sample, max_parallel_samples)
    """
    if available_threads is None:
        available_threads = multiprocessing.cpu_count()
    
    # Start with maximum parallelism
    max_parallel = min(num_samples, available_threads)
    threads_per_sample = max(min_threads_per_sample, available_threads // max_parallel)
    
    # Recalculate max_parallel based on threads_per_sample
    max_parallel = min(max_parallel, available_threads // threads_per_sample)
    
    return threads_per_sample, max_parallel

def get_memory_usage():
    """Get current memory usage in MB."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / (1024 * 1024)  # Convert to MB

def log_resource_usage(logger, sample_id=None):
    """Log current resource usage."""
    mem_usage = get_memory_usage()
    cpu_percent = psutil.cpu_percent()
    
    message = f"Resource usage: Memory: {mem_usage:.2f} MB, CPU: {cpu_percent}%"
    if sample_id:
        message = f"Sample {sample_id}: {message}"
        
    logger.info(message)

def estimate_memory_requirements(num_samples, tool="humann3"):
    """
    Estimate memory requirements for a tool based on number of samples.
    
    Args:
        num_samples: Number of samples to process
        tool: Tool name ("kneaddata" or "humann3")
        
    Returns:
        Estimated memory in MB
    """
    # These are rough estimates and should be adjusted based on real-world usage
    if tool.lower() == "kneaddata":
        # KneadData is less memory-intensive
        return 2000 + (500 * num_samples)  # Base 2GB + 500MB per sample
    elif tool.lower() == "humann3":
        # HUMAnN3 can be very memory-intensive
        return 6000 + (2000 * num_samples)  # Base 6GB + 2GB per sample
    else:
        return 4000 + (1000 * num_samples)  # Generic estimate

def check_resource_availability(required_memory, required_threads):
    """
    Check if the system has enough resources.
    
    Args:
        required_memory: Required memory in MB
        required_threads: Required CPU threads
        
    Returns:
        Tuple of (memory_ok, cpu_ok, available_memory, available_threads)
    """
    available_memory = psutil.virtual_memory().available / (1024 * 1024)  # MB
    available_threads = multiprocessing.cpu_count()
    
    memory_ok = available_memory >= required_memory
    cpu_ok = available_threads >= required_threads
    
    return memory_ok, cpu_ok, available_memory, available_threads