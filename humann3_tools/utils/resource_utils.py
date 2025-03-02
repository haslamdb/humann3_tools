# Add these to your existing resource_utils.py file

def monitor_memory_usage(logger, threshold_mb=1000, interval=60):
    """
    Start a background thread to monitor memory usage and log warnings 
    if it exceeds the threshold.
    
    Args:
        logger: Logger instance
        threshold_mb: Memory threshold in MB to trigger warnings
        interval: Check interval in seconds
        
    Returns:
        Thread object (can be used to stop monitoring)
    """
    import threading
    import time
    
    stop_event = threading.Event()
    
    def monitoring_worker():
        while not stop_event.is_set():
            mem_usage = get_memory_usage()
            if mem_usage > threshold_mb:
                logger.warning(f"High memory usage detected: {mem_usage:.2f} MB (threshold: {threshold_mb} MB)")
            time.sleep(interval)
    
    monitor_thread = threading.Thread(target=monitoring_worker, daemon=True)
    monitor_thread.start()
    
    return monitor_thread, stop_event

def stop_memory_monitoring(thread_data):
    """Stop memory monitoring thread."""
    if thread_data:
        thread, stop_event = thread_data
        stop_event.set()
        thread.join(timeout=1)

def track_peak_memory(func):
    """
    Decorator to track peak memory usage during function execution.
    
    Usage:
        @track_peak_memory
        def my_function(logger, ...):
            # function code
    """
    from functools import wraps
    
    @wraps(func)
    def wrapper(*args, **kwargs):
        # Find logger in args or kwargs
        logger = None
        for arg in args:
            if isinstance(arg, logging.Logger):
                logger = arg
                break
                
        if not logger and 'logger' in kwargs:
            logger = kwargs['logger']
        
        if not logger:
            logger = logging.getLogger('humann3_analysis')
        
        # Record starting memory
        start_mem = get_memory_usage()
        logger.info(f"Starting memory: {start_mem:.2f} MB")
        
        # Track peak memory
        peak_mem = start_mem
        
        def memory_checker():
            nonlocal peak_mem
            current = get_memory_usage()
            if current > peak_mem:
                peak_mem = current
                
        # Start periodic checking
        import threading
        import time
        
        stop_event = threading.Event()
        
        def check_loop():
            while not stop_event.is_set():
                memory_checker()
                time.sleep(1)
                
        monitor_thread = threading.Thread(target=check_loop, daemon=True)
        monitor_thread.start()
        
        try:
            # Run the original function
            result = func(*args, **kwargs)
            return result
        finally:
            # Stop monitoring
            stop_event.set()
            monitor_thread.join(timeout=1)
            
            # Log results
            end_mem = get_memory_usage()
            logger.info(f"Memory usage: start={start_mem:.2f} MB, peak={peak_mem:.2f} MB, end={end_mem:.2f} MB")
            
    return wrapper

def limit_memory_usage(max_memory_mb=None):
    """
    Try to limit memory usage of the current process.
    Note: This only works on Unix/Linux systems.
    
    Args:
        max_memory_mb: Maximum memory in MB (None = no limit)
    
    Returns:
        True if limit was set, False otherwise
    """
    if max_memory_mb is None:
        return False
        
    try:
        import resource
        
        # Convert MB to bytes
        max_memory_bytes = max_memory_mb * 1024 * 1024
        
        # Set soft limit
        resource.setrlimit(resource.RLIMIT_AS, (max_memory_bytes, max_memory_bytes))
        return True
    except (ImportError, ValueError, resource.error):
        return False