
import os
import subprocess
import logging
from humann3_tools.utils.cmd_utils import run_cmd
from humann3_tools.logger import log_print
from humann3_tools.utils.resource_utils import track_peak_memory 

def check_kneaddata_installation():
    """Check if KneadData is installed and available."""
    try:
        result = subprocess.run(["kneaddata", "--version"], 
                               capture_output=True, text=True, check=False)
        if result.returncode == 0:
            return True, result.stdout.strip()
        return False, "KneadData command exists but returned an error"
    except FileNotFoundError:
        return False, "KneadData not found in PATH"
    

def process_single_sample_kneaddata(input_file, sample_id=None, output_dir=None, 
                                   threads=1, reference_dbs=None, paired_file=None, 
                                   additional_options=None, logger=None):
    """Process a single sample with KneadData.
    
    Args:
        input_file: Input FASTQ file path
        sample_id: Sample identifier
        output_dir: Output directory path
        threads: Number of threads to use
        reference_dbs: List of reference database paths
        paired_file: Paired FASTQ file path (for paired-end data)
        additional_options: Dictionary of additional options
        logger: Logger instance
        
    Returns:
        List of output FASTQ files
    """
    if logger is None:
        logger = logging.getLogger('humann3_analysis')
    
    if sample_id is None:
        sample_id = os.path.basename(input_file).split('.')[0]
    
    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), "kneaddata_output", sample_id)
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Build command
    cmd = ["kneaddata"]
    
    # Add input file
    cmd.extend(["-i1", input_file])
    
    # Add paired file if provided
    if paired_file:
        cmd.extend(["-i2", paired_file])
    
    # Add output directory
    cmd.extend(["-o", output_dir])
    
    # Add threads (per sample)
    cmd.extend(["-t", str(threads)])
    
    # Add reference database(s) if provided
    if reference_dbs:
        # Handle both string and list inputs
        if isinstance(reference_dbs, str):
            cmd.extend(["-db", reference_dbs])
        else:
            # Add each reference database with its own -db flag
            for db in reference_dbs:
                cmd.extend(["-db", db])
    
    # Add additional options
    if additional_options:
        for key, value in additional_options.items():
            if value is True:
                cmd.append(f"--{key}")
            elif value is not None and value != "":
                cmd.extend([f"--{key}", str(value)])
    
    # Run KneadData
    logger.info(f"Running KneadData for sample {sample_id}: {' '.join(str(x) for x in cmd)}")
    success = run_cmd(cmd, exit_on_error=False)
    
    if not success:
        logger.error(f"KneadData run failed for sample {sample_id}")
        return None
    
    # Find output files - look for all non-contaminant FASTQ files
    output_files = []
    for file in os.listdir(output_dir):
        file_path = os.path.join(output_dir, file)
        # Check if it's a non-empty file
        if os.path.isfile(file_path) and os.path.getsize(file_path) > 0:
            # Check if it's a fastq file and not a contaminant file
            if file.endswith(".fastq") and "contam" not in file:
                output_files.append(file_path)
                logger.debug(f"Found KneadData output file: {file} (size: {os.path.getsize(file_path)} bytes)")
    
    logger.info(f"KneadData completed for sample {sample_id} with {len(output_files)} output files")
    return output_files

# Parallel processing
@track_peak_memory
def run_kneaddata_parallel(input_files, output_dir, threads=1, max_parallel=None, 
                          reference_dbs=None, paired=False, additional_options=None, 
                          logger=None):
    """
    Run KneadData on multiple samples in parallel.
    
    Args:
        input_files: List of input FASTQ files (single-end) or list of pairs for paired-end
        output_dir: Base directory for outputs
        threads: Number of threads per sample
        max_parallel: Maximum number of parallel samples (None = CPU count)
        reference_dbs: List of reference database paths
        paired: Whether input is paired-end
        additional_options: Dict of additional KneadData options
        logger: Logger instance
        
    Returns:
        Dict mapping sample IDs to output files
    """
    from humann3_tools.preprocessing.parallel import run_parallel
    
    if logger is None:
        logger = logging.getLogger('humann3_analysis')
    
    # Create sample list based on input type
    sample_list = []
    if paired:
        # Ensure we have pairs of files
        if len(input_files) % 2 != 0:
            logger.error("Paired mode requires an even number of input files")
            return {}
        
        # Pair files in order (assuming they're in order: R1, R2, R1, R2, etc.)
        for i in range(0, len(input_files), 2):
            r1_file = input_files[i]
            r2_file = input_files[i+1]
            # Try to extract sample name from filename
            sample_name = None
            r1_basename = os.path.basename(r1_file)
            
            # Try common naming patterns
            if "_R1" in r1_basename:
                sample_name = r1_basename.split("_R1")[0]
            elif "_1." in r1_basename:
                sample_name = r1_basename.split("_1.")[0]
            else:
                # Fallback to using the filename without extension
                sample_name = os.path.splitext(r1_basename)[0]
            
            # Store the paired files
            sample_list.append((sample_name, {"input_file": r1_file, "paired_file": r2_file}))
            logger.info(f"Paired sample {sample_name}: {r1_file} + {r2_file}")
    else:
        # Single-end reads
        for file in input_files:
            # Try to extract sample name from filename
            basename = os.path.basename(file)
            sample_name = os.path.splitext(basename)[0]
            sample_list.append((sample_name, {"input_file": file}))
            logger.info(f"Single-end sample {sample_name}: {file}")
    
    # Prepare common arguments for all samples
    kwargs = {
        'output_dir': output_dir,
        'threads': threads,
        'reference_dbs': reference_dbs,
        'additional_options': additional_options,
        'logger': logger
    }
    
    # Define a wrapper function to handle the different input format
    def process_sample_wrapper(sample_data, sample_id, **kwargs):
        input_file = sample_data["input_file"]
        paired_file = sample_data.get("paired_file")
        return process_single_sample_kneaddata(
            input_file=input_file, 
            paired_file=paired_file, 
            sample_id=sample_id, 
            **kwargs
        )
    
    # Run in parallel
    results = run_parallel(sample_list, process_sample_wrapper, 
                          max_workers=max_parallel, **kwargs)
    
    return results

def run_kneaddata(input_files, output_dir, threads=1, reference_dbs=None, 
                 paired=False, additional_options=None, logger=None):
    """
    Run KneadData on input sequence files.
    
    Args:
        input_files: List of input FASTQ files
        output_dir: Directory for KneadData output
        threads: Number of threads to use
        reference_dbs: List of reference database paths
        paired: Whether input files are paired
        additional_options: Dict of additional KneadData options
        logger: Logger instance
        
    Returns:
        List of output FASTQ files
    """
    if logger is None:
        logger = logging.getLogger('humann3_analysis')
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Basic command
    cmd = ["kneaddata"]
    
    # Handle paired vs single end using correct flags (-i1/-i2)
    if paired and len(input_files) >= 2:
        # For paired data
        cmd.extend(["-i1", input_files[0], "-i2", input_files[1]])
    elif len(input_files) >= 1:
        # For single-end data
        cmd.extend(["-i1", input_files[0]])
    else:
        logger.error("No input files provided for KneadData")
        return []
    
    # Add output directory
    cmd.extend(["-o", output_dir])
    
    # Add threads
    cmd.extend(["-t", str(threads)])
    
    # Add reference database(s) 
    if reference_dbs:
        # Handle both string and list inputs
        if isinstance(reference_dbs, str):
            cmd.extend(["-db", reference_dbs])
        else:
            # Add each reference database with its own -db flag
            for db in reference_dbs:
                cmd.extend(["-db", db])
    
    # Add any additional options
    if additional_options:
        for key, value in additional_options.items():
            if value is True:
                cmd.append(f"--{key}")
            elif value is not None:
                cmd.extend([f"--{key}", str(value)])
    
    # Run KneadData
    logger.info(f"Running KneadData: {' '.join(str(x) for x in cmd)}")
    success = run_cmd(cmd, exit_on_error=False)
    
    if not success:
        logger.error("KneadData run failed")
        return []
    
    # Find output files - look for all non-contaminant FASTQ files
    output_files = []
    for file in os.listdir(output_dir):
        file_path = os.path.join(output_dir, file)
        # Check if it's a non-empty file
        if os.path.isfile(file_path) and os.path.getsize(file_path) > 0:
            # Check if it's a fastq file and not a contaminant file
            if file.endswith(".fastq") and "contam" not in file:
                output_files.append(file_path)
                logger.debug(f"Found KneadData output file: {file} (size: {os.path.getsize(file_path)} bytes)")
    
    if not output_files:
        logger.warning("No valid output files found in KneadData output directory")
        # List all files in the directory for debugging
        for file in os.listdir(output_dir):
            file_path = os.path.join(output_dir, file)
            if os.path.isfile(file_path):
                logger.debug(f"File in output dir: {file} (size: {os.path.getsize(file_path)} bytes)")
    
    logger.info(f"KneadData completed with {len(output_files)} output files")
    return output_files