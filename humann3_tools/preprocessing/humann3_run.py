# humann3_tools/preprocessing/humann3_run.py
import os
import subprocess
import logging
from humann3_tools.utils.cmd_utils import run_cmd
from humann3_tools.logger import log_print
from humann3_tools.utils.resource_utils import track_peak_memory

def check_humann3_installation():
    """Check if HUMAnN3 is installed and available."""
    try:
        result = subprocess.run(["humann", "--version"], 
                               capture_output=True, text=True, check=False)
        if result.returncode == 0:
            return True, result.stdout.strip()
        return False, "HUMAnN3 command exists but returned an error"
    except FileNotFoundError:
        return False, "HUMAnN3 not found in PATH"
    

def process_single_sample_humann3(input_file, sample_id=None, output_dir=None, 
                                 threads=1, nucleotide_db=None, protein_db=None, 
                                 additional_options=None, logger=None):
    """Process a single sample with HUMAnN3."""
    if logger is None:
        logger = logging.getLogger('humann3_analysis')
    
    if sample_id is None:
        sample_id = os.path.basename(input_file).split('.')[0]
    
    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), "humann3_output", sample_id)
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Build command
    cmd = ["humann", "--input", input_file, "--output", output_dir]
    
    # Add threads (per sample)
    cmd.extend(["--threads", str(threads)])
    
    # Add database paths if provided
    if nucleotide_db:
        cmd.extend(["--nucleotide-database", nucleotide_db])
    if protein_db:
        cmd.extend(["--protein-database", protein_db])
    
    # Add additional options
    if additional_options:
        for key, value in additional_options.items():
            if value is True:
                cmd.append(f"--{key}")
            elif value is not None and value != "":
                cmd.extend([f"--{key}", str(value)])
    
    # Run HUMAnN3
    logger.info(f"Running HUMAnN3 for sample {sample_id}")
    success = run_cmd(cmd, exit_on_error=False)
    
    if not success:
        logger.error(f"HUMAnN3 run failed for sample {sample_id}")
        return None
    
    # Find output files
    output_files = {
        'genefamilies': None,
        'pathabundance': None,
        'pathcoverage': None
    }
    
    for file in os.listdir(output_dir):
        for output_type in output_files.keys():
            if output_type in file.lower() and file.endswith(".tsv"):
                output_files[output_type] = os.path.join(output_dir, file)
    
    logger.info(f"HUMAnN3 completed for sample {sample_id}")
    return output_files

# Add support for parallel processing
@track_peak_memory
def run_humann3_parallel(input_files, output_dir, threads=1, max_parallel=None,
                        nucleotide_db=None, protein_db=None, additional_options=None, 
                        logger=None):
    """
    Run HUMAnN3 on multiple samples in parallel.
    
    Args:
        input_files: List of input FASTQ files from KneadData
        output_dir: Base directory for outputs
        threads: Number of threads per sample
        max_parallel: Maximum number of parallel samples (None = CPU count)
        nucleotide_db: Path to nucleotide database
        protein_db: Path to protein database
        additional_options: Dict of additional HUMAnN3 options
        logger: Logger instance
        
    Returns:
        Dict mapping sample IDs to output files
    """
    from humann3_tools.preprocessing.parallel import run_parallel
    
    if logger is None:
        logger = logging.getLogger('humann3_analysis')
    
    # Create sample list
    sample_list = []
    for file in input_files:
        sample_name = os.path.basename(file).split('_')[0]
        sample_list.append((sample_name, file))
    
    # Prepare common arguments for all samples
    kwargs = {
        'output_dir': output_dir,
        'threads': threads,
        'nucleotide_db': nucleotide_db,
        'protein_db': protein_db,
        'additional_options': additional_options,
        'logger': logger
    }
    
    # Run in parallel
    results = run_parallel(sample_list, process_single_sample_humann3, 
                          max_workers=max_parallel, **kwargs)
    
    return results


def run_humann3(input_files, output_dir, threads=1, nucleotide_db=None, 
               protein_db=None, additional_options=None, logger=None):
    """
    Run HUMAnN3 on input sequence files.
    
    Args:
        input_files: List of input FASTQ files from KneadData
        output_dir: Directory for HUMAnN3 output
        threads: Number of threads to use
        nucleotide_db: Path to nucleotide database
        protein_db: Path to protein database
        additional_options: Dict of additional HUMAnN3 options
        logger: Logger instance
        
    Returns:
        Dict of output file paths by type
    """
    if logger is None:
        logger = logging.getLogger('humann3_analysis')
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Process each input file
    output_files = {}
    for input_file in input_files:
        # Extract sample name from file path
        sample_name = os.path.basename(input_file).split("_")[0]
        sample_output_dir = os.path.join(output_dir, sample_name)
        os.makedirs(sample_output_dir, exist_ok=True)
        
        # Basic command
        cmd = ["humann", "--input", input_file, "--output", sample_output_dir]
        
        # Add threads
        cmd.extend(["--threads", str(threads)])
        
        # Add database paths if provided
        if nucleotide_db:
            cmd.extend(["--nucleotide-database", nucleotide_db])
        if protein_db:
            cmd.extend(["--protein-database", protein_db])
        
        # Add any additional options
        if additional_options:
            for key, value in additional_options.items():
                if value is True:
                    cmd.append(f"--{key}")
                elif value is not None:
                    cmd.extend([f"--{key}", str(value)])
        
        # Run HUMAnN3
        logger.info(f"Running HUMAnN3 for sample {sample_name}: {' '.join(cmd)}")
        success = run_cmd(cmd, exit_on_error=False)
        
        if not success:
            logger.error(f"HUMAnN3 run failed for sample {sample_name}")
            continue
        
        # Find output files for this sample
        sample_outputs = {
            'genefamilies': None,
            'pathabundance': None,
            'pathcoverage': None
        }
        
        for file in os.listdir(sample_output_dir):
            for output_type in sample_outputs.keys():
                if output_type in file and file.endswith(".tsv"):
                    sample_outputs[output_type] = os.path.join(sample_output_dir, file)
        
        output_files[sample_name] = sample_outputs
    
    logger.info(f"HUMAnN3 completed for {len(output_files)} samples")
    return output_files