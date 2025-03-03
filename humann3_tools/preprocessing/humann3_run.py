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
        cmd.extend(["--nucleotide-database", str(nucleotide_db)])
    if protein_db:
        cmd.extend(["--protein-database", str(protein_db)])
    
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
        'pathcoverage': None,
        'metaphlan': None 
    }

    # First check the main directory for standard HUMAnN3 outputs
    for file in os.listdir(output_dir):
        for output_type in output_files.keys():
            if output_type in file.lower() and file.endswith(".tsv"):
                output_files[output_type] = os.path.join(output_dir, file)
        
        # Recursively search for metaphlan output in all subdirectories
        for root, dirs, files in os.walk(output_dir):
            for file in files:
                if "metaphlan_bugs_list" in file.lower() and file.endswith(".tsv"):
                    output_files['metaphlan'] = os.path.join(root, file)
                    logger.info(f"Found MetaPhlAn output file: {os.path.join(root, file)}")
                    break  # Stop once we find the first matching file
                    
        # Log which files were found and which are missing
        found_files = [k for k, v in output_files.items() if v is not None]
        missing_files = [k for k, v in output_files.items() if v is None]
        
        logger.info(f"Found output files: {', '.join(found_files)}")
        if missing_files:
            logger.warning(f"Missing output files: {', '.join(missing_files)}")
        
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
    
    # Post-process to ensure we found metaphlan files
    for sample_id, sample_outputs in results.items():
        # If we didn't find the metaphlan file in the main function, try to find it now
        if isinstance(sample_outputs, dict) and sample_outputs.get('metaphlan') is None:
            sample_dir = os.path.join(output_dir, sample_id)
            if os.path.exists(sample_dir):
                for root, dirs, files in os.walk(sample_dir):
                    for file in files:
                        if "metaphlan_bugs_list" in file.lower() and file.endswith(".tsv"):
                            sample_outputs['metaphlan'] = os.path.join(root, file)
                            logger.info(f"Post-process: Found MetaPhlAn output for {sample_id}: {file}")
                            break
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
        
        # Add threads - ensure it's a string
        cmd.extend(["--threads", str(threads)])
        
        # Add database paths if provided
        if nucleotide_db:
            cmd.extend(["--nucleotide-database", str(nucleotide_db)])
        if protein_db:
            cmd.extend(["--protein-database", str(protein_db)])
        
        # Add any additional options
        if additional_options:
            for key, value in additional_options.items():
                if value is True:
                    cmd.append(f"--{key}")
                elif value is not None:
                    # Convert value to string to ensure it's not a complex object
                    cmd.extend([f"--{key}", str(value)])
        
        # Ensure all command elements are strings before joining
        str_cmd = [str(item) for item in cmd]
        
        # Run HUMAnN3
        logger.info(f"Running HUMAnN3 for sample {sample_name}: {' '.join(str_cmd)}")
        success = run_cmd(cmd, exit_on_error=False)
        
        if not success:
            logger.error(f"HUMAnN3 run failed for sample {sample_name}")
            continue
        
        # Find output files for this sample
        sample_outputs = {
            'genefamilies': None,
            'pathabundance': None,
            'pathcoverage': None,
            'metaphlan': None
        }
        
        # Check main directory first
        for file in os.listdir(sample_output_dir):
            for output_type in sample_outputs.keys():
                if output_type in file.lower() and file.endswith(".tsv"):
                    sample_outputs[output_type] = os.path.join(sample_output_dir, file)
        
        # Recursively search for metaphlan output in subdirectories
        for root, dirs, files in os.walk(sample_output_dir):
            for file in files:
                if "metaphlan_bugs_list" in file.lower() and file.endswith(".tsv"):
                    sample_outputs['metaphlan'] = os.path.join(root, file)
                    logger.info(f"Found MetaPhlAn output for {sample_name}: {file}")
                    break  # Stop after finding the first match
        
        # Log which files were found and which are missing
        found_files = [k for k, v in sample_outputs.items() if v is not None]
        missing_files = [k for k, v in sample_outputs.items() if v is None]
        
        logger.info(f"Sample {sample_name} output files found: {', '.join(found_files)}")
        if missing_files:
            logger.warning(f"Sample {sample_name} missing files: {', '.join(missing_files)}")
            
        output_files[sample_name] = sample_outputs
    
    logger.info(f"HUMAnN3 completed for {len(output_files)} samples")
    return output_files