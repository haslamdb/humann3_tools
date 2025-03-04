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

import os
import logging

from humann3_tools.preprocessing.kneaddata import (
    run_kneaddata,
    check_kneaddata_installation
)
from humann3_tools.preprocessing.humann3_run import (
    run_humann3,
    check_humann3_installation
)


def run_preprocessing_pipeline(
    input_files,
    output_dir,
    threads=1,
    kneaddata_dbs=None,
    nucleotide_db=None,
    protein_db=None,
    kneaddata_options=None,
    humann3_options=None,
    paired=False,
    kneaddata_output_dir=None,
    humann3_output_dir=None,
    logger=None
):
    """
    Run the full preprocessing pipeline: KneadData → HUMAnN3.

    This version simplifies:
      1) Avoiding duplicates,
      2) Finding correct paired files, and
      3) Logging only the essentials.

    Returns a dictionary with:
      'kneaddata_files': list of all final kneaddata files
      'humann3_results': dict of HUMAnN3 outputs per sample
    """
    if logger is None:
        logger = logging.getLogger("humann3_analysis")


    # 1. Check installations

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


    # 2. Create output directories

    if kneaddata_output_dir is None:
        kneaddata_output_dir = os.path.join(output_dir, "kneaddata_output")
    if humann3_output_dir is None:
        humann3_output_dir = os.path.join(output_dir, "humann3_output")

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(kneaddata_output_dir, exist_ok=True)
    os.makedirs(humann3_output_dir, exist_ok=True)

    logger.info(f"Using KneadData output directory: {kneaddata_output_dir}")
    logger.info(f"Using HUMAnN3 output directory: {humann3_output_dir}")


    # 3. Run KneadData

    logger.info(f"Starting KneadData step (paired={paired})...")

    kneaddata_files = []
    if paired:
        # Ensure an even number of input files for paired reads
        if len(input_files) % 2 != 0:
            logger.error("Paired mode requires an even number of input files")
            return None

        # Process input files in pairs
        for i in range(0, len(input_files), 2):
            pair = [input_files[i], input_files[i + 1]]
            sample_name = os.path.basename(pair[0]).split("_")[0]
            sample_outdir = os.path.join(kneaddata_output_dir, sample_name)
            os.makedirs(sample_outdir, exist_ok=True)

            logger.info(f"KneadData on paired files for sample {sample_name}")
            results = run_kneaddata(
                input_files=pair,
                output_dir=sample_outdir,
                threads=threads,
                reference_dbs=kneaddata_dbs,
                paired=True,
                additional_options=kneaddata_options,
                logger=logger,
            )
            if results:
                kneaddata_files.extend(results)
            else:
                logger.warning(f"No KneadData output for sample {sample_name}")

    else:
        # Single-end mode
        for f in input_files:
            sample_name = os.path.basename(f).split("_")[0]
            sample_outdir = os.path.join(kneaddata_output_dir, sample_name)
            os.makedirs(sample_outdir, exist_ok=True)

            logger.info(f"KneadData on single-end file for sample {sample_name}")
            results = run_kneaddata(
                input_files=[f],
                output_dir=sample_outdir,
                threads=threads,
                reference_dbs=kneaddata_dbs,
                paired=False,
                additional_options=kneaddata_options,
                logger=logger,
            )
            if results:
                kneaddata_files.extend(results)
            else:
                logger.warning(f"No KneadData output for sample {sample_name}")

    if not kneaddata_files:
        logger.error("KneadData produced no output; cannot continue.")
        return None

    logger.info(f"KneadData completed with {len(kneaddata_files)} total output files")


    # 4. Organize KneadData outputs by sample

    logger.info("Preparing KneadData outputs for HUMAnN3...")
    sample_files = {}

    for f in kneaddata_files:
        fname = os.path.basename(f)
        # Only process .fastq files
        if not fname.endswith(".fastq"):
            continue

        # Identify if file is R1 or R2 based on naming
        if "_paired_1.fastq" in fname:
            sname = fname.replace("_paired_1.fastq", "")
        elif "_paired_2.fastq" in fname:
            sname = fname.replace("_paired_2.fastq", "")
        else:
            logger.debug(f"Skipping non-paired file: {fname}")
            continue

        # Remove trailing "_kneaddata" if it exists
        if sname.endswith("_kneaddata"):
            sname = sname[: -len("_kneaddata")]

        # Use a set to avoid duplicates
        if sname not in sample_files:
            sample_files[sname] = set()
        sample_files[sname].add(f)


    # 5. For each sample, find exactly one _paired_1.fastq and one _paired_2.fastq,
    #    then concatenate them for HUMAnN3.

    humann3_input_files = []

    for sample_name, file_set in sample_files.items():
        flist = sorted(list(file_set))
        logger.info(f"Sample {sample_name} has {len(flist)} candidate files")

        r1_list = [x for x in flist if "_paired_1.fastq" in os.path.basename(x)]
        r2_list = [x for x in flist if "_paired_2.fastq" in os.path.basename(x)]

        if len(r1_list) == 1 and len(r2_list) == 1:
            r1, r2 = r1_list[0], r2_list[0]
            logger.info(
                f"Concatenating R1 and R2 for sample {sample_name} -> single HUMAnN3 input"
            )
            out_dir = os.path.dirname(r1)
            concat_file = os.path.join(out_dir, f"{sample_name}_paired_concat.fastq")

            try:
                with open(concat_file, "w") as outfile:
                    for rf in (r1, r2):
                        size = os.path.getsize(rf)
                        logger.debug(f"  Adding {os.path.basename(rf)} ({size} bytes)")
                        with open(rf, "r") as infile:
                            outfile.write(infile.read())

                if os.path.getsize(concat_file) > 0:
                    logger.info(f"Created {os.path.basename(concat_file)} for {sample_name}")
                    humann3_input_files.append(concat_file)
                else:
                    logger.error(f"Concat file was empty for {sample_name}")
            except Exception as e:
                logger.error(f"Concat error for {sample_name}: {str(e)}")
        else:
            logger.warning(
                f"Sample {sample_name}: found {len(r1_list)} R1 files and {len(r2_list)} R2 files. Skipping."
            )

    logger.info(f"Prepared {len(humann3_input_files)} input files for HUMAnN3")

    if not humann3_input_files:
        logger.error("No valid input files after KneadData; cannot run HUMAnN3.")
        return None


    # 6. Run HUMAnN3
    logger.info("Starting HUMAnN3...")
    humann3_results = run_humann3(
        input_files=humann3_input_files,
        output_dir=humann3_output_dir,
        threads=threads,
        nucleotide_db=nucleotide_db,
        protein_db=protein_db,
        additional_options=humann3_options,
        logger=logger
    )
    if not humann3_results:
        logger.error("HUMAnN3 step failed.")
        return None

    logger.info(f"HUMAnN3 completed for {len(humann3_results)} samples.")

    # Return combined results
    return {
        "kneaddata_files": kneaddata_files,
        "humann3_results": humann3_results,
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