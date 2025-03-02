# humann3_tools/humann3_tools/humann3/gene_processing.py
import os
import logging

from humann3_tools.logger import log_print
from humann3_tools.utils.cmd_utils import run_cmd

def process_gene_families(valid_samples, gene_dir, output_dir, output_prefix, selected_columns=None):
    """
    Process gene families: copy, normalize, join, and split.
    Returns path to unstratified gene family file.
    
    Args:
        valid_samples: List of tuples (sample, filepath) for valid samples
        gene_dir: Directory containing original gene family files
        output_dir: Directory where processed files will be stored
        output_prefix: Prefix for output filenames
        selected_columns: Optional dict with sample key column selections
        
    Returns:
        Path to unstratified gene family file, or None if processing failed
    """
    log_print("PROCESSING GENE FAMILY FILES", level='info')
    gene_families_out = os.path.join(output_dir, "genes", output_prefix)
    gene_families_norm = os.path.join(gene_families_out, "Normalized")
    os.makedirs(gene_families_out, exist_ok=True)
    os.makedirs(gene_families_norm, exist_ok=True)
    
    if not valid_samples:
        log_print("No valid gene family files to process", level='warning')
        return None
    
    if selected_columns:
        column_info_file = os.path.join(gene_families_out, "column_selections.txt")
        try:
            with open(column_info_file, "w") as f:
                f.write(f"sample_id_column: {selected_columns['sample_id']}\n")
                f.write("grouping_columns:\n")
                for col in selected_columns.get('grouping', {}):
                    f.write(f"  - {col}\n")
        except Exception as e:
            log_print(f"Warning: Could not save column selection info: {e}", level='warning')
    
    processed_count = 0
    for (sample, src_path) in valid_samples:
        dst = os.path.join(gene_families_out, f"{sample}_genefamilies.tsv")
        if not run_cmd(["cp", src_path, dst], exit_on_error=False):
            continue
        
        out_cpm = os.path.join(gene_families_out, f"{sample}_genefamilies-cpm.tsv")
        if run_cmd([
            "humann_renorm_table",
            "--input", dst,
            "--output", out_cpm,
            "--units", "cpm",
            "--update-snames"
        ], exit_on_error=False):
            if run_cmd(["mv", out_cpm, gene_families_norm], exit_on_error=False):
                processed_count += 1
    
    if processed_count == 0:
        log_print("WARNING: No gene family files processed successfully", level='warning')
        return None
    
    norm_files = [f for f in os.listdir(gene_families_norm) if f.endswith("-cpm.tsv")]
    if not norm_files:
        log_print("WARNING: No normalized gene family files to join", level='warning')
        return None
    
    joined_output = os.path.join(gene_families_out, f"{output_prefix}_genefamilies-cpm.tsv")
    run_cmd([
        "humann_join_tables",
        "-i", gene_families_norm,
        "-o", joined_output
    ])
    
    if not os.path.exists(joined_output):
        log_print("WARNING: Joined gene families file not found", level='warning')
        return None
    
    run_cmd([
        "humann_split_stratified_table",
        "-i", joined_output,
        "-o", gene_families_out
    ])
    
    # Locate unstratified
    unstrat_file = None
    base_name = os.path.basename(joined_output).replace('.tsv', '')
    possible_patterns = [
        f"{base_name}_unstratified.tsv",
        f"{base_name}_unstratified.tsv",
        f"{base_name}-cpm_unstratified.tsv"
    ]
    
    for pattern in possible_patterns:
        test_path = os.path.join(gene_families_out, pattern)
        if os.path.isfile(test_path):
            unstrat_file = test_path
            break
    
    if not unstrat_file:
        for root, dirs, files in os.walk(gene_families_out):
            for fname in files:
                if 'unstratified' in fname.lower():
                    unstrat_file = os.path.join(root, fname)
                    break
            if unstrat_file:
                break
    
    if not unstrat_file:
        log_print("WARNING: Could not locate unstratified gene families file", level='warning')
        return None
    
    return unstrat_file