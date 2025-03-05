# humann3_tools/analysis/statistical_cli.py

import os
import argparse
import logging
import pandas as pd
import sys

from humann3_tools.logger import setup_logger, log_print
from humann3_tools.analysis.statistical import kruskal_wallis_dunn

def parse_args():
    """Parse command line arguments for the statistical testing CLI."""
    parser = argparse.ArgumentParser(
        description="Run statistical tests on HUMAnN3 output files"
    )
    
    # Required arguments
    parser.add_argument("--abundance-file", required=True, 
                      help="Path to the unstratified abundance file (pathway or gene family)")
    parser.add_argument("--metadata-file", required=True,
                      help="Path to sample metadata CSV file")
    
    # Analysis options
    parser.add_argument("--output-dir", default="./StatisticalTests",
                      help="Directory for output files")
    parser.add_argument("--feature-type", choices=["pathway", "gene"], default="pathway",
                      help="Type of features in the abundance file (pathway or gene)")
    parser.add_argument("--group-col", default="Group",
                      help="Column name in metadata for grouping samples")
    parser.add_argument("--sample-id-col", 
                      help="Column name in metadata for sample IDs (autodetected if not specified)")
    parser.add_argument("--alpha", type=float, default=0.05,
                      help="Significance threshold for statistical tests (default: 0.05)")
    
    # Additional options
    parser.add_argument("--log-file", default=None,
                      help="Path to log file")
    parser.add_argument("--log-level", default="INFO", 
                      choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                      help="Logging level (default: INFO)")
    
    return parser.parse_args()

def main():
    """Main function to run statistical tests."""
    args = parse_args()
    
    # Setup logging
    logger = setup_logger(log_file=args.log_file, 
                         log_level=getattr(logging, args.log_level.upper()))
    log_print("Starting statistical analysis", level="info")
    
    # Check if files exist
    if not os.path.exists(args.abundance_file):
        log_print(f"ERROR: Abundance file not found: {args.abundance_file}", level="error")
        sys.exit(1)
    
    if not os.path.exists(args.metadata_file):
        log_print(f"ERROR: Metadata file not found: {args.metadata_file}", level="error")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Read data
    try:
        abundance_df = pd.read_csv(args.abundance_file, sep="\t", index_col=0)
        log_print(f"Loaded abundance data with {abundance_df.shape[0]} features and {abundance_df.shape[1]} samples", 
                 level="info")
    except Exception as e:
        log_print(f"ERROR reading abundance file: {str(e)}", level="error")
        sys.exit(1)
    
    try:
        metadata_df = pd.read_csv(args.metadata_file)
        log_print(f"Loaded metadata with {metadata_df.shape[0]} samples and {metadata_df.shape[1]} columns", 
                 level="info")
    except Exception as e:
        log_print(f"ERROR reading metadata file: {str(e)}", level="error")
        sys.exit(1)
    
    # Find sample ID column in metadata if not specified
    sample_id_col = args.sample_id_col
    if not sample_id_col:
        # Try to auto-detect
        common_id_cols = ["SampleName", "Sample", "SampleID", "Sample_ID", "sample_name", "sample_id"]
        for col in common_id_cols:
            if col in metadata_df.columns:
                sample_id_col = col
                log_print(f"Auto-detected sample ID column: {sample_id_col}", level="info")
                break
        
        if not sample_id_col:
            log_print("ERROR: Could not auto-detect sample ID column in metadata", level="error")
            log_print(f"Available columns: {', '.join(metadata_df.columns)}", level="info")
            sys.exit(1)
    
    # Check if group column exists
    if args.group_col not in metadata_df.columns:
        log_print(f"ERROR: Group column '{args.group_col}' not found in metadata", level="error")
        log_print(f"Available columns: {', '.join(metadata_df.columns)}", level="info")
        sys.exit(1)
    
    # Convert to long format (needed for kruskal_wallis_dunn)
    long_df = abundance_df.reset_index().melt(
        id_vars=abundance_df.index.name, 
        var_name=sample_id_col, 
        value_name="Abundance"
    )
    
    # Rename feature column appropriately
    feature_col = "Pathway" if args.feature_type == "pathway" else "Gene_Family"
    long_df = long_df.rename(columns={abundance_df.index.name: feature_col})
    
    # Merge with metadata
    long_df = pd.merge(
        long_df, 
        metadata_df,
        on=sample_id_col,
        how="inner"
    )
    
    # Check if merge was successful
    if long_df.empty:
        log_print("ERROR: No matching samples between abundance data and metadata", level="error")
        sys.exit(1)
    
    log_print(f"Successfully merged abundance data with metadata. Working with {len(long_df)} rows.", level="info")
    
    # Run Kruskal-Wallis and Dunn's test
    kw_results, dunn_results = kruskal_wallis_dunn(
        df_long=long_df,
        group_col=args.group_col,
        feature_col=feature_col,
        abundance_col="Abundance",
        alpha=args.alpha,
        logger=logger
    )
    
    # Save results
    if not kw_results.empty:
        kw_path = os.path.join(args.output_dir, "kruskal_wallis_results.csv")
        kw_results.to_csv(kw_path, index=False)
        log_print(f"Saved Kruskal-Wallis results to {kw_path}", level="info")
        
        sig_count = sum(kw_results["Reject_H0"])
        log_print(f"{sig_count} significant features found after FDR correction", level="info")
        
        if dunn_results:
            dunn_dir = os.path.join(args.output_dir, "dunn_posthoc_tests")
            os.makedirs(dunn_dir, exist_ok=True)
            
            for feat, pdf in dunn_results.items():
                feat_safe = feat.replace("/", "_").replace(" ", "_")
                dunn_path = os.path.join(dunn_dir, f"dunn_{feat_safe}.csv")
                pdf.to_csv(dunn_path)
            
            log_print(f"Saved Dunn's post-hoc results for {len(dunn_results)} features", level="info")
    else:
        log_print("No valid Kruskal-Wallis results produced", level="warning")
    
    log_print("Statistical testing complete", level="info")
    return 0

if __name__ == "__main__":
    sys.exit(main())
