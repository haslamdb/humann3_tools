# humann3_tools/analysis/differential_abundance_cli.py

import os
import sys
import argparse
import logging
import pandas as pd

from humann3_tools.logger import setup_logger, log_print
from humann3_tools.analysis.differential_abundance import (
    aldex2_like,
    ancom,
    ancom_bc,
    run_differential_abundance_analysis
)
from humann3_tools.analysis.statistical import kruskal_wallis_dunn

def parse_args():
    """Parse command line arguments for the differential abundance CLI."""
    parser = argparse.ArgumentParser(
        description="Run differential abundance analysis on HUMAnN3 output files"
    )
    
    # Required arguments
    parser.add_argument("--abundance-file", required=True, 
                       help="Path to the unstratified abundance file (pathway or gene family)")
    parser.add_argument("--metadata-file", required=True,
                       help="Path to sample metadata CSV file")
    
    # Analysis options
    parser.add_argument("--output-dir", default="./DifferentialAbundance",
                       help="Directory for output files")
    parser.add_argument("--feature-type", choices=["pathway", "gene"], default="pathway",
                       help="Type of features in the abundance file (pathway or gene)")
    parser.add_argument("--group-col", default="Group",
                       help="Column name in metadata for grouping samples")
    parser.add_argument("--sample-id-col", 
                       help="Column name in metadata for sample IDs (autodetected if not specified)")
    parser.add_argument("--methods", default="aldex2,ancom,ancom-bc",
                       help="Comma-separated list of methods to use (aldex2,ancom,ancom-bc,kruskal)")
    parser.add_argument("--exclude-unmapped", action="store_true",
                       help="Exclude unmapped features from analysis")
    parser.add_argument("--filter-groups",
                       help="Comma-separated list of group names to include in the analysis. "
                            "For ALDEx2, exactly 2 groups must be specified.")
    
    # Additional options
    parser.add_argument("--log-file", default=None,
                       help="Path to log file")
    parser.add_argument("--log-level", default="INFO", 
                      choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                      help="Logging level (default: INFO)")
    
    return parser.parse_args()

def main():
    """Main function to run differential abundance analysis."""
    args = parse_args()
    
    # Setup logging
    logger = setup_logger(log_file=args.log_file, 
                         log_level=getattr(logging, args.log_level.upper()))
    log_print("Starting differential abundance analysis", level="info")
    
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
    
    # Set the sample ID as index
    metadata_df = metadata_df.set_index(sample_id_col)
    
    # Parse methods
    methods = [m.strip().lower() for m in args.methods.split(',')]
    log_print(f"Running {', '.join(methods)} on {args.feature_type} data", level="info")
    
    # Check if group column exists
    if args.group_col not in metadata_df.columns:
        log_print(f"ERROR: Group column '{args.group_col}' not found in metadata", level="error")
        log_print(f"Available columns: {', '.join(metadata_df.columns)}", level="info")
        sys.exit(1)
    
    # Process filter groups if specified
    filter_groups = None
    if args.filter_groups:
        filter_groups = [g.strip() for g in args.filter_groups.split(',')]
        log_print(f"Filtering to groups: {filter_groups}", level="info")
        
        # Verify these groups exist in the data
        unique_groups = metadata_df[args.group_col].unique()
        missing_groups = [g for g in filter_groups if g not in unique_groups]
        if missing_groups:
            log_print(f"ERROR: The following specified groups don't exist in the data: {missing_groups}", level="error")
            log_print(f"Available groups: {unique_groups}", level="info")
            sys.exit(1)
        
        # Check if ALDEx2 is selected and we have exactly 2 groups
        if "aldex2" in methods and len(filter_groups) != 2:
            log_print(f"WARNING: ALDEx2 requires exactly 2 groups, but {len(filter_groups)} were specified", level="warning")
            if "aldex2" in methods and len(methods) == 1:
                log_print("ERROR: Cannot proceed with ALDEx2 analysis without exactly 2 groups", level="error")
                sys.exit(1)
    
    # Create feature-specific output directory
    feature_dir = os.path.join(args.output_dir, args.feature_type.capitalize() + "s")
    os.makedirs(feature_dir, exist_ok=True)
    
    # Run Kruskal-Wallis test if selected
    if "kruskal" in methods:
        log_print("Running Kruskal-Wallis + Dunn's test...", level="info")
        
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
            metadata_df.reset_index(),
            on=sample_id_col,
            how="inner"
        )
        
        # Apply filtering if needed
        if filter_groups:
            long_df = long_df[long_df[args.group_col].isin(filter_groups)]
            if len(long_df) == 0:
                log_print(f"ERROR: No samples left after filtering to groups: {filter_groups}", level="error")
                sys.exit(1)
        
        # Run Kruskal-Wallis and Dunn's test
        kw_results, dunn_results = kruskal_wallis_dunn(
            df_long=long_df,
            group_col=args.group_col,
            feature_col=feature_col,
            abundance_col="Abundance",
            alpha=0.05,
            logger=logger
        )
        
        # Save results
        if not kw_results.empty:
            kw_path = os.path.join(feature_dir, "kruskal_wallis_results.csv")
            kw_results.to_csv(kw_path, index=False)
            log_print(f"Saved Kruskal-Wallis results to {kw_path}", level="info")
            
            sig_count = sum(kw_results["Reject_H0"])
            log_print(f"{sig_count} significant features found after FDR correction", level="info")
            
            if dunn_results:
                dunn_dir = os.path.join(feature_dir, "dunn_posthoc_tests")
                os.makedirs(dunn_dir, exist_ok=True)
                
                for feat, pdf in dunn_results.items():
                    feat_safe = feat.replace("/", "_").replace(" ", "_")
                    dunn_path = os.path.join(dunn_dir, f"dunn_{feat_safe}.csv")
                    pdf.to_csv(dunn_path)
                
                log_print(f"Saved Dunn's post-hoc results for {len(dunn_results)} features", level="info")
        else:
            log_print("No valid Kruskal-Wallis results produced", level="warning")
    
    # Remove kruskal from methods if present to avoid confusion with run_differential_abundance_analysis
    run_methods = [m for m in methods if m != "kruskal"]
    
    # Run other differential abundance methods if selected
    if run_methods:
        denom = "unmapped_excluded" if args.exclude_unmapped else "all"
        
        # Run differential abundance analysis
        results = run_differential_abundance_analysis(
            abundance_df=abundance_df,
            metadata_df=metadata_df,
            output_dir=feature_dir,
            group_col=args.group_col,
            methods=run_methods,
            denom=denom,
            filter_groups=filter_groups,
            logger=logger
        )
        
        if results:
            log_print("Differential abundance analysis completed successfully", level="info")
            # Print summary of significant features
            for method, result_df in results.items():
                if method == "aldex2":
                    sig_count = sum(result_df["q_value"] < 0.05)
                    log_print(f"ALDEx2: {sig_count} significant features (q < 0.05)", level="info")
                elif method == "ancom":
                    sig_count = sum(result_df["significant"])
                    log_print(f"ANCOM: {sig_count} significant features", level="info")
                elif method == "ancom_bc":
                    sig_count = sum(result_df["q_value"] < 0.05)
                    log_print(f"ANCOM-BC: {sig_count} significant features (q < 0.05)", level="info")
        else:
            log_print("Differential abundance analysis did not produce valid results", level="warning")
    
    log_print("Differential abundance analysis complete", level="info")
    return 0

if __name__ == "__main__":
    sys.exit(main())