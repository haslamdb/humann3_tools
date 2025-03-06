# humann3_tools/analysis/visualize_cli.py

import os
import argparse
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import sys

from humann3_tools.logger import setup_logger, log_print

def main():
    """Main function to create visualizations."""
    args = parse_args()
    
    # Setup logging
    logger = setup_logger(log_file=args.log_file, 
                         log_level=getattr(logging, args.log_level.upper()))
    log_print("Starting visualization generation", level="info")
    
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
    
    # Set the sample ID as index for metadata
    metadata_df = metadata_df.set_index(sample_id_col)
    
    # Check if group column exists
    if args.group_col not in metadata_df.columns:
        log_print(f"ERROR: Group column '{args.group_col}' not found in metadata", level="error")
        log_print(f"Available columns: {', '.join(metadata_df.columns)}", level="info")
        sys.exit(1)
    
    # Generate requested plots
    if args.pca:
        generate_pca_plot(abundance_df, metadata_df, args, logger)
    
    if args.barplot:
        generate_barplot(abundance_df, metadata_df, args, logger)
    
    if args.heatmap:
        generate_heatmap(abundance_df, metadata_df, args, logger)
    
    if args.abundance_hist:
        generate_abundance_histogram(abundance_df, metadata_df, args, logger)
    
    # Handle feature boxplot if specified
    if args.feature:
        generate_feature_boxplot(abundance_df, metadata_df, args, logger)
    
    # Handle top features boxplots if specified
    if args.box_top_n > 0:
        generate_top_features_boxplots(abundance_df, metadata_df, args, logger)
    
    log_print("Visualization generation complete", level="info")
    return 0

def generate_feature_boxplot(abundance_df, metadata_df, args, logger):
    """Generate boxplot for a specific feature."""
    feature = args.feature
    logger.info(f"Generating boxplot for feature: {feature}")
    
    # Check if feature exists in abundance data
    if feature not in abundance_df.index:
        logger.error(f"Feature '{feature}' not found in abundance data")
        logger.info(f"Available features include: {', '.join(abundance_df.index[:5])}...")
        return False
    
    # Get intersection of samples
    shared_samples = list(set(abundance_df.columns) & set(metadata_df.index))
    
    if not shared_samples:
        logger.error("No shared samples between abundance data and metadata")
        return False
    
    # Filter data to shared samples
    abundance = abundance_df.loc[feature, shared_samples]
    
    # Apply log transformation if requested
    if args.log_transform:
        abundance = np.log10(abundance + 1)
        y_label = f"log10({feature} + 1)"
        logger.info("Applied log10(x+1) transformation")
    else:
        y_label = feature
    
    # Create DataFrame for plotting
    plot_df = pd.DataFrame({
        'Abundance': abundance,
        args.group_col: [metadata_df.loc[sample, args.group_col] for sample in shared_samples]
    })
    
    # Create plot
    plt.figure(figsize=(10, 6))
    
    # Draw boxplot
    sns.boxplot(x=args.group_col, y='Abundance', data=plot_df)
    
    # Add individual points
    sns.stripplot(x=args.group_col, y='Abundance', data=plot_df, 
                 color='black', alpha=0.5, jitter=True)
    
    # Set labels and title
    plt.xlabel(args.group_col)
    plt.ylabel(y_label)
    
    feature_type = "Pathway" if args.feature_type == "pathway" else "Gene"
    plt.title(f"{feature_type} Abundance: {feature}")
    
    # Use tight layout to ensure everything fits
    plt.tight_layout()
    
    # Save plot
    feature_name_safe = feature.replace('/', '_').replace('\\', '_').replace(' ', '_')
    output_file = os.path.join(args.output_dir, f"boxplot_{feature_name_safe}.{args.format}")
    plt.savefig(output_file, dpi=args.dpi, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Boxplot saved to {output_file}")
    return True

def generate_top_features_boxplots(abundance_df, metadata_df, args, logger):
    """Generate boxplots for top N features."""
    logger.info(f"Generating boxplots for top {args.box_top_n} features...")
    
    # Get intersection of samples
    shared_samples = list(set(abundance_df.columns) & set(metadata_df.index))
    
    if not shared_samples:
        logger.error("No shared samples between abundance data and metadata")
        return False
    
    # Filter data to shared samples
    abundance = abundance_df[shared_samples]
    
    # Apply log transformation if requested for the selection of top features
    if args.log_transform:
        abundance_transformed = np.log10(abundance + 1)
    else:
        abundance_transformed = abundance
    
    # Calculate feature means and select top features
    feature_means = abundance_transformed.mean(axis=1)
    top_features = feature_means.sort_values(ascending=False).head(args.box_top_n).index
    
    # Create boxplot output directory
    boxplot_dir = os.path.join(args.output_dir, f"top{args.box_top_n}_boxplots")
    os.makedirs(boxplot_dir, exist_ok=True)
    
    # Generate boxplot for each top feature
    feature_type = "Pathway" if args.feature_type == "pathway" else "Gene"
    
    for i, feature in enumerate(top_features):
        logger.info(f"Generating boxplot for top feature {i+1}/{len(top_features)}: {feature}")
        
        # Extract feature abundance
        feature_abundance = abundance.loc[feature, shared_samples]
        
        # Apply log transformation if requested
        if args.log_transform:
            feature_abundance = np.log10(feature_abundance + 1)
            y_label = f"log10(Abundance + 1)"
        else:
            y_label = "Abundance"
        
        # Create DataFrame for plotting
        plot_df = pd.DataFrame({
            'Abundance': feature_abundance,
            args.group_col: [metadata_df.loc[sample, args.group_col] for sample in shared_samples]
        })
        
        # Create plot
        plt.figure(figsize=(10, 6))
        
        # Draw boxplot
        sns.boxplot(x=args.group_col, y='Abundance', data=plot_df)
        
        # Add individual points
        sns.stripplot(x=args.group_col, y='Abundance', data=plot_df, 
                     color='black', alpha=0.5, jitter=True)
        
        # Set labels and title
        plt.xlabel(args.group_col)
        plt.ylabel(y_label)
        
        plt.title(f"Top {i+1}: {feature}")
        
        # Use tight layout to ensure everything fits
        plt.tight_layout()
        
        # Save plot
        feature_name_safe = feature.replace('/', '_').replace('\\', '_').replace(' ', '_')
        output_file = os.path.join(boxplot_dir, f"boxplot_{i+1}_{feature_name_safe}.{args.format}")
        plt.savefig(output_file, dpi=args.dpi, bbox_inches='tight')
        plt.close()
    
    logger.info(f"Generated {len(top_features)} boxplots in {boxplot_dir}")
    return True

def generate_barplot(abundance_df, metadata_df, args, logger):
    """Generate barplot of top features by group."""
    logger.info("Generating barplot of top features...")
    
    # Get intersection of samples
    shared_samples = list(set(abundance_df.columns) & set(metadata_df.index))
    
    if not shared_samples:
        logger.error("No shared samples between abundance data and metadata")
        return False
    
    # Filter data to shared samples
    abundance = abundance_df[shared_samples]
    
    # Get group information for each sample
    sample_groups = metadata_df.loc[shared_samples, args.group_col]
    
    # Calculate mean abundance per group
    group_means = {}
    for group in sample_groups.unique():
        group_samples = sample_groups[sample_groups == group].index
        group_means[group] = abundance[group_samples].mean(axis=1)
    
    # Combine group means into a DataFrame
    mean_df = pd.DataFrame(group_means)
    
    # Calculate overall mean abundance for each feature
    mean_df['overall_mean'] = mean_df.mean(axis=1)
    
    # Sort by overall mean abundance and get top N features
    top_features = mean_df.sort_values('overall_mean', ascending=False).head(args.top_n).drop('overall_mean', axis=1)
    
    # Transpose for easier plotting
    plot_df = top_features.transpose()
    
    # Create bar plot
    plt.figure(figsize=(12, 8))
    plot_df.plot(kind='bar', ax=plt.gca())
    
    # Customize plot
    plt.xticks(rotation=45, ha='right')
    plt.xlabel('Group')
    plt.ylabel('Mean Abundance')
    
    feature_type = "Pathways" if args.feature_type == "pathway" else "Gene Families"
    plt.title(f'Top {args.top_n} {feature_type} by Mean Abundance')
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    # Save plot
    output_file = os.path.join(args.output_dir, f"barplot_top{args.top_n}_{args.feature_type}.{args.format}")
    plt.savefig(output_file, dpi=args.dpi, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Barplot saved to {output_file}")
    return True

def generate_heatmap(abundance_df, metadata_df, args, logger):
    """Generate heatmap of top features."""
    logger.info("Generating heatmap of top features...")
    
    # Get intersection of samples
    shared_samples = list(set(abundance_df.columns) & set(metadata_df.index))
    
    if not shared_samples:
        logger.error("No shared samples between abundance data and metadata")
        return False
    
    # Filter data to shared samples
    abundance = abundance_df[shared_samples]
    
    # Apply log transformation if requested
    if args.log_transform:
        abundance_transformed = np.log10(abundance + 1)
        logger.info("Applied log10(x+1) transformation")
    else:
        abundance_transformed = abundance
    
    # Calculate feature means and select top features
    feature_means = abundance_transformed.mean(axis=1)
    top_features = feature_means.sort_values(ascending=False).head(args.top_n).index
    
    # Filter to top features
    top_data = abundance_transformed.loc[top_features]
    
    # Get sample grouping
    sample_groups = metadata_df.loc[shared_samples, args.group_col]
    
    # Sort samples by group
    sorted_samples = sample_groups.sort_values().index
    
    # Prepare data for heatmap
    heatmap_data = top_data[sorted_samples]
    
    # Create row colors based on groups
    group_colors = sns.color_palette("husl", n_colors=len(sample_groups.unique()))
    group_color_map = dict(zip(sample_groups.unique(), group_colors))
    col_colors = pd.Series(sample_groups).map(group_color_map)
    
    # Create heatmap plot
    plt.figure(figsize=(14, 10))
    g = sns.clustermap(
        heatmap_data,
        col_cluster=False,
        yticklabels=True,
        cmap="viridis",
        col_colors=col_colors,
        figsize=(14, 10)
    )
    
    # Create legend for groups
    for group, color in group_color_map.items():
        g.ax_heatmap.bar(0, 0, color=color, label=group, alpha=0.8)
    g.ax_heatmap.legend(title=args.group_col, loc="center left", bbox_to_anchor=(1, 0.5))
    
    # Set title
    feature_type = "Pathways" if args.feature_type == "pathway" else "Gene Families"
    plt.suptitle(f'Heatmap of Top {args.top_n} {feature_type}', y=1.02)
    
    # Save plot
    output_file = os.path.join(args.output_dir, f"heatmap_top{args.top_n}_{args.feature_type}.{args.format}")
    plt.savefig(output_file, dpi=args.dpi, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Heatmap saved to {output_file}")
    return True

def generate_abundance_histogram(abundance_df, metadata_df, args, logger):
    """Generate histograms of abundance distributions."""
    logger.info("Generating abundance histograms...")
    
    # Get intersection of samples
    shared_samples = list(set(abundance_df.columns) & set(metadata_df.index))
    
    if not shared_samples:
        logger.error("No shared samples between abundance data and metadata")
        return False
    
    # Filter data to shared samples
    abundance = abundance_df[shared_samples]
    
    # Apply log transformation if requested
    if args.log_transform:
        abundance_transformed = np.log10(abundance + 1)
        transform_label = "log10(Abundance + 1)"
        logger.info("Applied log10(x+1) transformation")
    else:
        abundance_transformed = abundance
        transform_label = "Abundance"
    
    # Get group information
    sample_groups = metadata_df.loc[shared_samples, args.group_col]
    
    # Create plot
    plt.figure(figsize=(12, 8))
    
    # Plot histograms for each group
    for group in sample_groups.unique():
        group_samples = sample_groups[sample_groups == group].index
        group_data = abundance_transformed[group_samples].values.flatten()
        
        # Filter out zeros and NaNs
        group_data = group_data[~np.isnan(group_data)]
        group_data = group_data[group_data > 0]
        
        sns.histplot(group_data, kde=True, label=group, alpha=0.5)
    
    # Set plot attributes
    plt.xlabel(transform_label)
    plt.ylabel('Count')
    
    feature_type = "Pathways" if args.feature_type == "pathway" else "Gene Families"
    plt.title(f'Distribution of {feature_type} Abundance by Group')
    
    plt.legend(title=args.group_col)
    plt.tight_layout()
    
    # Save plot
    output_file = os.path.join(args.output_dir, f"histogram_{args.feature_type}.{args.format}")
    plt.savefig(output_file, dpi=args.dpi, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Abundance histogram saved to {output_file}")
    return True

def parse_args():
    """Parse command line arguments for the visualization CLI."""
    parser = argparse.ArgumentParser(
        description="Create visualizations for HUMAnN3 output files"
    )
    
    # Required arguments
    parser.add_argument("--abundance-file", required=True, 
                      help="Path to the unstratified abundance file (pathway or gene family)")
    parser.add_argument("--metadata-file", required=True,
                      help="Path to sample metadata CSV file")
    
    # Visualization options
    parser.add_argument("--output-dir", default="./Visualizations",
                      help="Directory for output files")
    parser.add_argument("--feature-type", choices=["pathway", "gene"], default="pathway",
                      help="Type of features in the abundance file (pathway or gene)")
    parser.add_argument("--group-col", default="Group",
                      help="Column name in metadata for coloring points")
    parser.add_argument("--shape-col", 
                      help="Column name in metadata for point shapes")
    parser.add_argument("--sample-id-col", 
                      help="Column name in metadata for sample IDs (autodetected if not specified)")
    parser.add_argument("--top-n", type=int, default=25,
                      help="Number of top features to include in bar plots (default: 25)")
    parser.add_argument("--format", default="svg", choices=["svg", "png", "pdf"],
                      help="Output format for plots (default: svg)")
    parser.add_argument("--dpi", type=int, default=300,
                      help="DPI for raster formats like PNG (default: 300)")
    
    # Plot selection
    parser.add_argument("--pca", action="store_true", default=True,
                      help="Generate PCA plot")
    parser.add_argument("--heatmap", action="store_true",
                      help="Generate heatmap of top features")
    parser.add_argument("--barplot", action="store_true", default=True,
                      help="Generate barplot of top features by group")
    parser.add_argument("--abundance-hist", action="store_true",
                      help="Generate histograms of abundance distributions")
    
    # Boxplot specific options
    parser.add_argument("--feature", 
                      help="Generate boxplot for a specific feature (pathway or gene)")
    parser.add_argument("--box-top-n", type=int, default=0,
                      help="Generate boxplots for top N features (default: 0, disabled)")
    
    # Additional options
    parser.add_argument("--log-transform", action="store_true", default=True,
                      help="Apply log10(x+1) transformation to abundance data")
    parser.add_argument("--log-file", default=None,
                      help="Path to log file")
    parser.add_argument("--log-level", default="INFO", 
                      choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                      help="Logging level (default: INFO)")
    
    return parser.parse_args()

def generate_pca_plot(abundance_df, metadata_df, args, logger):
    """Generate PCA plot from abundance data."""
    logger.info("Generating PCA plot...")
    
    # Get intersection of samples
    shared_samples = list(set(abundance_df.columns) & set(metadata_df.index))
    
    if not shared_samples:
        logger.error("No shared samples between abundance data and metadata")
        return False
    
    # Filter data to shared samples
    abundance = abundance_df[shared_samples]
    
    # Apply log transformation if requested
    if args.log_transform:
        abundance_transformed = np.log10(abundance + 1)
        logger.info("Applied log10(x+1) transformation")
    else:
        abundance_transformed = abundance
    
    # Scale data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(abundance_transformed.T)
    
    # Run PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_data)
    
    # Create DataFrame with PCA results
    pca_df = pd.DataFrame(
        data=pca_result,
        columns=['PC1', 'PC2'],
        index=shared_samples
    )
    
    # Merge with metadata
    pca_merged = pca_df.join(metadata_df[metadata_df.index.isin(shared_samples)])
    
    # Create plot
    plt.figure(figsize=(10, 8))
    
    if args.shape_col and args.shape_col in metadata_df.columns:
        # Use both color and shape for grouping
        sns.scatterplot(
            data=pca_merged,
            x='PC1',
            y='PC2',
            hue=args.group_col,
            style=args.shape_col,
            s=100,
            alpha=0.8
        )
    else:
        # Use only color for grouping
        sns.scatterplot(
            data=pca_merged,
            x='PC1',
            y='PC2',
            hue=args.group_col,
            s=100,
            alpha=0.8
        )
    
    # Add variance explained
    variance_explained = pca.explained_variance_ratio_ * 100
    plt.xlabel(f'PC1 ({variance_explained[0]:.1f}%)')
    plt.ylabel(f'PC2 ({variance_explained[1]:.1f}%)')
    
    # Add title
    feature_type = "Pathways" if args.feature_type == "pathway" else "Gene Families"
    plt.title(f'PCA of {feature_type}')
    
    # Move legend outside the plot if needed
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    # Save plot
    output_file = os.path.join(args.output_dir, f"pca_{args.feature_type}.{args.format}")
    plt.savefig(output_file, dpi=args.dpi, bbox_inches='tight')
    plt.close()
    
    logger.info(f"PCA plot saved to {output_file}")
    return True

if __name__ == "__main__":
    sys.exit(main())