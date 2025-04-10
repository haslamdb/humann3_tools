def parse_args():
    """Parse command line arguments for the Differential Abundance module."""
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
                      help="Comma-separated list of methods to use (aldex2,ancom,ancom-bc)")
    parser.add_argument("--exclude-unmapped", action="store_true",
                      help="Exclude unmapped features from analysis")
    parser.add_argument("--filter-groups",
                      help="Comma-separated list of group names to include in the analysis. "
                            "For ALDEx2, exactly 2 groups must be specified.")
    parser.add_argument("--alpha", type=float, default=0.05,
                      help="Significance threshold for statistical tests (default: 0.05)")
    
    # Additional options
    parser.add_argument("--log-file", 
                      help="Path to log file")
    parser.add_argument("--log-level", default="INFO", 
                      choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                      help="Logging level")
    
    return parser.parse_args()

@track_peak_memory
def main():
    """Main function to run differential abundance analysis."""
    # Parse arguments
    args = parse_args()
    
    # Setup logging
    log_level = getattr(logging, args.log_level.upper())
    setup_logger(args.log_file, log_level)
    
    logger.info("Starting HUMAnN3 Tools Differential Abundance Module")
    start_time = time.time()
    
    # Process methods
    methods = [m.strip().lower() for m in args.methods.split(',')]
    logger.info(f"Using methods: {', '.join(methods)}")
    
    # Process filter groups
    filter_groups = None
    if args.filter_groups:
        filter_groups = [g.strip() for g in args.filter_groups.split(',')]
        logger.info(f"Filtering to groups: {filter_groups}")
    
    # Check if files exist
    for file_path, desc in [(args.abundance_file, "Abundance file"), (args.metadata_file, "Metadata file")]:
        if not os.path.exists(file_path):
            logger.error(f"ERROR: {desc} not found: {file_path}")
            return 1
    
    # Set denominator based on exclude_unmapped flag
    denom = "unmapped_excluded" if args.exclude_unmapped else "all"
    
    # Run differential abundance analysis
    success = run_differential_abundance_analysis(
        abundance_file=args.abundance_file,
        metadata_file=args.metadata_file,
        output_dir=args.output_dir,
        group_col=args.group_col,
        methods=methods,
        feature_type=args.feature_type,
        denom=denom,
        filter_groups=filter_groups,
        sample_id_col=args.sample_id_col,
        alpha=args.alpha
    )
    
    if not success:
        logger.error("Differential abundance analysis failed")
        return 1
    
    # Print elapsed time
    elapsed_time = time.time() - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    
    logger.info(f"Total processing time: {int(minutes)}m {int(seconds)}s")
    
    # Print next steps
    logger.info("\nNext Steps:")
    logger.info("  For visualizations:")
    logger.info(f"  humann3-tools viz --abundance-file {args.abundance_file} --metadata-file {args.metadata_file}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())    """
    A Python implementation similar to ALDEx2 for differential abundance testing
    
    Parameters:
        abundance_df: Feature table with samples as columns and features as rows
        metadata_df: Metadata with sample IDs as index and metadata as columns
        group_col: Column name in metadata_df that contains the grouping variable
        mc_samples: Number of Monte Carlo samples to generate
        denom: Features to use as denominator: "all" for all features, 
               "unmapped_excluded" to exclude unmapped features
        filter_groups: List of group names to include in the analysis. If provided, 
                      only these groups will be used. Must contain exactly 2 groups for ALDEx2.
        
    Returns:
        pandas DataFrame with test results
    """
    # Make sure metadata and abundance data have matching samples
    shared_samples = list(set(abundance_df.columns) & set(metadata_df.index))
    if len(shared_samples) == 0:
        logger.error("No shared samples between abundance data and metadata")
        return pd.DataFrame()
    
    abundance = abundance_df[shared_samples].copy()
    metadata = metadata_df.loc[shared_samples].copy()
    
    # Handle unmapped reads based on denom parameter
    if denom == "unmapped_excluded" and "UNMAPPED" in abundance.index:
        logger.info("Excluding unmapped reads from denominator")
        abundance = abundance.drop("UNMAPPED", axis=0)
    
    # Replace zeros with small value (pseudocount)
    min_val = abundance[abundance > 0].min().min() / 2
    abundance = abundance.replace(0, min_val)
    
    # Get group information
    groups = metadata[group_col]
    unique_groups = groups.unique()
    
    # Apply group filtering if specified
    if filter_groups is not None:
        if not isinstance(filter_groups, list):
            filter_groups = [filter_groups]
        
        # Verify the specified groups exist in the data
        missing_groups = [g for g in filter_groups if g not in unique_groups]
        if missing_groups:
            logger.error(f"The following specified groups don't exist in the data: {missing_groups}")
            return pd.DataFrame()
        
        # Filter metadata and abundance to only include samples from specified groups
        valid_samples = metadata[metadata[group_col].isin(filter_groups)].index
        if len(valid_samples) == 0:
            logger.error(f"No samples found for groups: {filter_groups}")
            return pd.DataFrame()
        
        metadata = metadata.loc[valid_samples]
        abundance = abundance[valid_samples]
        groups = metadata[group_col]
        unique_groups = groups.unique()
        
        logger.info(f"Filtered to {len(unique_groups)} groups: {unique_groups}")
    
    # ALDEx2 requires exactly 2 groups for comparison
    if len(unique_groups) != 2:
        logger.error(f"ALDEx2 implementation requires exactly 2 groups for comparison, found {len(unique_groups)}: {unique_groups}")
        if filter_groups:
            logger.error(f"Please filter to exactly 2 groups using --filter-groups option")
        else:
            logger.error(f"Please use --filter-groups option to select 2 groups from: {unique_groups}")
        return pd.DataFrame()
        
    group1_samples = groups[groups == unique_groups[0]].index
    group2_samples = groups[groups == unique_groups[1]].index
    
    logger.info(f"Running ALDEx2 analysis with {len(group1_samples)} samples in group '{unique_groups[0]}' and "
                f"{len(group2_samples)} samples in group '{unique_groups[1]}'")
    
    # Results dataframe
    results = pd.DataFrame(index=abundance.index)
    results['feature'] = abundance.index
    
    # Monte Carlo sampling of CLR transformations
    all_clrs = []
    
    for i in range(mc_samples):
        # Generate Dirichlet Monte Carlo instance
        mc_instance = pd.DataFrame(index=abundance.index, columns=abundance.columns)
        
        for col in abundance.columns:
            # Add random noise according to Dirichlet distribution
            mc_instance[col] = np.random.dirichlet(abundance[col], 1)[0] * abundance[col].sum()
        
        # CLR transformation
        clr_data = clr_transform(mc_instance.T + 0.5).T
        all_clrs.append(clr_data)
    
    # Calculate effect sizes and p-values across MC instances
    effect_sizes = []
    pvals = []
    
    for clr_data in all_clrs:
        # For each feature, compare between groups
        for feature in clr_data.index:
            group1_values = clr_data.loc[feature, group1_samples]
            group2_values = clr_data.loc[feature, group2_samples]
            
            # Calculate effect size (difference of means)
            effect = group1_values.mean() - group2_values.mean()
            effect_sizes.append((feature, effect))
            
            # Welch's t-test
            t_stat, p_val = stats.ttest_ind(
                group1_values, group2_values, equal_var=False
            )
            pvals.append((feature, p_val))
    
    # Aggregate results
    effect_df = pd.DataFrame(effect_sizes, columns=['feature', 'effect'])
    pval_df = pd.DataFrame(pvals, columns=['feature', 'pval'])
    
    # Group by feature and calculate median effect and p-value
    median_effects = effect_df.groupby('feature')['effect'].median()
    median_pvals = pval_df.groupby('feature')['pval'].median()
    
    # Add to results
    results['effect_size'] = median_effects
    results['p_value'] = median_pvals
    
    # Multiple testing correction
    results['q_value'] = multipletests(results['p_value'], method='fdr_bh')[1]
    
    # Add mean abundance information
    results['mean_abundance_group1'] = abundance[group1_samples].mean(axis=1)
    results['mean_abundance_group2'] = abundance[group2_samples].mean(axis=1)
    
    return results.sort_values('q_value')

def ancom(
    abundance_df: pd.DataFrame, 
    metadata_df: pd.DataFrame, 
    group_col: str, 
    alpha: float = 0.05, 
    denom: str = "all", 
    filter_groups: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    ANCOM for differential abundance testing
    
    Parameters:
        abundance_df: Feature table with samples as columns and features as rows
        metadata_df: Metadata with sample IDs as index and metadata as columns
        group_col: Column name in metadata_df that contains the grouping variable
        alpha: Significance level for tests
        denom: Features to use as denominator: "all" for all features, 
               "unmapped_excluded" to exclude unmapped features
        filter_groups: List of group names to include in the analysis.
        
    Returns:
        pandas DataFrame with test results
    """
    # Make sure metadata and abundance data have matching samples
    shared_samples = list(set(abundance_df.columns) & set(metadata_df.index))
    if len(shared_samples) == 0:
        logger.error("No shared samples between abundance data and metadata")
        return pd.DataFrame()
    
    abundance = abundance_df[shared_samples].copy()
    metadata = metadata_df.loc[shared_samples].copy()
    
    # Handle unmapped reads based on denom parameter
    if denom == "unmapped_excluded" and "UNMAPPED" in abundance.index:
        logger.info("Excluding unmapped reads from denominator")
        abundance = abundance.drop("UNMAPPED", axis=0)
    
    # Replace zeros with small value (pseudocount)
    min_val = 0.5
    abundance = abundance.replace(0, min_val)
    
    # Get group information
    groups = metadata[group_col]
    unique_groups = groups.unique()
    
    # Apply group filtering if specified
    if filter_groups is not None:
        if not isinstance(filter_groups, list):
            filter_groups = [filter_groups]
        
        # Verify the specified groups exist in the data
        missing_groups = [g for g in filter_groups if g not in unique_groups]
        if missing_groups:
            logger.error(f"The following specified groups don't exist in the data: {missing_groups}")
            return pd.DataFrame()
        
        # Filter metadata and abundance to only include samples from specified groups
        valid_samples = metadata[metadata[group_col].isin(filter_groups)].index
        if len(valid_samples) == 0:
            logger.error(f"No samples found for groups: {filter_groups}")
            return pd.DataFrame()
        
        metadata = metadata.loc[valid_samples]
        abundance = abundance[valid_samples]
        groups = metadata[group_col]
        unique_groups = groups.unique()
        
        logger.info(f"Filtered to {len(unique_groups)} groups: {unique_groups}")
    
    # ANCOM needs at least 2 groups
    if len(unique_groups) < 2:
        logger.error(f"ANCOM requires at least 2 groups for comparison, found {len(unique_groups)}")
        return pd.DataFrame()
    
    logger.info(f"Running ANCOM analysis with {len(unique_groups)} groups: {unique_groups}")
    
    # Results dataframe
    results = pd.DataFrame(index=abundance.index)
    results['feature'] = abundance.index
    
    # For each feature, compare against all other features
    n_features = len(abundance.index)
    W = {feature: 0 for feature in abundance.index}
    
    # For each feature (i)
    for i, feature_i in enumerate(abundance.index):
        if i % 50 == 0:  # Progress update for large datasets
            logger.info(f"ANCOM: Processing feature {i+1}/{n_features}")
            
        # Compare with all other features (j)
        for j, feature_j in enumerate(abundance.index):
            if i != j:
                # Log ratio
                log_ratio = np.log(abundance.loc[feature_i] / abundance.loc[feature_j])
                
                # Perform statistical test on log ratio between groups
                group_values = {}
                for group in groups.unique():
                    group_values[group] = log_ratio[groups == group]
                
                if len(groups.unique()) == 2:
                    # Two groups: t-test
                    group_list = list(group_values.values())
                    t_stat, p_val = stats.ttest_ind(
                        group_list[0], group_list[1], equal_var=False
                    )
                else:
                    # Multiple groups: ANOVA
                    anova_groups = []
                    for group, values in group_values.items():
                        anova_groups.append(values)
                    f_stat, p_val = stats.f_oneway(*anova_groups)
                
                # Count significant tests
                if p_val < alpha:
                    W[feature_i] += 1
    
    # Calculate W statistic and add to results
    results['W'] = [W[feature] for feature in results['feature']]
    results['W_ratio'] = results['W'] / (n_features - 1)
    
    # Add detection threshold (cutoff)
    cutoff = 0.7  # Can be adjusted based on desired sensitivity
    results['significant'] = results['W_ratio'] > cutoff
    
    # Add mean abundance information
    for group in groups.unique():
        group_samples = groups[groups == group].index
        results[f'mean_abundance_{group}'] = abundance[group_samples].mean(axis=1)
    
    return results.sort_values('W', ascending=False)

def ancom_bc(
    abundance_df: pd.DataFrame, 
    metadata_df: pd.DataFrame, 
    group_col: str, 
    formula: Optional[str] = None, 
    denom: str = "all", 
    filter_groups: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    ANCOM-BC for differential abundance testing
    
    Parameters:
        abundance_df: Feature table with samples as columns and features as rows
        metadata_df: Metadata with sample IDs as index and metadata as columns
        group_col: Column name in metadata_df that contains the grouping variable
        formula: R-style formula for the model (e.g., "~ Group + Covariate")
                If None, will use simple one-way formula with group_col
        denom: Features to use as denominator: "all" for all features, 
               "unmapped_excluded" to exclude unmapped features
        filter_groups: List of group names to include in the analysis.
        
    Returns:
        pandas DataFrame with test results
    """
    try:
        import statsmodels.api as sm
        from statsmodels.formula.api import ols
    except ImportError:
        logger.error("statsmodels is required for ANCOM-BC")
        return pd.DataFrame()
    
    # Make sure metadata and abundance data have matching samples
    shared_samples = list(set(abundance_df.columns) & set(metadata_df.index))
    if len(shared_samples) == 0:
        logger.error("No shared samples between abundance data and metadata")
        return pd.DataFrame()
    
    abundance = abundance_df[shared_samples].copy()
    metadata = metadata_df.loc[shared_samples].copy()
    
    # Handle unmapped reads based on denom parameter
    if denom == "unmapped_excluded" and "UNMAPPED" in abundance.index:
        logger.info("Excluding unmapped reads from denominator")
        abundance = abundance.drop("UNMAPPED", axis=0)
    
    # Replace zeros with small value (pseudocount)
    min_val = 0.5
    abundance = abundance.replace(0, min_val)
    
    # Get group information
    groups = metadata[group_col]
    unique_groups = groups.unique()
    
    # Apply group filtering if specified
    if filter_groups is not None:
        if not isinstance(filter_groups, list):
            filter_groups = [filter_groups]
        
        # Verify the specified groups exist in the data
        missing_groups = [g for g in filter_groups if g not in unique_groups]
        if missing_groups:
            logger.error(f"The following specified groups don't exist in the data: {missing_groups}")
            return pd.DataFrame()
        
        # Filter metadata and abundance to only include samples from specified groups
        valid_samples = metadata[metadata[group_col].isin(filter_groups)].index
        if len(valid_samples) == 0:
            logger.error(f"No samples found for groups: {filter_groups}")
            return pd.DataFrame()
        
        metadata = metadata.loc[valid_samples]
        abundance = abundance[valid_samples]
        groups = metadata[group_col]
        unique_groups = groups.unique()
        
        logger.info(f"Filtered to {len(unique_groups)} groups: {unique_groups}")
    
    # ANCOM-BC needs at least 2 groups
    if len(unique_groups) < 2:
        logger.error(f"ANCOM-BC requires at least 2 groups for comparison, found {len(unique_groups)}")
        return pd.DataFrame()
    
    logger.info(f"Running ANCOM-BC analysis with {len(unique_groups)} groups: {unique_groups}")
    
    # Results dataframe
    results = pd.DataFrame(index=abundance.index)
    results['feature'] = abundance.index
    
    # Estimate sequencing depth using ANCOM-BC procedure
    # (This is a simplified version of the bias correction step)
    
    # 1. Log transform the data
    log_abundance = np.log(abundance)
    
    # 2. Get sample-wise means (log geometric means)
    sample_means = log_abundance.mean(axis=0)
    
    # 3. Center log-ratio transform to remove compositional bias
    clr_abundance = log_abundance.sub(sample_means, axis=1)
    
    # 4. Transpose for regression (samples as rows, features as columns)
    clr_abundance_t = clr_abundance.T
    
    # Set up the model formula if not provided
    if formula is None:
        formula = f"feature ~ C({group_col})"
        logger.info(f"Using formula: {formula}")
    
    # Run models for each feature
    pvals = []
    effects = []
    
    for i, feature in enumerate(clr_abundance_t.columns):
        if i % 100 == 0:  # Progress update for large datasets
            logger.info(f"ANCOM-BC: Processing feature {i+1}/{len(clr_abundance_t.columns)}")
            
        # Create a temporary dataframe for regression
        temp_df = pd.DataFrame({
            'feature': clr_abundance_t[feature],
            **metadata
        })
        
        try:
            # Fit the model
            model = ols(formula, data=temp_df).fit()
            
            # Extract p-values for the group effect
            for term in model.pvalues.index:
                if group_col in term and term != 'Intercept':
                    pvals.append((feature, model.pvalues[term]))
                    effects.append((feature, model.params[term]))
                    break
        except Exception as e:
            logger.warning(f"Error in ANCOM-BC regression for feature {feature}: {str(e)}")
    
    # Compile results
    pval_df = pd.DataFrame(pvals, columns=['feature', 'p_value'])
    effect_df = pd.DataFrame(effects, columns=['feature', 'effect_size'])
    
    # Add to results
    results = results.merge(pval_df, on='feature', how='left')
    results = results.merge(effect_df, on='feature', how='left')
    
    # Multiple testing correction
    results['q_value'] = multipletests(results['p_value'], method='fdr_bh')[1]
    
    # Add mean abundance information
    for group in groups.unique():
        group_samples = groups[groups == group].index
        results[f'mean_abundance_{group}'] = abundance[group_samples].mean(axis=1)
    
    return results.sort_values('q_value')

def run_differential_abundance_analysis(
    abundance_file: str,
    metadata_file: str,
    output_dir: str,
    group_col: str = "Group",
    methods: List[str] = None,
    feature_type: str = "pathway",
    denom: str = "all",
    filter_groups: Optional[List[str]] = None,
    sample_id_col: Optional[str] = None,
    alpha: float = 0.05
) -> bool:
    """
    Run differential abundance analysis using multiple methods.
    
    Args:
        abundance_file: Path to abundance file (unstratified)
        metadata_file: Path to metadata file
        output_dir: Directory for output files
        group_col: Column in metadata for grouping
        methods: List of methods to run ("aldex2", "ancom", "ancom-bc")
        feature_type: Type of features ("pathway" or "gene")
        denom: How to handle "UNMAPPED" features ("all" or "unmapped_excluded")
        filter_groups: Groups to include in analysis
        sample_id_col: Column in metadata for sample IDs
        alpha: Significance threshold
        
    Returns:
        Boolean indicating success
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Set default methods if not provided
    if methods is None:
        methods = ["aldex2", "ancom", "ancom-bc"]
    
    # Load abundance data
    try:
        abundance_df = pd.read_csv(abundance_file, sep='\t', index_col=0)
        logger.info(f"Loaded abundance data with {abundance_df.shape[0]} features and {abundance_df.shape[1]} samples")
    except Exception as e:
        logger.error(f"Error reading abundance file: {str(e)}")
        return False
    
    # Load metadata
    try:
        metadata_df = pd.read_csv(metadata_file)
        logger.info(f"Loaded metadata with {metadata_df.shape[0]} samples and {metadata_df.shape[1]} columns")
    except Exception as e:
        logger.error(f"Error reading metadata file: {str(e)}")
        return False
    
    # Auto-detect sample ID column if not specified
    if not sample_id_col:
        common_id_cols = ["SampleName", "Sample", "SampleID", "Sample_ID", "sample_name", "sample_id"]
        for col in common_id_cols:
            if col in metadata_df.columns:
                sample_id_col = col
                logger.info(f"Auto-detected sample ID column: {sample_id_col}")
                break
        
        if not sample_id_col:
            sample_id_col = metadata_df.columns[0]
            logger.warning(f"Could not auto-detect sample ID column, using the first column: {sample_id_col}")
    
    # Verify columns exist
    if sample_id_col not in metadata_df.columns:
        logger.error(f"Sample ID column '{sample_id_col}' not found in metadata")
        return False
    
    if group_col not in metadata_df.columns:
        logger.error(f"Group column '{group_col}' not found in metadata")
        return False
    
    # Set sample ID as index
    metadata_df = metadata_df.set_index(sample_id_col)
    
    # Process filter_groups if provided as string
    if isinstance(filter_groups, str):
        filter_groups = [g.strip() for g in filter_groups.split(',')]
    
    # Get unique groups in the full dataset
    unique_groups = metadata_df[group_col].unique()
    logger.info(f"Full dataset contains {len(unique_groups)} groups: {unique_groups}")
    
    # Check if filter_groups are valid
    if filter_groups:
        logger.info(f"Filtering to specified groups: {filter_groups}")
        missing_groups = [g for g in filter_groups if g not in unique_groups]
        if missing_groups:
            logger.error(f"The following specified groups don't exist in the data: {missing_groups}")
            logger.error(f"Available groups: {unique_groups}")
            return False
        
        # For ALDEx2, check if we have exactly 2 groups after filtering
        if "aldex2" in methods and len(filter_groups) != 2:
            logger.warning(f"ALDEx2 requires exactly 2 groups, but {len(filter_groups)} were specified")
            if "aldex2" in methods and len(methods) == 1:
                logger.error("Cannot proceed with ALDEx2 analysis without exactly 2 groups")
                return False
    
    # Results to track
    results = {}
    significant_features = {}
    
    # Run each method
    for method in methods:
        method_lower = method.lower()
        logger.info(f"Running {method_lower.upper()} analysis...")
        
        method_dir = os.path.join(output_dir, method_lower)
        os.makedirs(method_dir, exist_ok=True)
        
        if method_lower == "aldex2":
            if filter_groups and len(filter_groups) != 2:
                logger.warning(f"Skipping ALDEx2: found {len(filter_groups)} groups after filtering, but ALDEx2 requires exactly 2")
                continue
                
            try:
                result_df = aldex2_like(
                    abundance_df, 
                    metadata_df, 
                    group_col=group_col, 
                    denom=denom,
                    filter_groups=filter_groups
                )
                
                if result_df.empty:
                    logger.warning("ALDEx2 analysis did not produce results")
                    continue
                    
                # Save results
                output_file = os.path.join(method_dir, "aldex2_results.csv")
                result_df.to_csv(output_file)
                logger.info(f"Saved ALDEx2 results to {output_file}")
                
                # Track significant features
                sig_features = result_df[result_df['q_value'] < alpha]['feature'].tolist()
                significant_features["aldex2"] = set(sig_features)
                logger.info(f"ALDEx2: {len(sig_features)} significant features (q < {alpha})")
                
                # Create volcano plot
                plt.figure(figsize=(10, 6))
                plt.scatter(
                    result_df['effect_size'], 
                    -np.log10(result_df['p_value']),
                    alpha=0.7
                )
                # Highlight significant features
                sig_df = result_df[result_df['q_value'] < alpha]
                plt.scatter(
                    sig_df['effect_size'], 
                    -np.log10(sig_df['p_value']),
                    color='red',
                    alpha=0.7
                )
                plt.axhline(-np.log10(alpha), linestyle='--', color='gray')
                plt.axvline(0, linestyle='--', color='gray')
                plt.xlabel('Effect Size')
                plt.ylabel('-log10(p-value)')
                plt.title('ALDEx2 Volcano Plot')
                plt.savefig(os.path.join(method_dir, "aldex2_volcano.png"), dpi=300, bbox_inches='tight')
                plt.close()
                
                results["aldex2"] = result_df
                
            except Exception as e:
                logger.error(f"Error in ALDEx2 analysis: {str(e)}")
                
        elif method_lower == "ancom":
            try:
                result_df = ancom(
                    abundance_df, 
                    metadata_df, 
                    group_col=group_col, 
                    denom=denom,
                    filter_groups=filter_groups
                )
                
                if result_df.empty:
                    logger.warning("ANCOM analysis did not produce results")
                    continue
                    
                # Save results
                output_file = os.path.join(method_dir, "ancom_results.csv")
                result_df.to_csv(output_file)
                logger.info(f"Saved ANCOM results to {output_file}")
                
                # Track significant features
                sig_features = result_df[result_df['significant']]['feature'].tolist()
                significant_features["ancom"] = set(sig_features)
                logger.info(f"ANCOM: {len(sig_features)} significant features")
                
                # Bar plot for top ANCOM features
                top_ancom = result_df.head(20)
                plt.figure(figsize=(12, 8))
                plt.barh(top_ancom['feature'], top_ancom['W_ratio'])
                plt.axvline(0.7, linestyle='--', color='red', label='Significance threshold')
                plt.xlabel('W ratio')
                plt.ylabel('Feature')
                plt.title('Top 20 Features by ANCOM W-ratio')
                plt.tight_layout()
                plt.savefig(os.path.join(method_dir, "ancom_top_features.png"), dpi=300, bbox_inches='tight')
                plt.close()
                
                results["ancom"] = result_df
                
            except Exception as e:
                logger.error(f"Error in ANCOM analysis: {str(e)}")
                
        elif method_lower == "ancom-bc":
            try:
                result_df = ancom_bc(
                    abundance_df, 
                    metadata_df, 
                    group_col=group_col, 
                    denom=denom,
                    filter_groups=filter_groups
                )
                
                if result_df.empty:
                    logger.warning("ANCOM-BC analysis did not produce results")
                    continue
                    
                # Save results
                output_file = os.path.join(method_dir, "ancom_bc_results.csv")
                result_df.to_csv(output_file)
                logger.info(f"Saved ANCOM-BC results to {output_file}")
                
                # Track significant features
                sig_features = result_df[result_df['q_value'] < alpha]['feature'].tolist()
                significant_features["ancom-bc"] = set(sig_features)
                logger.info(f"ANCOM-BC: {len(sig_features)} significant features (q < {alpha})")
                
                # Create volcano plot
                plt.figure(figsize=(10, 6))
                plt.scatter(
                    result_df['effect_size'], 
                    -np.log10(result_df['p_value']),
                    alpha=0.7
                )
                # Highlight significant features
                sig_df = result_df[result_df['q_value'] < alpha]
                plt.scatter(
                    sig_df['effect_size'], 
                    -np.log10(sig_df['p_value']),
                    color='red',
                    alpha=0.7
                )
                plt.axhline(-np.log10(alpha), linestyle='--', color='gray')
                plt.axvline(0, linestyle='--', color='gray')
                plt.xlabel('Effect Size')
                plt.ylabel('-log10(p-value)')
                plt.title('ANCOM-BC Volcano Plot')
                plt.savefig(os.path.join(method_dir, "ancom_bc_volcano.png"), dpi=300, bbox_inches='tight')
                plt.close()
                
                results["ancom-bc"] = result_df
                
            except Exception as e:
                logger.error(f"Error in ANCOM-BC analysis: {str(e)}")
        
        else:
            logger.warning(f"Unknown method: {method}")
    
    # Compare results across methods if we have more than one
    if len(significant_features) > 1:
        logger.info("Comparing results across methods...")
        
        # Log the comparison information
        comparison_log = ["Overlap between significant features:"]
        for method, features in significant_features.items():
            comparison_log.append(f"{method.upper()} significant features: {len(features)}")
        
        # Pairwise comparisons
        methods = list(significant_features.keys())
        for i in range(len(methods)):
            for j in range(i+1, len(methods)):
                method1, method2 = methods[i], methods[j]
                overlap = len(significant_features[method1].intersection(significant_features[method2]))
                comparison_log.append(f"Overlap between {method1.upper()} and {method2.upper()}: {overlap}")
        
        # Three-way comparison if applicable
        if len(methods) >= 3:
            overlap = len(significant_features[methods[0]].intersection(
                significant_features[methods[1]]).intersection(
                significant_features[methods[2]]))
            comparison_log.append(f"Overlap between all three methods: {overlap}")
        
        # Log the comparison and save to file
        for line in comparison_log:
            logger.info(line)
            
        with open(os.path.join(output_dir, "method_comparison.txt"), "w") as f:
            f.write("\n".join(comparison_log))
        
        # Optionally, create a Venn diagram if matplotlib-venn is available
        try:
            from matplotlib_venn import venn2, venn3
            
            if len(methods) == 2:
                plt.figure(figsize=(8, 6))
                venn2([significant_features[methods[0]], significant_features[methods[1]]],
                    [methods[0].upper(), methods[1].upper()])
                plt.title("Overlap of Significant Features")
                plt.savefig(os.path.join(output_dir, "venn_diagram.png"), dpi=300, bbox_inches='tight')
                plt.close()
            elif len(methods) == 3:
                plt.figure(figsize=(8, 6))
                venn3([significant_features[methods[0]], significant_features[methods[1]], significant_features[methods[2]]],
                    [methods[0].upper(), methods[1].upper(), methods[2].upper()])
                plt.title("Overlap of Significant Features")
                plt.savefig(os.path.join(output_dir, "venn_diagram.png"), dpi=300, bbox_inches='tight')
                plt.close()
        except ImportError:
            logger.info("matplotlib-venn not available; skipping Venn diagram")
    
    return len(results) > 0#!/usr/bin/env python3
"""
HUMAnN3 Tools Differential Abundance Module

This module performs differential abundance analysis using methods that 
account for the compositional nature of microbiome data:
- ALDEx2-like analysis (effects + Welch's t-test)
- ANCOM (Analysis of Composition of Microbiomes)
- ANCOM-BC (Analysis of Composition of Microbiomes with Bias Correction)

Examples:
  # Basic usage:
  humann3-tools diff --abundance-file joined_output/pathway_abundance_cpm_unstratified.tsv \
                   --metadata-file metadata.csv \
                   --output-dir ./differential_abundance
  
  # Specify group column and methods:
  humann3-tools diff --abundance-file joined_output/gene_families_cpm_unstratified.tsv \
                   --metadata-file metadata.csv \
                   --group-col Treatment \
                   --methods aldex2,ancom \
                   --output-dir ./differential_abundance
                    
  # Filter groups for comparison (required for ALDEx2):
  humann3-tools diff --abundance-file joined_output/pathway_abundance_cpm_unstratified.tsv \
                   --metadata-file metadata.csv \
                   --methods aldex2 \
                   --filter-groups Control,Treatment \
                   --output-dir ./differential_abundance
                    
  # Exclude unmapped features:
  humann3-tools diff --abundance-file joined_output/pathway_abundance_cpm_unstratified.tsv \
                   --metadata-file metadata.csv \
                   --exclude-unmapped \
                   --output-dir ./differential_abundance
"""

import os
import sys
import time
import argparse
import logging
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple, Union
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests

# Import internal modules
try:
    from humann3_tools.utils.resource_utils import track_peak_memory
except ImportError:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
    from humann3_tools.utils.resource_utils import track_peak_memory

# Set up logging
logger = logging.getLogger('humann3_tools')

def setup_logger(log_file=None, log_level=logging.INFO):
    """Set up the logger with console and optional file output."""
    # Remove any existing handlers to avoid duplication
    logger.handlers = []
    logger.setLevel(log_level)

    # Format for logs
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # Optional file handler
    if log_file:
        os.makedirs(os.path.dirname(os.path.abspath(log_file)), exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger

# Custom CLR implementation to avoid skbio dependency
def clr_transform(data_matrix):
    """
    Compute the centered log-ratio (CLR) transformation.
    
    Parameters:
    -----------
    data_matrix : numpy.ndarray or pandas.DataFrame
        The data matrix to transform, with features as columns
        
    Returns:
    --------
    numpy.ndarray or pandas.DataFrame
        The CLR-transformed data matrix
    """
    if isinstance(data_matrix, pd.DataFrame):
        # For pandas DataFrame
        log_data = np.log(data_matrix)
        geometric_means = log_data.mean(axis=1)
        clr_data = log_data.subtract(geometric_means, axis=0)
        return clr_data
    else:
        # For numpy array
        log_data = np.log(data_matrix)
        geometric_means = np.mean(log_data, axis=1, keepdims=True)
        return log_data - geometric_means

def aldex2_like(
    abundance_df: pd.DataFrame, 
    metadata_df: pd.DataFrame, 
    group_col: str, 
    mc_samples: int = 128, 
    denom: str = "all", 
    filter_groups: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    A Python implementation similar to ALDEx2 for differential abundance testing
    
    Parameters:
        abundance_df: Feature table with samples as columns and features as rows
        metadata_df: Metadata with sample IDs as index and metadata as columns
        group_col: Column name in metadata_df that contains the grouping variable
        mc_samples: Number of Monte Carlo samples to generate
        denom: Features to use as denominator: "all" for all features, 
               "unmapped_excluded" to exclude unmapped features
        filter_groups: List of group names to include in the analysis. If provided, 
                      only these groups will be used. Must contain exactly 2 groups for ALDEx2.
        
    Returns:
        pandas DataFrame with test results
    """