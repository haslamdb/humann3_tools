# Installing the Updated HUMAnN3 Tools Package

This guide explains how to install the updated HUMAnN3 Tools package with both the original pipeline and the new modular command-line tools.

## Installation

1. Update your `pyproject.toml` file with the changes provided.

2. Reinstall the package with pip:

```bash
# From the package directory
pip install -e .
```

This will install the package in "editable" mode, allowing you to make further changes without reinstalling.

## Verifying Installation

After installation, you should have access to these command-line tools:

```bash
# Original tools
humann3-tools --help
join_unstratify_humann_output --help

# New modular tools
humann3-join --help
humann3-diff --help
humann3-stats --help
humann3-viz --help
```

## Using the New Modular Tools

### 1. Join and Normalize HUMAnN3 Output Files

```bash
humann3-join --input-dir /path/to/humann3_output \
  --pathabundance \
  --output-dir ./processed_data \
  --units cpm
```

### 2. Run Differential Abundance Analysis

```bash
humann3-diff --abundance-file ./processed_data/pathabundance_cpm_unstratified.tsv \
  --metadata-file ./sample_metadata.csv \
  --output-dir ./differential_abundance \
  --methods aldex2,ancom,ancom-bc,kruskal
```

### 3. Run Statistical Tests

```bash
humann3-stats --abundance-file ./processed_data/pathabundance_cpm_unstratified.tsv \
  --metadata-file ./sample_metadata.csv \
  --output-dir ./statistical_tests \
  --group-col Group
```

### 4. Generate Visualizations

```bash
humann3-viz --abundance-file ./processed_data/pathabundance_cpm_unstratified.tsv \
  --metadata-file ./sample_metadata.csv \
  --output-dir ./visualizations \
  --group-col Group \
  --pca --heatmap --barplot
```

## Using Original Pipeline

You can still use the original pipeline for full workflow execution:

```bash
humann3-tools --run-preprocessing \
  --input-fastq /path/to/fastq_files/*.fastq.gz \
  --sample-key metadata.csv \
  --output-dir output_directory
```

Or just the join and unstratify functionality:

```bash
join_unstratify_humann_output \
  --sample-key metadata.csv \
  --pathway-dir /path/to/pathway_files \
  --gene-dir /path/to/gene_files \
  --output-dir output_directory
```

## Mixing Approaches

You can mix both approaches in your workflow:

```bash
# Run the full pipeline for preprocessing and basic analysis
humann3-tools --run-preprocessing --input-fastq *.fastq.gz --sample-key metadata.csv

# Then use the modular tools for more specialized analysis
humann3-diff --abundance-file output/pathways/pathway_abundance_cpm_unstratified.tsv \
  --metadata-file metadata.csv \
  --output-dir differential_analysis \
  --methods aldex2,ancom,kruskal

# Generate custom visualizations
humann3-viz --abundance-file output/pathways/pathway_abundance_cpm_unstratified.tsv \
  --metadata-file metadata.csv \
  --output-dir visualizations \
  --heatmap --top-n 50
```

This gives you the flexibility to use both pipeline approaches as needed for your analysis.
