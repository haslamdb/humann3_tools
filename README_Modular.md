# HUMAnN3 Tools - Modular CLI Usage Guide

This guide explains how to use the redesigned modular command-line tools for processing and analyzing HUMAnN3 output.

## Overview 

The HUMAnN3 tools package has been redesigned with a modular architecture that allows you to:

1. Run each step independently
2. Process data without needing to rerun previous steps
3. Mix and match analysis methods based on your needs

The redesigned tools include:

- `humann3-join`: Join, normalize, and unstratify HUMAnN3 output files
- `humann3-diff`: Run differential abundance analysis on processed files
- `humann3-stats`: Run statistical tests like Kruskal-Wallis and Dunn's post-hoc tests
- `humann3-viz`: Generate visualizations from processed data

## Installation

```bash
# Clone the repository
git clone https://github.com/your-repo/humann3_tools.git

# Install the package
cd humann3_tools
pip install -e .
```

## Usage Examples

### 1. Join and Normalize HUMAnN3 Output Files

This step can be run independently after HUMAnN3 has generated output files:

```bash
# Process pathway abundance files
humann3-join --input-dir /path/to/humann3_output \
  --pathabundance \
  --output-dir ./processed_data \
  --units cpm \
  --update-snames

# Process gene family files
humann3-join --input-dir /path/to/humann3_output \
  --genefamilies \
  --output-dir ./processed_data \
  --units cpm \
  --update-snames
```

### 2. Run Differential Abundance Analysis

After joining and unstratifying your files, you can run differential abundance analysis:

```bash
# For pathway abundance data
humann3-diff --abundance-file ./processed_data/pathabundance_cpm_unstratified.tsv \
  --metadata-file ./sample_metadata.csv \
  --output-dir ./differential_abundance \
  --group-col Group \
  --methods aldex2,ancom,ancom-bc,kruskal

# For gene family data
humann3-diff --abundance-file ./processed_data/genefamilies_cpm_unstratified.tsv \
  --metadata-file ./sample_metadata.csv \
  --output-dir ./differential_abundance \
  --feature-type gene \
  --group-col Group \
  --methods aldex2,ancom,ancom-bc
```

### 3. Run Statistical Tests

If you're specifically interested in Kruskal-Wallis and Dunn's post-hoc tests:

```bash
humann3-stats --abundance-file ./processed_data/pathabundance_cpm_unstratified.tsv \
  --metadata-file ./sample_metadata.csv \
  --output-dir ./statistical_tests \
  --group-col Group \
  --alpha 0.05
```

### 4. Generate Visualizations

Create visualizations from your processed data:

```bash
humann3-viz --abundance-file ./processed_data/pathabundance_cpm_unstratified.tsv \
  --metadata-file ./sample_metadata.csv \
  --output-dir ./visualizations \
  --group-col Group \
  --shape-col Treatment \
  --pca --heatmap --barplot --abundance-hist \
  --format svg \
  --top-n 25
```

## Command Options

### humann3-join

```
usage: humann3-join [-h] --input-dir INPUT_DIR
                    (--pathabundance | --pathcoverage | --genefamilies)
                    [--output-dir OUTPUT_DIR] [--output-basename OUTPUT_BASENAME]
                    [--units {cpm,relab}] [--update-snames]
                    [--file-pattern FILE_PATTERN] [--log-file LOG_FILE]
                    [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

Join, normalize, and unstratify HUMAnN3 output files

required arguments:
  --input-dir INPUT_DIR    Directory containing HUMAnN3 output files
  --pathabundance          Process pathabundance files
  --pathcoverage           Process pathcoverage files
  --genefamilies           Process genefamilies files

optional arguments:
  --output-dir OUTPUT_DIR  Directory for output files (default: ./joined_output)
  --output-basename OUTPUT_BASENAME
                           Base filename for output (default derived from file type)
  --units {cpm,relab}      Units for normalization (default: cpm)
  --update-snames          Update sample names during normalization
  --file-pattern FILE_PATTERN
                           Glob pattern for input files (default based on file type)
  --log-file LOG_FILE      Path to log file
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                           Logging level (default: INFO)
```

### humann3-diff

```
usage: humann3-diff [-h] --abundance-file ABUNDANCE_FILE --metadata-file METADATA_FILE
                   [--output-dir OUTPUT_DIR] [--feature-type {pathway,gene}]
                   [--group-col GROUP_COL] [--sample-id-col SAMPLE_ID_COL]
                   [--methods METHODS] [--exclude-unmapped]
                   [--log-file LOG_FILE]
                   [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

Run differential abundance analysis on HUMAnN3 output files

required arguments:
  --abundance-file ABUNDANCE_FILE
                           Path to the unstratified abundance file (pathway or gene family)
  --metadata-file METADATA_FILE
                           Path to sample metadata CSV file

optional arguments:
  --output-dir OUTPUT_DIR  Directory for output files (default: ./DifferentialAbundance)
  --feature-type {pathway,gene}
                           Type of features in the abundance file (default: pathway)
  --group-col GROUP_COL    Column name in metadata for grouping samples (default: Group)
  --sample-id-col SAMPLE_ID_COL
                           Column name in metadata for sample IDs (autodetected if not specified)
  --methods METHODS        Comma-separated list of methods to use (default: aldex2,ancom,ancom-bc)
  --exclude-unmapped       Exclude unmapped features from analysis
  --log-file LOG_FILE      Path to log file
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                           Logging level (default: INFO)
```

### humann3-stats

```
usage: humann3-stats [-h] --abundance-file ABUNDANCE_FILE --metadata-file METADATA_FILE
                    [--output-dir OUTPUT_DIR] [--feature-type {pathway,gene}]
                    [--group-col GROUP_COL] [--sample-id-col SAMPLE_ID_COL]
                    [--alpha ALPHA] [--log-file LOG_FILE]
                    [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

Run statistical tests on HUMAnN3 output files

required arguments:
  --abundance-file ABUNDANCE_FILE
                           Path to the unstratified abundance file (pathway or gene family)
  --metadata-file METADATA_FILE
                           Path to sample metadata CSV file

optional arguments:
  --output-dir OUTPUT_DIR  Directory for output files (default: ./StatisticalTests)
  --feature-type {pathway,gene}
                           Type of features in the abundance file (default: pathway)
  --group-col GROUP_COL    Column name in metadata for grouping samples (default: Group)
  --sample-id-col SAMPLE_ID_COL
                           Column name in metadata for sample IDs (autodetected if not specified)
  --alpha ALPHA            Significance threshold for statistical tests (default: 0.05)
  --log-file LOG_FILE      Path to log file
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                           Logging level (default: INFO)
```

### humann3-viz

```
usage: humann3-viz [-h] --abundance-file ABUNDANCE_FILE --metadata-file METADATA_FILE
                  [--output-dir OUTPUT_DIR] [--feature-type {pathway,gene}]
                  [--group-col GROUP_COL] [--shape-col SHAPE_COL]
                  [--sample-id-col SAMPLE_ID_COL] [--top-n TOP_N]
                  [--format {svg,png,pdf}] [--dpi DPI] [--pca] [--heatmap]
                  [--barplot] [--abundance-hist] [--log-transform]
                  [--log-file LOG_FILE]
                  [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

Create visualizations for HUMAnN3 output files

required arguments:
  --abundance-file ABUNDANCE_FILE
                          Path to the unstratified abundance file (pathway or gene family)
  --metadata-file METADATA_FILE
                          Path to sample metadata CSV file

optional arguments:
  --output-dir OUTPUT_DIR  Directory for output files (default: ./Visualizations)
  --feature-type {pathway,gene}
                          Type of features in the abundance file (default: pathway)
  --group-col GROUP_COL    Column name in metadata for coloring points (default: Group)
  --shape-col SHAPE_COL    Column name in metadata for point shapes
  --sample-id-col SAMPLE_ID_COL
                          Column name in metadata for sample IDs (autodetected if not specified)
  --top-n TOP_N            Number of top features to include in bar plots (default: 25)
  --format {svg,png,pdf}   Output format for plots (default: svg)
  --dpi DPI                DPI for raster formats like PNG (default: 300)
  --pca                    Generate PCA plot (default: True)
  --heatmap                Generate heatmap of top features
  --barplot                Generate barplot of top features by group (default: True)
  --abundance-hist         Generate histograms of abundance distributions
  --log-transform          Apply log10(x+1) transformation to abundance data (default: True)
  --log-file LOG_FILE      Path to log file
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                          Logging level (default: INFO)
```

## Workflow Examples

### Complete Analysis Workflow

Here's a complete workflow from joining to visualization:

```bash
# 1. Join and normalize pathway abundance files
humann3-join --input-dir /path/to/humann3_output \
  --pathabundance \
  --output-dir ./processed_data \
  --units cpm

# 2. Run differential abundance analysis
humann3-diff --abundance-file ./processed_data/pathabundance_cpm_unstratified.tsv \
  --metadata-file ./sample_metadata.csv \
  --output-dir ./results/differential_abundance \
  --methods aldex2,ancom,ancom-bc,kruskal

# 3. Generate visualizations
humann3-viz --abundance-file ./processed_data/pathabundance_cpm_unstratified.tsv \
  --metadata-file ./sample_metadata.csv \
  --output-dir ./results/visualizations \
  --pca --heatmap --barplot
```

### Working with Both Pathway and Gene Family Data

Process both data types in parallel:

```bash
# Join and normalize both data types
humann3-join --input-dir /path/to/humann3_output --pathabundance --output-dir ./processed_data
humann3-join --input-dir /path/to/humann3_output --genefamilies --output-dir ./processed_data

# Run differential abundance analysis on both
humann3-diff --abundance-file ./processed_data/pathabundance_cpm_unstratified.tsv \
  --metadata-file ./sample_metadata.csv --output-dir ./results/differential_abundance/pathways

humann3-diff --abundance-file ./processed_data/genefamilies_cpm_unstratified.tsv \
  --metadata-file ./sample_metadata.csv --output-dir ./results/differential_abundance/genes \
  --feature-type gene

# Create visualizations for both
humann3-viz --abundance-file ./processed_data/pathabundance_cpm_unstratified.tsv \
  --metadata-file ./sample_metadata.csv --output-dir ./results/visualizations/pathways

humann3-viz --abundance-file ./processed_data/genefamilies_cpm_unstratified.tsv \
  --metadata-file ./sample_metadata.csv --output-dir ./results/visualizations/genes \
  --feature-type gene
```

## Tips and Best Practices

1. **File Organization**: Keep your data organized with clear directory structures. Separate raw data, processed files, and results.

2. **Metadata Preparation**: Ensure your metadata CSV file has clear sample IDs that match those in your HUMAnN3 output files. Include all relevant grouping variables.

3. **Log Files**: Always use the `--log-file` option to save processing logs, which will help with troubleshooting.

4. **Visualization Format**: Use SVG format for publication-quality figures (scalable), and PNG for presentations or web display.

5. **Parallel Processing**: To process multiple datasets, you can run these tools in parallel batch jobs on a cluster.

6. **Memory Usage**: For large datasets, be aware of memory usage. The differential abundance analysis, especially ANCOM, can be memory-intensive.