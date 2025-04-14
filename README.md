# HUMAnN3 Tools

A comprehensive Python package for assigning raw metagenomic sequence reads to microbial gene and pathway databases using HUMAnN3, followed by downstream processing and analysis.

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
  - [Step 1: Set up bioBakery3 environment](#step-1-set-up-biobakery3-environment)
  - [Step 2: Install HUMAnN3 Tools](#step-2-install-humann3-tools)
  - [Step 3: Verify Installation](#step-3-verify-installation)
- [Workflow Overview](#workflow-overview)
- [Command Line Interface](#command-line-interface)
  - [1. KneadData](#1-kneaddata)
  - [2. HUMAnN3](#2-humann3)
  - [3. Join and Normalize](#3-join-and-normalize)
  - [4. Statistical Testing](#4-statistical-testing)
  - [5. Differential Abundance](#5-differential-abundance)
  - [6. Visualization](#6-visualization)
- [Input Methods](#input-methods)
- [Example Workflow](#example-workflow)
- [Troubleshooting](#troubleshooting)
- [Getting Help](#getting-help)
- [License](#license)

## Introduction

HUMAnN3 Tools provides a complete workflow for metagenomic analysis, from quality control to visualization. The package is designed to simplify and standardize the analysis of metagenomic data using the HUMAnN3 pipeline, with modular components that can be run individually or as a complete workflow.

## Installation

### Step 1: Set up bioBakery3 environment

**IMPORTANT**: HUMAnN3 Tools requires the bioBakery3 conda environment, which includes properly configured HUMAnN3, KneadData, and other dependencies. You must activate this environment before using HUMAnN3 Tools.

```bash
# Install bioBakery3 environment (if not already installed)
conda create -n biobakery3 python=3.12
conda activate biobakery3
conda install -c biobakery -c conda-forge -c bioconda humann=3.9 kneaddata=0.10.0

# Download and install necessary databases
# MetaPhlAn database (required for HUMAnN3)
metaphlan --install

# Install HUMAnN databases (if not already installed)
# humann_databases --download chocophlan full /path/to/databases
# humann_databases --download uniref uniref90_diamond /path/to/databases
# humann_databases --download utility_mapping full /path/to/databases
```

### Step 2: Install HUMAnN3 Tools

**ALWAYS ENSURE biobakery3 ENVIRONMENT IS ACTIVATED BEFORE RUNNING ANY COMMANDS**

```bash
# Make sure biobakery3 is activated
conda activate biobakery3

# Install HUMAnN3 Tools from the repository
pip install git+https://github.com/haslamdb/humann3_tools.git
```

For development installation:

```bash
# Clone the repository
git clone https://github.com/haslamdb/humann3_tools.git
cd humann3_tools

# Make sure biobakery3 is activated
conda activate biobakery3

# Install in development mode
pip install -e .
```

After installation, you can use either:
- Individual commands: `humann3-kneaddata`, `humann3-humann3`, etc.
- The unified interface: `humann3-tools kneaddata`, `humann3-tools humann3`, etc.

### Step 3: Verify Installation

To verify that HUMAnN3 Tools is installed correctly:

```bash
# Make sure biobakery3 is activated
conda activate biobakery3

# Check version
humann3-tools --version

# List available commands
humann3-tools --help
```

## Workflow Overview

HUMAnN3 Tools provides a modular workflow for metagenomic analysis:

1. **KneadData**: Quality control and host depletion
2. **HUMAnN3**: Process cleaned sequences through HUMAnN3
3. **Join & Normalize**: Combine, normalize, and split HUMAnN3 output files
4. **Statistical Testing**: Perform statistical tests across groups
5. **Differential Abundance**: Apply methods like ALDEx2, ANCOM, and ANCOM-BC
6. **Visualization**: Generate plots and figures for publication

![Workflow Diagram](workflow_diagram.png)

## Command Line Interface

HUMAnN3 Tools provides a consistent command-line interface for each step of the workflow. Each command can be accessed either through the main `humann3-tools` command or as individual commands.

**IMPORTANT: Always make sure the biobakery3 environment is activated before running these commands.**

```bash
conda activate biobakery3
```

### 1. KneadData

Quality control and host depletion:

```bash
# Using the main command
humann3-tools kneaddata --input-files sample_R1.fastq.gz sample_R2.fastq.gz --paired --reference-dbs human_db --output-dir kneaddata_output

# Using the standalone command
humann3-kneaddata --input-files sample_R1.fastq.gz sample_R2.fastq.gz --paired --reference-dbs human_db --output-dir kneaddata_output
```

Key options:
- `--input-files`: Input FASTQ file(s)
- `--paired`: Flag if input files are paired-end reads
- `--reference-dbs`: Path to reference database(s) for decontamination
- `--output-dir`: Directory for output files
- `--threads`: Number of threads to use

### 2. HUMAnN3

Run HUMAnN3 on cleaned sequence files:

```bash
# Using the main command
humann3-tools humann3 --input-dir kneaddata_output --output-dir humann3_output --nucleotide-db chocophlan --protein-db uniref

# Using the standalone command
humann3-humann3 --input-dir kneaddata_output --output-dir humann3_output --nucleotide-db chocophlan --protein-db uniref
```

Key options:
- `--input-dir`: Directory containing KneadData output files
- `--input-files`: Alternatively, specify cleaned files directly
- `--nucleotide-db`: Path to nucleotide database
- `--protein-db`: Path to protein database
- `--output-dir`: Directory for output files
- `--threads`: Number of threads to use
- `--use-parallel`: Process multiple samples in parallel
- `--bypass-prescreen`: Skip MetaPhlAn taxonomic prescreen (useful if MetaPhlAn database isn't installed)

### 3. Join and Normalize

Join, normalize, and unstratify HUMAnN3 output files:

```bash
# Using the main command
humann3-tools join --input-dir humann3_output/PathwayAbundance --pathabundance --output-dir joined_output --units cpm

# Using the standalone command
humann3-join --input-dir humann3_output/PathwayAbundance --pathabundance --output-dir joined_output --units cpm
```

Key options:
- `--input-dir`: Directory with HUMAnN3 output files
- One of: `--pathabundance`, `--pathcoverage`, or `--genefamilies` to specify file type
- `--output-dir`: Directory for output files
- `--units`: Units for normalization (cpm or relab)
- `--update-snames`: Update sample names during normalization

### 4. Statistical Testing

Run statistical tests on processed data:

```bash
# Using the main command
humann3-tools stats --abundance-file joined_output/pathabundance_cpm_unstratified.tsv --metadata-file metadata.csv --output-dir statistical_results --group-col Treatment

# Using the standalone command
humann3-stats --abundance-file joined_output/pathabundance_cpm_unstratified.tsv --metadata-file metadata.csv --output-dir statistical_results --group-col Treatment
```

Key options:
- `--abundance-file`: Path to unstratified abundance file
- `--metadata-file`: Path to metadata CSV file
- `--output-dir`: Directory for output files
- `--group-col`: Column name for grouping samples
- `--feature-type`: Type of features (pathway or gene)
- `--alpha`: Significance threshold (default: 0.05)

### 5. Differential Abundance

Run differential abundance analysis:

```bash
# Using the main command
humann3-tools diff --abundance-file joined_output/pathabundance_cpm_unstratified.tsv --metadata-file metadata.csv --output-dir diff_abundance --group-col Treatment --methods aldex2,ancom,ancom-bc

# Using the standalone command
humann3-diff --abundance-file joined_output/pathabundance_cpm_unstratified.tsv --metadata-file metadata.csv --output-dir diff_abundance --group-col Treatment --methods aldex2,ancom,ancom-bc
```

Key options:
- `--abundance-file`: Path to unstratified abundance file
- `--metadata-file`: Path to metadata CSV file
- `--output-dir`: Directory for output files
- `--group-col`: Column name for grouping samples
- `--methods`: Methods to use (aldex2, ancom, ancom-bc)
- `--filter-groups`: Filter groups for comparison (required for ALDEx2)
- `--exclude-unmapped`: Exclude unmapped features from analysis

### 6. Visualization

Create visualizations from processed data:

```bash
# Using the main command
humann3-tools viz --abundance-file joined_output/pathabundance_cpm_unstratified.tsv --metadata-file metadata.csv --output-dir visualizations --pca --heatmap --barplot

# Using the standalone command
humann3-viz --abundance-file joined_output/pathabundance_cpm_unstratified.tsv --metadata-file metadata.csv --output-dir visualizations --pca --heatmap --barplot
```

Key options:
- `--abundance-file`: Path to unstratified abundance file
- `--metadata-file`: Path to metadata CSV file
- `--output-dir`: Directory for output files
- `--group-col`: Column name for coloring points
- Plot selection: `--pca`, `--heatmap`, `--barplot`, `--abundance-hist`
- `--feature`: Generate boxplot for a specific feature
- `--format`: Output format (svg, png, pdf)

## Input Methods

HUMAnN3 Tools supports three different input methods across all commands:

1. **Direct File Input**: Specify the input files directly with `--input-files`
2. **Sample List File**: Provide a tab-delimited file with sample IDs and file paths using `--samples-file`
3. **Metadata-Driven**: Use a metadata CSV file to locate sequence files with `--metadata-file` and `--seq-dir`

Example for metadata-driven workflow:

```bash
humann3-tools kneaddata --metadata-file metadata.csv --seq-dir /path/to/sequences --r1-suffix _R1.fastq.gz --r2-suffix _R2.fastq.gz --paired --reference-dbs human_db
```

## Example Workflow

Here's a complete example workflow:

```bash
# Activate bioBakery3 environment (REQUIRED)
conda activate biobakery3

# 1. Quality control with KneadData
humann3-tools kneaddata --input-files sample1_R1.fastq.gz sample1_R2.fastq.gz --paired --reference-dbs human_db --output-dir kneaddata_output --threads 8

# 2. HUMAnN3 functional profiling
humann3-tools humann3 --input-dir kneaddata_output --output-dir humann3_output --nucleotide-db chocophlan --protein-db uniref --threads 8

# 3. Join and normalize pathway abundance files
humann3-tools join --input-dir humann3_output/PathwayAbundance --pathabundance --output-dir joined_output --units cpm

# 4. Join and normalize gene family files
humann3-tools join --input-dir humann3_output/GeneFamilies --genefamilies --output-dir joined_output --units cpm

# 5. Statistical testing
humann3-tools stats --abundance-file joined_output/pathabundance_cpm_unstratified.tsv --metadata-file metadata.csv --output-dir statistical_results --group-col Treatment

# 6. Differential abundance analysis
humann3-tools diff --abundance-file joined_output/pathabundance_cpm_unstratified.tsv --metadata-file metadata.csv --output-dir diff_abundance --group-col Treatment

# 7. Visualization
humann3-tools viz --abundance-file joined_output/pathabundance_cpm_unstratified.tsv --metadata-file metadata.csv --output-dir visualizations --pca --heatmap --barplot
```

## Troubleshooting

### Common Issues

1. **Environment Issues**:
   - Make sure you've activated the biobakery3 environment: `conda activate biobakery3`
   - If you see errors about missing databases or tools, this is often because the biobakery3 environment is not activated

2. **Missing Databases**:
   - MetaPhlAn error: Run `metaphlan --install` to install the required database
   - HUMAnN3 database errors: Make sure to install the ChocoPhlAn and UniRef databases

3. **Missing Input Files**: 
   - Check file paths and naming patterns
   - Verify that you're using the correct directory structure

4. **Sample Key Issues**:
   - Ensure sample identifiers in the CSV match the file names
   - Check for duplicate sample IDs
   - Verify CSV encoding (use UTF-8)

5. **KneadData Errors**:
   - Ensure reference databases are properly built and indexed
   - For paired-end issues, try different `--decontaminate-pairs` options

6. **HUMAnN3 Errors**:
   - If you get a "No MetaPhlAn BowTie2 database found" error, run `metaphlan --install`
   - As a workaround, you can use `--bypass-prescreen` to skip the MetaPhlAn step
   - Ensure nucleotide and protein databases are correctly installed
   - Check that input files are in the correct format

7. **Statistical Test Errors**:
   - ALDEx2 requires exactly two groups
   - Check for missing values in abundance data

8. **Memory Issues with Large Datasets**:
   - Run steps separately to manage memory usage
   - Use `--threads` to control CPU usage

## Getting Help

If you encounter issues not covered in this documentation, please:

- Check the log file for detailed error messages
- Set `--log-level DEBUG` for more verbose output
- Verify you're using the biobakery3 environment
- Open an issue on the GitHub repository with a description of the problem and relevant log entries

## License

This project is licensed under the MIT License - see the LICENSE file for details.
