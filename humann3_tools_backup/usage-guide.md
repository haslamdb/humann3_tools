# HUMAnN3 Tools - Modular CLI Usage Guide

This guide explains how to use the redesigned modular command-line tools for processing and analyzing HUMAnN3 output, including metadata-driven workflows.

## Overview 

The HUMAnN3 tools package has been redesigned with a modular architecture that allows you to:

1. Run each step independently
2. Process data without needing to rerun previous steps
3. Mix and match analysis methods based on your needs
4. Use metadata files to drive your workflows

The redesigned tools include:

- `humann3-tools`: Main workflow command for running the complete pipeline
- `humann3-kneaddata`: Run KneadData for quality control and host depletion
- `humann3-run`: Run HUMAnN3 on processed sequence files
- `humann3-preprocess`: Run combined KneadData and HUMAnN3 steps
- `humann3-join`: Join, normalize, and unstratify HUMAnN3 output files
- `humann3-diff`: Run differential abundance analysis on processed files
- `humann3-stats`: Run statistical tests like Kruskal-Wallis and Dunn's post-hoc tests
- `humann3-viz`: Generate visualizations from processed data
- `join_unstratify_humann_output`: Combined join and unstratify operation

## Installation

```bash
# Clone the repository
git clone https://github.com/dhaslam/humann3_tools.git

# Install the package
cd humann3_tools
pip install -e .
```

## Input File Approaches

HUMAnN3 Tools supports multiple ways to specify input files:

1. **Direct file specification**: Provide individual file paths
2. **Metadata-driven workflow**: Use a metadata CSV file with sample information
3. **Samples file**: Use a tab-delimited file that maps sample IDs to file paths
4. **Directory with pattern matching**: Process files in a directory that match a specific pattern

## Metadata-Driven Workflows

### Using a Metadata CSV File

A metadata CSV file typically contains:
- A column for sample IDs
- Optional columns for file paths or suffixes
- Grouping variables for statistical analysis

Example metadata.csv:
```
SampleID,Group,Treatment,R1_file,R2_file
sample1,Control,Drug1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,Treatment,Drug1,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
sample3,Control,Drug2,/path/to/sample3_R1.fastq.gz,/path/to/sample3_R2.fastq.gz
```

### Using a Samples File

A tab-delimited samples file maps sample IDs to their sequence files:

Example samples.txt:
```
sample1    /path/to/sample1_R1.fastq.gz /path/to/sample1_R2.fastq.gz
sample2    /path/to/sample2_R1.fastq.gz /path/to/sample2_R2.fastq.gz
sample3    /path/to/sample3_R1.fastq.gz /path/to/sample3_R2.fastq.gz
```

## Usage Examples

### 1. Complete Metadata-Driven Workflow

Process all samples from a metadata file in a single command:

```bash
# Run complete workflow using metadata
humann3-tools --run-preprocessing --use-metadata \
    --sample-key metadata.csv \
    --seq-dir /path/to/sequence/files \
    --paired \
    --r1-suffix "_R1.fastq.gz" --r2-suffix "_R2.fastq.gz" \
    --kneaddata-dbs /path/to/kneaddata_db \
    --humann3-nucleotide-db /path/to/chocophlan \
    --humann3-protein-db /path/to/uniref \
    --output-dir ./complete_results \
    --group-col "Treatment" \
    --run-diff-abundance \
    --threads 8
```

### 2. Preprocessing with Metadata

```bash
# Run only preprocessing steps using metadata
humann3-preprocess --use-metadata \
    --metadata-file metadata.csv \
    --seq-dir /path/to/sequence/files \
    --sample-col "SampleID" \
    --r1-suffix "_R1.fastq.gz" --r2-suffix "_R2.fastq.gz" \
    --paired \
    --kneaddata-dbs /path/to/kneaddata_db \
    --humann3-nucleotide-db /path/to/chocophlan \
    --humann3-protein-db /path/to/uniref \
    --output-dir ./preprocessing_results \
    --threads 8
```

### 3. KneadData with Metadata

```bash
# Run only KneadData using metadata
humann3-kneaddata --use-metadata \
    --metadata-file metadata.csv \
    --seq-dir /path/to/sequence/files \
    --sample-col "SampleID" \
    --r1-suffix "_R1.fastq.gz" --r2-suffix "_R2.fastq.gz" \
    --paired \
    --reference-dbs /path/to/kneaddata_db \
    --output-dir ./kneaddata_output \
    --threads 8
```

### 4. HUMAnN3 Run with Sample List

```bash
# Run HUMAnN3 on KneadData outputs
humann3-run --input-files \
    kneaddata_output/sample1_paired_concat.fastq \
    kneaddata_output/sample2_paired_concat.fastq \
    kneaddata_output/sample3_paired_concat.fastq \
    --nucleotide-db /path/to/chocophlan \
    --protein-db /path/to/uniref \
    --output-dir ./humann3_output \
    --threads 8
```

### 5. Using a Samples File

```bash
# Using a tab-delimited samples file
humann3-tools --run-preprocessing \
    --samples-file samples.txt \
    --paired \
    --kneaddata-dbs /path/to/kneaddata_db \
    --humann3-nucleotide-db /path/to/chocophlan \
    --humann3-protein-db /path/to/uniref \
    --output-dir ./results \
    --sample-key metadata.csv \
    --group-col "Treatment" \
    --threads 8
```

### 6. Join and Normalize HUMAnN3 Output Files

After preprocessing, join the output files:

```bash
# Process pathway abundance files from a directory
humann3-join --input-dir ./humann3_output \
    --pathabundance \
    --output-dir ./processed_data \
    --units cpm

# Alternative using join_unstratify_humann_output with sample_key
join_unstratify_humann_output \
    --sample-key metadata.csv \
    --pathway-dir ./humann3_output/PathwayAbundance \
    --gene-dir ./humann3_output/GeneFamilies \
    --output-dir ./processed_data \
    --units cpm
```

### 7. Run Differential Abundance Analysis

Using metadata file for grouping:

```bash
# For pathway abundance data
humann3-diff --abundance-file ./processed_data/pathway_abundance-cpm_unstratified.tsv \
    --metadata-file metadata.csv \
    --output-dir ./differential_abundance \
    --group-col Treatment \
    --methods aldex2,ancom,ancom-bc
```

#### Comparing specific groups

You can use the `--filter-groups` option to perform comparisons between specific groups:

```bash
# ALDEx2 requires exactly 2 groups
humann3-diff --abundance-file ./processed_data/pathway_abundance-cpm_unstratified.tsv \
    --metadata-file metadata.csv \
    --output-dir ./differential_abundance/control_vs_treatment1 \
    --group-col Treatment \
    --methods aldex2 \
    --filter-groups Control,Treatment1

# ANCOM and ANCOM-BC can compare any number of groups
humann3-diff --abundance-file ./processed_data/pathway_abundance-cpm_unstratified.tsv \
    --metadata-file metadata.csv \
    --output-dir ./differential_abundance/treatment_comparison \
    --group-col Treatment \
    --methods ancom,ancom-bc \
    --filter-groups Treatment1,Treatment2,Treatment3
```

### humann3-diff

```
usage: humann3-diff --abundance-file ABUNDANCE_FILE 
                   --metadata-file METADATA_FILE
                   [--group-col GROUP_COL]
                   [--sample-id-col SAMPLE_ID_COL]
                   [--methods METHODS]
                   [--filter-groups GROUPS]
```

Key parameters for humann3-diff:
```
Required:
  --abundance-file        Path to unstratified abundance file
  --metadata-file         Path to metadata CSV file

Options:
  --group-col             Column in metadata for grouping samples (default: "Group")
  --methods               Comma-separated list of methods (default: "aldex2,ancom,ancom-bc")
  --filter-groups         Comma-separated list of group names to include in the analysis.
                          For ALDEx2, exactly 2 groups must be specified.
  --exclude-unmapped      Exclude unmapped features from the analysis
```

### 8. Run Statistical Tests

Using metadata for statistics:

```bash
humann3-stats --abundance-file ./processed_data/pathway_abundance-cpm_unstratified.tsv \
    --metadata-file metadata.csv \
    --output-dir ./statistical_tests \
    --group-col Treatment \
    --alpha 0.05
```

### 9. Generate Visualizations

Create visualizations using metadata groups:

```bash
humann3-viz --abundance-file ./processed_data/pathway_abundance-cpm_unstratified.tsv \
    --metadata-file metadata.csv \
    --output-dir ./visualizations \
    --group-col Treatment \
    --shape-col Group \
    --pca --heatmap --barplot \
    --format svg
```

## File Pattern Matching

When using metadata, you can specify how to find files with these options:

1. **Explicit file columns**: Use `--r1-col` and `--r2-col` to specify columns containing file paths
2. **File pattern**: Use `--file-pattern "{sample}_S*_L001_R*.fastq.gz"` with sample placeholder
3. **File suffixes**: Use `--r1-suffix` and `--r2-suffix` (e.g., "_R1.fastq.gz", "_R2.fastq.gz")

## Command Options

### humann3-tools (with metadata options)

Key metadata options:
```
Metadata-driven workflow options:
  --use-metadata          Read samples and file paths from metadata
  --seq-dir SEQ_DIR       Directory containing sequence files
  --sample-col SAMPLE_COL Column name for sample IDs
  --r1-col R1_COL         Column name for R1 sequence file paths
  --r2-col R2_COL         Column name for R2 sequence file paths
  --file-pattern PATTERN  File pattern to match for samples (can use {sample})
  --r1-suffix R1_SUFFIX   Suffix to append for R1 sequence file paths
  --r2-suffix R2_SUFFIX   Suffix to append for R2 sequence file paths
  --samples-file FILE     Tab-delimited file with sample IDs and sequence file paths
```

### humann3-join

```
usage: humann3-join --input-dir INPUT_DIR
                    (--pathabundance | --pathcoverage | --genefamilies)
                    [--output-dir OUTPUT_DIR] [--output-basename OUTPUT_BASENAME]
                    [--units {cpm,relab}]
```

### join_unstratify_humann_output

```
usage: join_unstratify_humann_output --sample-key SAMPLE_KEY
                                     --pathway-dir PATHWAY_DIR
                                     --gene-dir GENE_DIR
                                     [--output-dir OUTPUT_DIR]
                                     [--units {cpm,relab}]
```

### humann3-diff

```
usage: humann3-diff --abundance-file ABUNDANCE_FILE 
                   --metadata-file METADATA_FILE
                   [--group-col GROUP_COL]
                   [--sample-id-col SAMPLE_ID_COL]
                   [--methods METHODS]
```

### humann3-stats

```
usage: humann3-stats --abundance-file ABUNDANCE_FILE 
                    --metadata-file METADATA_FILE
                    [--group-col GROUP_COL]
                    [--sample-id-col SAMPLE_ID_COL]
```

### humann3-viz

```
usage: humann3-viz --abundance-file ABUNDANCE_FILE 
                  --metadata-file METADATA_FILE
                  [--group-col GROUP_COL]
                  [--shape-col SHAPE_COL]
```

## Notes on Sample ID Handling

1. For preprocessing steps, the sample IDs are extracted from sequence file names or specified in metadata.
2. For analysis steps, sample IDs are typically in the column headers of abundance files and should match the IDs in the metadata file.
3. When using `--sample-id-col` in analysis commands, it specifies which column in the metadata contains the same identifiers that are used in the abundance file's columns.
