# HUMAnN3 Tools

A comprehensive Python package for assigning raw metagenomic sequence reads to microbial gene and pathway databases using HUMAnN3, followed by downstream processing and analysis.

## Table of Contents
- [Installation](#installation)
- [Workflow Overview](#workflow-overview)
- [Complete Workflow](#complete-workflow)
- [Step-by-Step Workflow](#step-by-step-workflow)
  - [Step 1: KneadData Processing](#step-1-kneaddata-processing)
  - [Step 2: HUMAnN3 Analysis](#step-2-humann3-analysis) 
  - [Step 3: Preprocessing (Combined Steps 1 & 2)](#step-3-preprocessing-combined-steps-1--2)
  - [Step 4: Join and Normalize](#step-4-join-and-normalize)
  - [Step 5: Statistical Testing](#step-5-statistical-testing)
  - [Step 6: Differential Abundance Analysis](#step-6-differential-abundance-analysis)
  - [Step 7: Visualization](#step-7-visualization)
- [Metadata-Driven Workflow](#metadata-driven-workflow)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)
- [Getting Help](#getting-help)
- [License](#license)
- [Citation](#citation)

---

## Installation

HUMAnN3 Tools offers multiple installation methods to suit different user preferences and environments.

### Conda Environment Setup (Recommended)

It is **highly recommended** to use a dedicated Conda environment to manage the dependencies for `humann3_tools`.

#### Option 1: Using the Python Setup Script

```bash
# Clone the repository
git clone https://github.com/dhaslam/humann3_tools.git

# Navigate to package directory
cd humann3_tools 

# Run the setup script
python conda_setup.py

# Activate the environment
conda activate humann3-tools

# Verify the installation
humann3-tools --help
```

#### Option 2: Using the Shell Script (Unix/Linux/macOS)

```bash
# Clone the repository
git clone https://github.com/dhaslam/humann3_tools.git
cd humann3_tools

# Make the script executable and run it
chmod +x conda_setup.sh
./conda_setup.sh

# Activate the environment
conda activate humann3-tools
```

#### For Customized Installation

Both setup scripts support various options:
```bash
# Example with custom environment name and Python version
python conda_setup.py --name my-humann-env --python 3.8
# Or with shell script
./conda_setup.sh --name my-humann-env --python 3.8
```

### Pip Installation

If you prefer not to use Conda, you can install directly using pip:

```bash
# Clone the repository
git clone https://github.com/dhaslam/humann3_tools.git
cd humann3_tools

# Install in development mode
pip install -e .
```

**Note**: When using pip, you'll need to manually install HUMAnN3, KneadData, and MetaPhlAn3 from the biobakery channel.

### Detailed Installation Instructions

For detailed installation instructions, troubleshooting tips, and advanced options, see the [Installation Guide](install-guide.md).

---

## Workflow Overview

HUMAnN3 Tools provides two ways to analyze your metagenomic data:

1. **Complete Workflow**: Run the entire analysis pipeline from raw reads to visualizations
2. **Step-by-Step Workflow**: Run individual steps separately, with the ability to skip or modify steps as needed

The overall workflow consists of these steps:

1. **KneadData Processing**: Quality control and host depletion 
2. **HUMAnN3 Analysis**: Process cleaned sequences through HUMAnN3
3. **Join & Normalize**: Combine, normalize, and split HUMAnN3 output files
4. **Statistical Testing**: Perform statistical tests across groups
5. **Differential Abundance**: Apply methods like ALDEx2, ANCOM, and ANCOM-BC
6. **Visualization**: Generate plots and figures for publication

![Workflow Diagram](workflow_diagram.png)

---

## Complete Workflow

The `humann3-tools` command provides a way to run the full analysis pipeline in a single command:

```bash
humann3-tools --run-preprocessing \
    --input-fastq reads_1.fastq reads_2.fastq \
    --paired \
    --kneaddata-dbs /path/to/kneaddata_db \
    --humann3-nucleotide-db /path/to/chocophlan \
    --humann3-protein-db /path/to/uniref \
    --sample-key metadata.csv \
    --output-dir results_directory \
    --group-col "Treatment" \
    --run-diff-abundance \
    --threads 8
```

### Key Parameters for Complete Workflow

| Parameter | Description |
|-----------|-------------|
| `--run-preprocessing` | Flag to run KneadData and HUMAnN3 on raw sequence files |
| `--input-fastq` | Input FASTQ file(s) for preprocessing |
| `--paired` | Flag indicating input files are paired-end reads |
| `--kneaddata-dbs` | Path(s) to KneadData reference database(s) |
| `--humann3-nucleotide-db` | Path to HUMAnN3 nucleotide database (ChocoPhlAn) |
| `--humann3-protein-db` | Path to HUMAnN3 protein database (UniRef) |
| `--sample-key` | Path to CSV file with sample metadata |
| `--output-dir` | Directory where output files will be saved |
| `--group-col` | Column name to use for grouping in stats (default: 'Group') |
| `--run-diff-abundance` | Flag to run differential abundance analysis |
| `--threads` | Number of threads to use (default: 1) |

For a complete list of options, run:
```bash
humann3-tools --help
```

---

## Step-by-Step Workflow

If you prefer to run each step individually, or need to rerun specific steps, you can use the specialized commands for each part of the workflow.

### Step 1: KneadData Processing

The `humann3-kneaddata` command allows you to run quality control and host depletion on raw sequence files:

```bash
humann3-kneaddata --input-fastq reads_1.fastq reads_2.fastq \
    --paired \
    --reference-dbs /path/to/kneaddata_db \
    --output-dir kneaddata_output \
    --threads 8 \
    --decontaminate-pairs strict
```

Key options:
- `--input-fastq`: Input FASTQ file(s)
- `--paired`: Flag if input files are paired-end reads
- `--reference-dbs`: Path to reference database(s) for decontamination
- `--output-dir`: Directory for output files
- `--threads`: Number of threads to use
- `--decontaminate-pairs`: Method for paired reads (strict, lenient, or unpaired)

For detailed options:
```bash
humann3-kneaddata --help
```

### Step 2: HUMAnN3 Analysis

The `humann3-run` command performs functional profiling on cleaned sequence files:

```bash
humann3-run --input-files kneaddata_output/sample_paired_1.fastq \
    --nucleotide-db /path/to/chocophlan \
    --protein-db /path/to/uniref \
    --output-dir humann3_output \
    --threads 8
```

Key options:
- `--input-files`: KneadData output file(s)
- `--nucleotide-db`: Path to nucleotide database
- `--protein-db`: Path to protein database
- `--output-dir`: Directory for output files
- `--threads`: Number of threads to use

For detailed options:
```bash
humann3-run --help
```

### Step 3: Preprocessing (Combined Steps 1 & 2)

If you want to run both KneadData and HUMAnN3 together but still separately from downstream analysis, use `humann3-preprocess`:

```bash
humann3-preprocess --input-fastq reads_1.fastq reads_2.fastq \
    --paired \
    --kneaddata-dbs /path/to/kneaddata_db \
    --humann3-nucleotide-db /path/to/chocophlan \
    --humann3-protein-db /path/to/uniref \
    --output-dir preprocessing_output \
    --threads 8
```

Key options:
- Combines options from both `humann3-kneaddata` and `humann3-run`
- `--skip-kneaddata`: Skip KneadData and use existing clean reads
- `--kneaddata-output-files`: Existing KneadData output files when skipping KneadData

For detailed options:
```bash
humann3-preprocess --help
```

### Step 4: Join and Normalize

After running HUMAnN3, the next step is to join, normalize and unstratify the output files. Two commands are available:

#### Option A: Process a specific file type with `humann3-join`:

```bash
humann3-join --input-dir humann3_output \
    --pathabundance \
    --output-dir joined_output \
    --units cpm
```

Key options:
- `--input-dir`: Directory with HUMAnN3 output files
- File type: `--pathabundance`, `--pathcoverage`, or `--genefamilies`
- `--output-dir`: Directory for output files
- `--units`: Units for normalization (cpm or relab)

For detailed options:
```bash
humann3-join --help
```

#### Option B: Process both pathway and gene files with `join_unstratify_humann_output`:

```bash
join_unstratify_humann_output \
    --sample-key metadata.csv \
    --pathway-dir humann3_output/pathabundance \
    --gene-dir humann3_output/genefamilies \
    --output-dir joined_output \
    --units cpm
```

Key options:
- `--sample-key`: Path to metadata CSV file
- `--pathway-dir`: Directory with pathabundance files
- `--gene-dir`: Directory with genefamilies files
- `--output-dir`: Directory for output files
- `--units`: Units for normalization (cpm or relab)

For detailed options:
```bash
join_unstratify_humann_output --help
```

### Step 5: Statistical Testing

The `humann3-stats` command performs statistical tests (Kruskal-Wallis and Dunn's post-hoc tests) on processed HUMAnN3 files:

```bash
humann3-stats --abundance-file joined_output/pathway_abundance-cpm_unstratified.tsv \
    --metadata-file metadata.csv \
    --group-col Treatment \
    --output-dir statistical_results \
    --feature-type pathway
```

Key options:
- `--abundance-file`: Path to unstratified abundance file
- `--metadata-file`: Path to metadata CSV file
- `--group-col`: Column name for grouping samples
- `--output-dir`: Directory for output files
- `--feature-type`: Type of features (pathway or gene)

For detailed options:
```bash
humann3-stats --help
```

### Step 6: Differential Abundance Analysis

The `humann3-diff` command performs differential abundance analysis using methods that account for the compositional nature of microbiome data:

```bash
humann3-diff --abundance-file joined_output/pathway_abundance-cpm_unstratified.tsv \
    --metadata-file metadata.csv \
    --group-col Treatment \
    --output-dir differential_abundance \
    --methods aldex2,ancom,ancom-bc \
    --feature-type pathway
```

Key options:
- `--abundance-file`: Path to unstratified abundance file
- `--metadata-file`: Path to metadata CSV file
- `--group-col`: Column name for grouping samples
- `--methods`: Comma-separated list of methods
- `--output-dir`: Directory for output files
- `--feature-type`: Type of features (pathway or gene)

For detailed options:
```bash
humann3-diff --help
```

### Step 7: Visualization

The `humann3-viz` command creates visualizations from processed HUMAnN3 data:

```bash
humann3-viz --abundance-file joined_output/pathway_abundance-cpm_unstratified.tsv \
    --metadata-file metadata.csv \
    --output-dir visualizations \
    --group-col Treatment \
    --pca --heatmap --barplot
```

Key options:
- `--abundance-file`: Path to unstratified abundance file
- `--metadata-file`: Path to metadata CSV file
- `--group-col`: Column name for coloring points
- `--output-dir`: Directory for output files
- Plot selection: `--pca`, `--heatmap`, `--barplot`, `--abundance-hist`
- `--feature`: Generate boxplot for a specific feature

For detailed options:
```bash
humann3-viz --help
```

---

## Metadata-Driven Workflow

HUMAnN3 Tools supports metadata-driven workflows that automatically locate and use sequence files based on sample information in your metadata file:

```bash
humann3-tools --run-preprocessing --use-metadata \
    --sample-key metadata.csv \
    --seq-dir /path/to/sequence/files \
    --r1-suffix "_R1.fastq.gz" --r2-suffix "_R2.fastq.gz" --paired \
    --kneaddata-dbs /path/to/kneaddata_db \
    --humann3-nucleotide-db /path/to/chocophlan \
    --humann3-protein-db /path/to/uniref \
    --output-dir results_directory \
    --group-col "Treatment"
```

Key options:
- `--use-metadata`: Use metadata file to locate sequence files
- `--seq-dir`: Directory containing sequence files
- `--sample-col`: Column name for sample IDs (autodetected if not specified)
- `--r1-suffix`, `--r2-suffix`: Suffixes for paired-end files

---

## Examples

### Example 1: Complete Workflow

Note that metadata-driven inputs can be used instead of --input-fastq with the following flags:
  --use-metadata \
  --sample-key metadata.csv \

```bash
# Run the full pipeline from raw reads to visualization
humann3-tools --run-preprocessing \
    --input-fastq sample1_R1.fastq.gz sample1_R2.fastq.gz sample2_R1.fastq.gz sample2_R2.fastq.gz \
    --paired \
    --kneaddata-dbs /path/to/human_db /path/to/contaminants_db \
    --humann3-nucleotide-db /path/to/chocophlan \
    --humann3-protein-db /path/to/uniref \
    --sample-key metadata.csv \
    --output-dir full_analysis_results \
    --group-col "Treatment" \
    --run-diff-abundance \
    --threads 8
```

### Example 2: Step-by-Step Workflow

```bash
# Step 1: KneadData processing
humann3-kneaddata --input-fastq sample1_R1.fastq.gz sample1_R2.fastq.gz \
    --paired \
    --reference-dbs /path/to/human_db \
    --output-dir kneaddata_output \
    --threads 8

# Step 2: HUMAnN3 analysis
humann3-run --input-files kneaddata_output/sample1_paired_1.fastq kneaddata_output/sample1_paired_2.fastq \
    --nucleotide-db /path/to/chocophlan \
    --protein-db /path/to/uniref \
    --output-dir humann3_output \
    --threads 8

# Step 3: Join and normalize
humann3-join --input-dir humann3_output \
    --pathabundance \
    --output-dir joined_output \
    --units cpm

# Step 4: Statistical testing
humann3-stats --abundance-file joined_output/pathabundance_cpm_unstratified.tsv \
    --metadata-file metadata.csv \
    --group-col Treatment \
    --output-dir statistical_results

# Step 5: Differential abundance
humann3-diff --abundance-file joined_output/pathabundance_cpm_unstratified.tsv \
    --metadata-file metadata.csv \
    --group-col Treatment \
    --methods aldex2,ancom,ancom-bc \
    --output-dir differential_abundance

# Step 6: Visualization
humann3-viz --abundance-file joined_output/pathabundance_cpm_unstratified.tsv \
    --metadata-file metadata.csv \
    --group-col Treatment \
    --pca --heatmap --barplot \
    --output-dir visualizations
```

### Example 3: Skip KneadData and Start with HUMAnN3

```bash
# If you already have quality-filtered reads:
humann3-run --input-files clean_reads/sample1.fastq clean_reads/sample2.fastq \
    --nucleotide-db /path/to/chocophlan \
    --protein-db /path/to/uniref \
    --output-dir humann3_output \
    --threads 8
```

### Example 4: Skip Preprocessing and Start with Existing HUMAnN3 Outputs

```bash
# Only join and normalize existing HUMAnN3 outputs
join_unstratify_humann_output \
    --sample-key metadata.csv \
    --pathway-dir existing_outputs/pathways \
    --gene-dir existing_outputs/genes \
    --output-dir joined_output \
    --units cpm
```

---

## Troubleshooting

### Common Issues

1. **Missing Input Files**  
   - Check file paths and naming patterns.  
   - Use `--list-files` to see what files are detected.  

2. **Sample Key Issues**  
   - Ensure sample identifiers in the CSV match the file names.  
   - Check for duplicate sample IDs.  
   - Verify CSV encoding (use UTF-8).

3. **KneadData Errors**
   - Ensure reference databases are properly built and indexed.
   - For paired-end issues, try different `--decontaminate-pairs` options.

4. **HUMAnN3 Errors**
   - Ensure nucleotide and protein databases are correctly installed.
   - Check that input files are in the correct format.

5. **Statistical Test Errors**  
   - ALDEx2 requires exactly two groups.  
   - Check for missing values in abundance data.

6. **Memory Issues with Large Datasets**  
   - Run steps separately to manage memory usage.  
   - Use `--threads` to control CPU usage.

---

## Getting Help

If you encounter issues not covered in this documentation, please:

- Check the log file for detailed error messages  
- Set `--log-level DEBUG` for more verbose output  
- Open an issue on the GitHub repository with a description of the problem and relevant log entries

---

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use HUMAnN3 Tools in your research, please cite:

- The original HUMAnN3 paper: Franzosa EA, et al. (2018). Species-level functional profiling of metagenomes and metatranscriptomes. Nature Methods, 15(11), 962-968.
- KneadData: The Huttenhower Lab (https://github.com/biobakery/kneaddata)
- This tool: Haslam, D. (2025). HUMAnN3 Tools: A comprehensive framework for metagenomic analysis.
