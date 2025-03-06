# HUMAnN3 Tools

A comprehensive Python package for assigning raw metagenomic sequence reads to microbial gene and pathway databases using HUMAnN3, followed by downstream processing and analysis.

## Table of Contents
- [Installation](#installation)
    - [Conda Environment Setup](#conda-environment-setup)
    - [Standard Installation](#standard-installation)
- [Overview](#overview)
- [Command-Line Usage](#command-line-usage)
    - [Main Pipeline](#main-pipeline)
    - [Modular Commands](#modular-commands)
    - [Required and Optional Parameters](#required-and-optional-parameters)
- [Modular Analysis Commands](#modular-analysis-commands)
    - [Preprocessing Pipeline (`humann3-preprocess`)](#preprocessing-pipeline)
    - [Join and Unstratify (`humann3-join` and `join_unstratify`)](#join-and-unstratify)
    - [Statistical Testing (`humann3-stats`)](#statistical-testing)
    - [Differential Abundance Analysis (`humann3-dif`)](#differential-abundance-analysis)
    - [Visualizations (`humann3-viz`)](#visualizations)
- [Python API Usage](#python-api-usage)
- [Metadata-Driven Workflow](#metadata-driven-workflow)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)
- [Getting Help](#getting-help)
- [License](#license)
- [Citation](#citation)

---

## Installation

This section describes how to install and setup the `humann3_tools` package.

### Conda Environment Setup

It is **highly recommended** to use a Conda environment to manage the dependencies for `humann3_tools`. This ensures a clean and isolated environment, avoiding potential conflicts with other Python packages you might have installed.

**Steps:**

1.  **Navigate to the Package Directory:**
    Before running the setup script, make sure you are in the top-level directory of the `humann3_tools` package (the directory containing `conda_setup.py`, `README.md`, `pyproject.toml` etc). If you cloned the repository, this would be:

    ```bash
    cd humann3_tools 
    ```

2.  **Run the Setup Script:**
    Use the provided `conda_setup.py` script to create and activate a new Conda environment, and then install the `humann3_tools` package and its dependencies.

    ```bash
    python conda_setup.py
    ```

    This script will:
    *   Create a new Conda environment named `humann3_tools` (you can modify this in the script if needed).
    *   Install all required and recommended dependencies into this environment.
    *   Install the `humann3_tools` package in development mode (`-e .`).
    *   Activate the new environment

3.  **Activate the Environment (if not done automatically):**
    If the environment is not activated automatically, you can activate it manually:

    ```bash
    conda activate humann3_tools
    ```
    
4. **Verify the installation**
    You can then check the installation by using:
    ```bash
    humann3-tools --help
    ```

**Note:** You only need to run `conda_setup.py` once. To use `humann3_tools` in the future, just activate the environment using `conda activate humann3_tools`.

### Standard Installation

If you prefer not to use Conda, you can install directly using pip:

```bash
# Clone the repository
git clone https://github.com/dhaslam/humann3_tools.git
cd humann3_tools

# Install in development mode
pip install -e .
```

Or install directly using pip (once published to PyPI):
```bash
pip install humann3-tools
```

**Required dependencies** (will install automatically if not present):
- pandas  
- numpy  
- scipy  
- scikit-bio  
- scikit-learn  
- scikit-posthocs  
- statsmodels  
- matplotlib  
- seaborn  
- matplotlib-venn (optional, for visualization of method comparisons)

---

## Overview

**HUMAnN3 Tools** provides an end-to-end solution for metagenomic analysis:

1. **Preprocessing**: Run quality control and host depletion with KneadData
2. **Functional Profiling**: Process raw sequences through HUMAnN3
3. **HUMAnN3 Processing**: Normalize, join, and split HUMAnN3 output files
4. **Downstream Analysis**: Perform statistical tests, PCA, and visualization
5. **Differential Abundance Analysis**: Apply methods like ALDEx2, ANCOM, and ANCOM-BC

The package can be used as either an integrated pipeline or through modular commands that can be run separately for more focused analyses.

---

## Command-Line Usage

### Main Pipeline

The package provides a command-line tool `humann3-tools` that can be used to run the full analysis pipeline.

```bash
humann3-tools --run-preprocessing --input-fastq reads_1.fastq reads_2.fastq --paired \
    --kneaddata-dbs /path/to/kneaddata_db1 /path/to/kneaddata_db2 /path/to/kneaddata_db3 \
    --humann3-nucleotide-db /path/to/chocophlan \
    --humann3-protein-db /path/to/uniref \
    --sample-key /path/to/metadata.csv \
    --output-dir /path/to/output \
    --group-col "Group" \
    --threads 8
```

### Modular Commands

The package now includes several modular commands that allow for more focused and efficient analysis:

1. **`humann3-preprocess`**: Run KneadData and HUMAnN3 on raw sequence files
2. **`humann3-join`**: Join, normalize, and unstratify HUMAnN3 output files
3. **`join_unstratify_humann_output`**: Simplified version for joining and unstratifying HUMAnN3 outputs
4. **`humann3-stats`**: Run statistical tests on processed HUMAnN3 files
5. **`humann3-dif`**: Perform differential abundance analysis
6. **`humann3-viz`**: Generate visualizations from HUMAnN3 output files

### Required and Optional Parameters

Below is a list of key parameters for the `humann3-tools` command, indicating which are required and which are optional:

| Parameter | Required? | Description |
|-----------|-----------|-------------|
| `--sample-key` | Required | Path to CSV file with sample metadata |
| `--output-dir` | Required | Directory where output files will be saved |
| `--pathway-dir` | Required* | Directory containing HUMAnN3 pathway abundance files |
| `--gene-dir` | Required* | Directory containing HUMAnN3 gene family files |
| `--input-fastq` | Required** | Input FASTQ file(s) when running preprocessing |
| `--run-preprocessing` | Optional | Flag to run KneadData and HUMAnN3 on raw sequence files |
| `--kneaddata-dbs` | Required*** | Path(s) to KneadData reference database(s) |
| `--humann3-nucleotide-db` | Optional | Path to HUMAnN3 nucleotide database (ChocoPhlAn) |
| `--humann3-protein-db` | Optional | Path to HUMAnN3 protein database (UniRef) |
| `--paired` | Optional | Flag indicating input files are paired-end reads |
| `--decontaminate-pairs` | Optional | Select level of decontamination of paired reads (default 'strict') |
| `--threads` | Optional | Number of threads to use (default: 1) |
| `--group-col` | Optional | Column name to use for grouping in stats (default: 'Group') |
| `--skip-pathway` | Optional | Skip HUMAnN3 pathway processing |
| `--skip-gene` | Optional | Skip HUMAnN3 gene family processing |
| `--skip-downstream` | Optional | Skip downstream analysis |
| `--run-diff-abundance` | Optional | Run differential abundance analysis |
| `--diff-methods` | Optional | Comma-separated list of methods to use (default: aldex2,ancom,ancom-bc) |
| `--exclude-unmapped` | Optional | Exclude unmapped features from differential abundance analysis |
| `--use-parallel` | Optional | Use parallel processing for preprocessing steps |
| `--threads-per-sample` | Optional | Number of threads to use per sample (with parallel processing) |
| `--max-parallel` | Optional | Maximum number of samples to process in parallel |
| `--log-file` | Optional | Path to combined log file |
| `--log-level` | Optional | Logging level (default: INFO) |
| `--no-interactive` | Optional | Non-interactive mode for sample key column selection |
| `--join-only` | Optional | Run only join and unstratify operations |
| `--units` | Optional | Units for normalization (default: cpm, can be relab) |

\* Required for processing existing HUMAnN3 output files  
\*\* Required when `--run-preprocessing` is used  
\*\*\* Required when `--run-preprocessing` is used

---

## Modular Analysis Commands

HUMAnN3 Tools now includes several specialized commands to run specific parts of the pipeline independently for greater flexibility.

### Preprocessing Pipeline

The `humann3-preprocess` command provides a standalone tool for running KneadData and HUMAnN3 preprocessing steps on raw sequence files.

```bash
# Basic usage with paired-end reads
humann3-preprocess --input-fastq sample1_R1.fastq.gz sample1_R2.fastq.gz --paired \
    --kneaddata-dbs /path/to/kneaddata_db \
    --humann3-nucleotide-db /path/to/chocophlan \
    --humann3-protein-db /path/to/uniref \
    --output-dir /path/to/output \
    --threads 8
```

**Key Features:**
- Run only KneadData and HUMAnN3 steps without downstream analysis
- Process multiple samples in one command
- Skip KneadData and use existing outputs
- Detailed logging and progress reporting
- Organize outputs for easy downstream analysis

### Join and Unstratify

Two commands are available for joining and processing HUMAnN3 output files:

#### 1. `humann3-join` Command

The `humann3-join` command is specialized for joining, normalizing, and unstratifying a specific type of HUMAnN3 output file.

```bash
humann3-join --input-dir /path/to/humann3_outputs \
    --pathabundance \
    --output-dir /path/to/output \
    --units cpm
```

**Key Features:**
- Focus on one file type (pathabundance, pathcoverage, or genefamilies)
- Normalize to CPM or relative abundance
- Generate unstratified outputs automatically
- Simple and fast processing

#### 2. `join_unstratify_humann_output` Command

This command provides a simplified interface for joining and unstratifying both pathway and gene family files in one step.

```bash
join_unstratify_humann_output \
    --sample-key metadata.csv \
    --pathway-dir path/to/pathabundance/files \
    --gene-dir path/to/genefamilies/files \
    --output-dir processed_output
```

**Key Features:**
- Process both pathway and gene files in a single command
- Connect with sample metadata
- Skip pathway or gene processing optionally
- Choose normalization units (CPM or relative abundance)
- Non-interactive mode for automated workflows

### Statistical Testing

The `humann3-stats` command performs statistical tests on processed HUMAnN3 data to identify significant differences between groups. 

```bash
humann3-stats --abundance-file pathway_abundance_unstratified.tsv \
    --metadata-file metadata.csv \
    --group-col Treatment \
    --output-dir statistical_results \
    --feature-type pathway
```

**Key Features:**
- Run Kruskal-Wallis tests across groups
- Perform Dunn's post-hoc tests for significant features
- Apply multiple testing correction (Benjamini-Hochberg FDR)
- Detailed CSV outputs for results
- Support for both pathway and gene data

### Differential Abundance Analysis

The `humann3-dif` command implements advanced methods for differential abundance analysis that account for the compositional nature of microbiome data.

```bash
humann3-dif --abundance-file pathway_abundance_unstratified.tsv \
    --metadata-file metadata.csv \
    --group-col Treatment \
    --output-dir differential_abundance_results \
    --methods aldex2,ancom,ancom-bc \
    --feature-type pathway
```

**Key Features:**
- Multiple method support: ALDEx2, ANCOM, ANCOM-BC, and Kruskal-Wallis
- Method comparison with Venn diagrams (requires matplotlib-venn)
- Volcano plot for ALDEx2 results
- Top features visualization for ANCOM
- Option to exclude unmapped reads from analysis

### Visualizations

The `humann3-viz` command creates publication-quality visualizations from processed HUMAnN3 data.

```bash
humann3-viz --abundance-file pathway_abundance_unstratified.tsv \
    --metadata-file metadata.csv \
    --output-dir visualizations \
    --group-col Treatment \
    --pca --heatmap --barplot
```

**Key Features:**
- PCA plots for sample relationships
- Heatmaps of top features across samples
- Bar plots for comparing group abundance patterns
- Abundance distribution histograms
- Feature-specific boxplots
- Multiple output formats (SVG, PNG, PDF)

For detailed information about available visualizations:
```bash
humann3-viz --help-info
```

---

## Python API Usage

You can also use the Python API for more flexibility and integration with your own scripts.

### Run End-to-End Preprocessing and Analysis

```python
from humann3_tools import run_preprocessing_and_analysis

pathway_file, gene_file, success = run_preprocessing_and_analysis(
    input_fastq=["reads_1.fastq", "reads_2.fastq"],
    sample_key="metadata.csv",
    output_dir="results",
    paired=True,
    threads=8,
    kneaddata_db="/path/to/kneaddata_db",
    nucleotide_db="/path/to/chocophlan",
    protein_db="/path/to/uniref",
    group_col="Group"
)
```

### Process and Analyze Pathway_abundance and Genefamilies Output Files

```python
from humann3_tools import run_full_pipeline

pathway_file, gene_file, success = run_full_pipeline(
    sample_key="/path/to/sample_key.csv",
    pathway_dir="/path/to/pathway_dir",
    gene_dir="/path/to/gene_dir",
    output_dir="/path/to/output",
    group_col="Group",
    run_diff_abundance=True,
    log_file="humann3_analysis.log"
)

if success:
    print("Analysis completed successfully!")
    print(f"Pathway file: {pathway_file}")
    print(f"Gene family file: {gene_file}")
```

### Processing HUMAnN3 Files Only

```python
from humann3_tools import process_humann3_files_only

pathway_file, gene_file = process_humann3_files_only(
    sample_key="/path/to/sample_key.csv",
    pathway_dir="/path/to/pathway_dir",
    gene_dir="/path/to/gene_dir",
    output_dir="/path/to/output",
    log_file="humann3_processing.log"
)
```

### Running Downstream Analysis on Existing Files

```python
from humann3_tools import analyze_existing_humann3_files

success = analyze_existing_humann3_files(
    pathway_file="/path/to/pathway_abundance.tsv",
    gene_file="/path/to/gene_families_unstratified.tsv",
    sample_key="/path/to/sample_key.csv",
    output_dir="/path/to/analysis_results",
    group_col="Treatment",
    log_file="downstream_analysis.log"
)
```

### Running Only Differential Abundance Analysis

```python
from humann3_tools import run_pathway_differential_abundance

results = run_pathway_differential_abundance(
    pathway_file="/path/to/pathway_abundance.tsv",
    sample_key="/path/to/sample_key.csv",
    output_dir="/path/to/results",
    group_col="Treatment",
    methods=["aldex2", "ancom-bc"],
    include_unmapped=False
)

# Access the results
if 'aldex2' in results:
    significant = results['aldex2'][results['aldex2']['q_value'] < 0.05]
    print(f"Found {len(significant)} significant features with ALDEx2")
```

---

## Metadata-Driven Workflow

HUMAnN3 Tools supports metadata-driven workflows that automatically locate and use sequence files based on sample information in your metadata file. This eliminates the need to manually specify each input file on the command line.

```bash
humann3-tools --run-preprocessing --use-metadata \
    --sample-key /path/to/metadata.csv \
    --seq-dir /path/to/sequence/files \
    --r1-suffix "_R1.fastq.gz" --r2-suffix "_R2.fastq.gz" --paired \
    --kneaddata-dbs /path/to/kneaddata_db1 /path/to/kneaddata_db2 \
    --humann3-nucleotide-db /path/to/chocophlan \
    --humann3-protein-db /path/to/uniref \
    --output-dir /path/to/output \
    --group-col "Group"
```

See the full documentation for more information on metadata-driven workflows.

---

## Examples

### Example 1: Using Modular Commands for a Complete Analysis

```bash
# Step 1: Preprocess raw sequence files
humann3-preprocess --input-fastq sample1_R1.fastq.gz sample1_R2.fastq.gz --paired \
    --kneaddata-dbs /path/to/kneaddata_db \
    --humann3-nucleotide-db /path/to/chocophlan \
    --humann3-protein-db /path/to/uniref \
    --output-dir preprocessing_output \
    --threads 8

# Step 2: Join and unstratify HUMAnN3 outputs
join_unstratify_humann_output \
    --sample-key metadata.csv \
    --pathway-dir preprocessing_output/humann3_output/PathwayAbundance \
    --gene-dir preprocessing_output/humann3_output/GeneFamilies \
    --output-dir processed_output

# Step 3: Run statistical tests
humann3-stats \
    --abundance-file processed_output/pathways/ProcessedFiles/pathway_abundance-cpm_unstratified.tsv \
    --metadata-file metadata.csv \
    --group-col Treatment \
    --output-dir statistical_results

# Step 4: Run differential abundance analysis
humann3-dif \
    --abundance-file processed_output/pathways/ProcessedFiles/pathway_abundance-cpm_unstratified.tsv \
    --metadata-file metadata.csv \
    --group-col Treatment \
    --methods aldex2,ancom,ancom-bc \
    --output-dir differential_abundance_results

# Step 5: Create visualizations
humann3-viz \
    --abundance-file processed_output/pathways/ProcessedFiles/pathway_abundance-cpm_unstratified.tsv \
    --metadata-file metadata.csv \
    --group-col Treatment \
    --pca --heatmap --barplot \
    --output-dir visualizations
```

### Example 2: End-to-End Analysis with the Main Pipeline

```bash
# Complete analysis using the main pipeline
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

### Example 3: Processing Existing HUMAnN3 Output Files

```bash
# Process existing HUMAnN3 output files
humann3-tools --sample-key metadata.csv \
    --pathway-dir existing_outputs/pathways \
    --gene-dir existing_outputs/genes \
    --output-dir processed_results \
    --group-col "Group" \
    --run-diff-abundance
```

### Example 4: Join and Unstratify Only

```bash
# Just join and unstratify existing HUMAnN3 outputs
humann3-tools --join-only \
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

3. **Statistical Test Errors**  
   - Minimum number of samples per group may be required.  
   - ALDEx2 requires exactly two groups.  
   - Check for missing values in abundance data.

4. **Memory Issues with Large Datasets**  
   - Consider running the pipeline in stages using the modular commands.  
   - Use `--skip-pathway` or `--skip-gene` to process fewer data types at once.

### Error Messages and Solutions

- **"No valid pathway/gene files found for any samples"**  
  Check naming patterns and directory structure; use `--list-files`.

- **"No shared samples between abundance data and metadata"**  
  Sample IDs in metadata do not match file names. Check for case sensitivity.

- **"This implementation only supports two groups for comparison"**  
  ALDEx2 is two-group only; try ANCOM or ANCOM-BC for multiple groups.

- **"No significant pathways found after FDR correction"**  
  - Try a less stringent significance threshold.  
  - Exclude unmapped reads with `--exclude-unmapped`.  
  - Use alternative methods.

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
