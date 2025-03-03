# HUMAnN3 Tools

A comprehensive Python package for assigning raw metagenomic sequence reads to microbial gene and pathway databases using HUMAnN3, followed by downstream processing and analysis.

## Table of Contents
- [Installation](#installation)
    - [Conda Environment Setup](#conda-environment-setup)
    - [Standard Installation](#standard-installation)
- [Overview](#overview)
- [Command-Line Usage](#command-line-usage)
- [Python API Usage](#python-api-usage)
- [HUMAnN3 Processing](#humann3-processing)
- [Downstream Analysis](#downstream-analysis)
- [Differential Abundance Analysis](#differential-abundance-analysis)
- [Visualization](#visualization)
- [Statistical Testing](#statistical-testing)
- [Working with Unmapped Reads](#working-with-unmapped-reads)
- [Advanced Configuration](#advanced-configuration)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)
- [Getting Help](#getting-help)

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

---

## Command-Line Usage

The package provides a command-line tool `humann3-tools` that can be used to run the full analysis pipeline.

### End-to-End Pipeline (Raw Sequences to Analysis)

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

### Analysis of Existing HUMAnN Output Files (pathway-abundance and gene-family files)

```bash
# Full pipeline with defaults
humann3-tools --sample-key /path/to/SampleKey.csv \
    --pathway-dir /path/to/PathwayAbundance \
    --gene-dir /path/to/GeneFamilies \
    --output-dir /path/to/Output \
    --group-col "Group"
```


### Common Options

```bash
# Skip pathway or gene processing
humann3-tools --sample-key samples.csv --pathway-dir pathways/ --gene-dir genes/ \
    --output-dir results/ --skip-pathway

# Skip downstream analysis (only process HUMAnN3 files)
humann3-tools --sample-key samples.csv --pathway-dir pathways/ --gene-dir genes/ \
    --output-dir results/ --skip-downstream

# Use a specific column for grouping in statistical tests
humann3-tools --sample-key samples.csv --pathway-dir pathways/ --gene-dir genes/ \
    --output-dir results/ --group-col "Treatment"

# Enable differential abundance analysis
humann3-tools --sample-key samples.csv --pathway-dir pathways/ --gene-dir genes/ \
    --output-dir results/ --run-diff-abundance
```

### Parallel Processing for Large Datasets

```bash
humann3-tools --run-preprocessing --input-fastq reads_*.fastq --paired \
    --kneaddata-dbs /path/to/kneaddata_db1 /path/to/kneaddata_db2 /path/to/kneaddata_db3 \
    --humann3-nucleotide-db /path/to/chocophlan \
    --humann3-protein-db /path/to/uniref \
    --sample-key /path/to/metadata.csv \
    --output-dir /path/to/output \
    --group-col "Group" \
    --use-parallel \
    --threads-per-sample 4 \
    --max-parallel 8
```

### Run Differential Abundance Analysis

```bash
humann3-tools --sample-key /path/to/metadata.csv \
    --pathway-dir /path/to/pathways \
    --gene-dir /path/to/genes \
    --output-dir /path/to/output \
    --run-diff-abundance \
    --diff-methods aldex2,ancom,ancom-bc \
    --group-col "Group"
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

## HUMAnN3 Processing

HUMAnN3 Tools automates several key processing steps for HUMAnN3 output files:

- **Normalization**: Converts abundance values to counts per million (CPM) using HUMAnN3's `humann_renorm_table` utility  
- **Joining**: Combines files from multiple samples into a single table using `humann_join_tables`  
- **Stratification**: Splits tables into stratified and unstratified versions using `humann_split_stratified_table`

The output directory structure will look like:

```
output_dir/
├── pathways/
│   ├── ProcessedFiles/
│   │   ├── Normalized/
│   │   │   ├── sample1_pathabundance-cpm.tsv
│   │   │   ├── sample2_pathabundance-cpm.tsv
│   │   │   └── ...
│   │   ├── ProcessedFiles_pathabundance-cpm.tsv
│   │   ├── ProcessedFiles_pathabundance-cpm_stratified.tsv
│   │   ├── ProcessedFiles_pathabundance-cpm_unstratified.tsv
│   │   └── pathway_abundance.tsv 
│   └── ...
├── genes/
│   ├── ProcessedFiles/
│   │   ├── Normalized/
│   │   │   ├── sample1_genefamilies-cpm.tsv
│   │   │   ├── sample2_genefamilies-cpm.tsv
│   │   │   └── ...
│   │   ├── ProcessedFiles_genefamilies-cpm.tsv
│   │   ├── ProcessedFiles_genefamilies-cpm_stratified.tsv
│   │   └── ProcessedFiles_genefamilies-cpm_unstratified.tsv
│   └── ...
└── ...
```

---

## Downstream Analysis

The downstream analysis module performs several types of analyses on the processed HUMAnN3 data:

- **Principal Component Analysis (PCA)**: Visualizes similarities and differences between samples  
- **Bar Plots**: Shows abundance patterns across different groups  
- **Statistical Tests**: Performs Kruskal-Wallis tests followed by Dunn's post-hoc tests  
- **Differential Abundance Analysis**: Identifies significantly different features between groups (if enabled)

Results are organized in the `DownstreamAnalysis` directory:

```
output_dir/
├── DownstreamAnalysis/
│   ├── gene_families_bar.svg
│   ├── gene_families_pca.svg
│   ├── pathways_bar.svg
│   ├── pathways_pca.svg
│   ├── kruskal_wallis_pathways.csv
│   ├── dunn_posthoc_tests/
│   │   ├── dunn_pathway1.csv
│   │   ├── dunn_pathway2.csv
│   │   └── ...
│   └── ...
└── ...
```

---

## Differential Abundance Analysis

HUMAnN3 Tools includes implementations of three popular differential abundance methods for microbiome data:

- **ALDEx2**: Uses a Dirichlet-multinomial model to account for compositional data; works only with two-group comparisons  
- **ANCOM**: Analysis of Composition of Microbiomes, based on log-ratio testing; works with multiple groups  
- **ANCOM-BC**: ANCOM with bias correction, addressing limitations of the original ANCOM; works with multiple groups

### Running Differential Abundance Analysis

```bash
# Basic usage
humann3-tools --sample-key samples.csv --pathway-dir pathways/ --gene-dir genes/ \
    --output-dir results/ --group-col Group --run-diff-abundance

# Specify methods and exclude unmapped features
humann3-tools --sample-key samples.csv --pathway-dir pathways/ --gene-dir genes/ \
    --output-dir results/ --group-col Group --run-diff-abundance \
    --diff-methods aldex2,ancom --exclude-unmapped
```

### Output Files

Differential abundance analysis produces the following files in `output_dir/DifferentialAbundance/`:

- `Pathways/aldex2_results.csv`: Results from ALDEx2 analysis on pathways  
- `Pathways/ancom_results.csv`: Results from ANCOM on pathways  
- `Pathways/ancom_bc_results.csv`: Results from ANCOM-BC on pathways  
- `Pathways/method_comparison.txt`: Comparison of significant features across methods  
- `Pathways/venn_diagram.png`: Venn diagram showing overlap of significant features  
- `Pathways/aldex2_volcano.png`: Volcano plot of ALDEx2 results  
- `Pathways/ancom_top_features.png`: Bar plot of top features by ANCOM W-ratio  

Similar files are produced for gene families in the `Genes/` subdirectory.

### Interpreting Results

- **ALDEx2**: Features with *q-value* < 0.05 are considered significant  
- **ANCOM**: Features with *W-ratio* > 0.7 are considered significant  
- **ANCOM-BC**: Features with *q-value* < 0.05 are considered significant

---

## Visualization

HUMAnN3 Tools generates several visualization types:

- **PCA Plots**: Show the relationships between samples based on their functional profiles  
- **Bar Plots**: Display the mean abundance of pathways or gene families across different groups  
- **Volcano Plots** (ALDEx2): Relationship between effect size (x-axis) and significance (y-axis)  
- **Venn Diagrams** (multiple methods): Overlap of significant features among different methods

---

## Statistical Testing

### Kruskal-Wallis Test

A non-parametric alternative to ANOVA for comparing the abundance of features across multiple groups. Results include:

- Test statistic  
- p-value  
- Adjusted p-value (q-value) using Benjamini-Hochberg FDR correction  
- Binary indicator of significance (Reject_H0)

### Dunn's Post-hoc Test

For features with significant Kruskal-Wallis results, a Dunn's post-hoc test identifies which specific group pairs differ significantly.

---

## Working with Unmapped Reads

HUMAnN3 output often includes `UNMAPPED` reads, i.e., sequences not assigned to known functional categories. These can significantly affect differential abundance analysis if they vary between groups.

### Options for Handling Unmapped Reads

1. **Include Unmapped Reads (Default)**  
   ```bash
   humann3-tools --sample-key samples.csv --pathway-dir pathways/ --gene-dir genes/ \
       --output-dir results/ --run-diff-abundance
   ```

2. **Exclude Unmapped Reads**  
   ```bash
   humann3-tools --sample-key samples.csv --pathway-dir pathways/ --gene-dir genes/ \
       --output-dir results/ --run-diff-abundance --exclude-unmapped
   ```

**Recommended Approach**:  
If the proportion of unmapped reads varies significantly between groups, it might be biologically meaningful. Often, you should:

- First analyze with unmapped reads included  
- Then analyze with unmapped reads excluded  
- Compare results to identify robust findings

---

## Advanced Configuration

### Sample Key Format

The sample key CSV file should contain:

- A column with sample identifiers matching the file names in the input directories  
- Additional columns for grouping and other metadata

**Example**:
```csv
SampleName,Group,Treatment,Site
sample1,Control,Placebo,Gut
sample2,Treatment,Drug,Gut
sample3,Control,Placebo,Skin
```

### Resource Management

Optimize memory and CPU usage:

```bash
humann3-tools --run-preprocessing --input-fastq reads_*.fastq \
    --max-memory 32000 \        # Limit memory usage to 32GB
    --threads-per-sample 4 \    # Threads per sample
    --max-parallel 6            # Max samples to process in parallel
```

### Non-interactive Mode

For automated pipelines or batch processing, use `--no-interactive` to skip user prompts and auto-select the sample identifier column:

```bash
humann3-tools --sample-key samples.csv --pathway-dir pathways/ --gene-dir genes/ \
    --output-dir results/ --no-interactive
```

### Custom Logging

Configure logging level and output file:

```bash
humann3-tools --sample-key samples.csv --pathway-dir pathways/ --gene-dir genes/ \
    --output-dir results/ --log-level DEBUG --log-file analysis.log
```

---

# Metadata-Driven Workflow

HUMAnN3 Tools now supports metadata-driven workflows that automatically locate and use sequence files based on sample information in your metadata file. This eliminates the need to manually specify each input file on the command line.

## Using Metadata to Find Sequence Files

You can use the sample information in your metadata file to automatically find and process sequence files:

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

This command will:
1. Read your metadata file to get sample IDs
2. Find sequence files in the specified directory using the provided suffixes
3. Process all found samples through the pipeline

## Metadata with File Paths

If your metadata already contains columns with file paths:

```bash
humann3-tools --run-preprocessing --use-metadata \
    --sample-key /path/to/metadata.csv \
    --seq-dir /path/to/sequence/files \
    --sample-col "SampleID" --r1-col "ForwardRead" --r2-col "ReverseRead" --paired \
    --kneaddata-dbs /path/to/kneaddata_db1 /path/to/kneaddata_db2 \
    --humann3-nucleotide-db /path/to/chocophlan \
    --humann3-protein-db /path/to/uniref \
    --output-dir /path/to/output
```

In this case, the tool will use the file paths specified in the "ForwardRead" and "ReverseRead" columns of your metadata file.

## Using File Patterns

You can also specify a pattern to find files:

```bash
humann3-tools --run-preprocessing --use-metadata \
    --sample-key /path/to/metadata.csv \
    --seq-dir /path/to/sequence/files \
    --file-pattern "{sample}_S*_L001_R*_001.fastq.gz" --paired \
    --kneaddata-dbs /path/to/kneaddata_db1 /path/to/kneaddata_db2 \
    --output-dir /path/to/output
```

Replace `{sample}` with the sample ID from your metadata file to generate file paths.

## Using a Pre-Generated Samples File

For more control, you can generate a samples file separately and use it with the workflow:

```bash
# First, generate the samples file
python -m humann3_tools.utils.metadata_collect \
    --metadata /path/to/metadata.csv \
    --seq-dir /path/to/sequence/files \
    --r1-suffix "_R1.fastq.gz" --r2-suffix "_R2.fastq.gz" --paired \
    --output samples.txt

# Then use it in the workflow
humann3-tools --run-preprocessing \
    --samples-file samples.txt --paired \
    --kneaddata-dbs /path/to/kneaddata_db1 /path/to/kneaddata_db2 \
    --output-dir /path/to/output
```

The samples file is a tab-delimited file with sample IDs and their corresponding file paths.

## Metadata Format Example

Your metadata file should include a column for sample identifiers and optionally columns for file paths:

```csv
SampleID,Group,Treatment,ForwardRead,ReverseRead
sample1,Control,Placebo,sample1_R1.fastq.gz,sample1_R2.fastq.gz
sample2,Treatment,Drug,sample2_R1.fastq.gz,sample2_R2.fastq.gz
sample3,Control,Placebo,sample3_R1.fastq.gz,sample3_R2.fastq.gz
```

The tool will automatically detect common sample ID column names (SampleID, Sample_ID, SampleName, etc.) if not specified.

## Python API for Metadata-Driven Workflow

```python
from humann3_tools.utils.metadata_utils import collect_samples_from_metadata
from humann3_tools import run_preprocessing_and_analysis

# Collect samples from metadata
samples = collect_samples_from_metadata(
    metadata_file="metadata.csv",
    seq_dir="/path/to/sequence/files",
    r1_suffix="_R1.fastq.gz",
    r2_suffix="_R2.fastq.gz",
    paired=True
)

# Flatten the list of files for processing
input_files = []
for sample_files in samples.values():
    input_files.extend(sample_files)

# Run the analysis pipeline
pathway_file, gene_file, success = run_preprocessing_and_analysis(
    input_fastq=input_files,
    sample_key="metadata.csv",
    output_dir="results",
    paired=True,
    threads=8,
    kneaddata_dbs=["/path/to/kneaddata_db1", "/path/to/kneaddata_db2"],
    nucleotide_db="/path/to/chocophlan",
    protein_db="/path/to/uniref",
    group_col="Group"
)
```

## Examples

### Example 1: Basic Analysis Workflow

```bash
# Process HUMAnN3 files and run standard downstream analysis
humann3-tools --sample-key metadata.csv \
    --pathway-dir raw_pathways/ \
    --gene-dir raw_genes/ \
    --output-dir results/ \
    --group-col "DiseaseStatus"
```

### Example 2: End-to-End Pipeline Using Metadata-Driven Workflow

```bash
# Process raw sequence files using sample information from metadata.csv
humann3-tools --run-preprocessing --use-metadata \
    --sample-key metadata.csv \
    --seq-dir /path/to/sequence/files \
    --r1-suffix "_R1.fastq.gz" --r2-suffix "_R2.fastq.gz" --paired \
    --kneaddata-dbs /path/to/human_db /path/to/contaminants_db \
    --humann3-nucleotide-db /path/to/chocophlan \
    --humann3-protein-db /path/to/uniref \
    --output-dir results/ \
    --group-col "DiseaseStatus" \
    --threads 8
```

### Example 3: Processing a Large Dataset with Parallel Execution

```bash
# Automatically identify and process samples in parallel from metadata
humann3-tools --run-preprocessing --use-metadata \
    --sample-key metadata.csv \
    --seq-dir /path/to/sequence/files \
    --file-pattern "{sample}_S*_L001_R*_001.fastq.gz" --paired \
    --kneaddata-dbs /path/to/human_db /path/to/mouse_db \
    --humann3-nucleotide-db /path/to/chocophlan \
    --humann3-protein-db /path/to/uniref \
    --output-dir results/ \
    --group-col "Group" \
    --use-parallel \
    --threads-per-sample 4 \
    --max-parallel 6
```

### Example 4: Focused Differential Abundance Analysis

```bash
# Run only specific differential abundance methods on pathways, excluding unmapped reads
humann3-tools --sample-key metadata.csv \
    --pathway-dir raw_pathways/ \
    --gene-dir raw_genes/ \
    --output-dir results/ \
    --skip-gene \
    --group-col "DiseaseStatus" \
    --run-diff-abundance \
    --diff-methods aldex2,ancom-bc \
    --exclude-unmapped
```

### Example 5: Python API for Custom Analysis with Metadata

```python
import pandas as pd
from humann3_tools import process_humann3_files_only, run_pathway_differential_abundance
from humann3_tools.utils.metadata_utils import collect_samples_from_metadata

# First, collect sample information from metadata
samples = collect_samples_from_metadata(
    metadata_file="metadata.csv",
    seq_dir="/path/to/sequence/files",
    r1_suffix="_R1.fastq.gz",
    r2_suffix="_R2.fastq.gz",
    paired=True
)

print(f"Found {len(samples)} samples with valid sequence files")

# Process the raw HUMAnN3 files
pathway_file, gene_file = process_humann3_files_only(
    sample_key="metadata.csv",
    pathway_dir="pathways/",
    gene_dir="genes/",
    output_dir="processed/"
)

# Run a custom differential abundance analysis
# For a subgroup of samples
metadata = pd.read_csv("metadata.csv")
responders = metadata[metadata["Response"] == "Yes"]["SampleID"].tolist()

# Create a filtered metadata file
responder_metadata = metadata[metadata["SampleID"].isin(responders)]
responder_metadata.to_csv("responders.csv", index=False)

# Run differential abundance on just the responders
results = run_pathway_differential_abundance(
    pathway_file=pathway_file,
    sample_key="responders.csv",
    output_dir="responder_analysis/",
    group_col="Timepoint",
    methods=["ancom", "ancom-bc"]
)

# Print the top significant features (ANCOM example)
if "ancom" in results:
    top_features = results["ancom"].sort_values("W_ratio", ascending=False).head(10)
    print("Top features by ANCOM:")
    for _, row in top_features.iterrows():
        print(f"- {row['feature']}: W-ratio = {row['W_ratio']:.2f}")
```

### Example 6: Processing Samples from Pre-Generated Sample File

```bash
# First, generate a sample file with sequence file paths
python -m humann3_tools.utils.metadata_collect \
    --metadata metadata.csv \
    --seq-dir /path/to/sequence/files \
    --r1-suffix "_R1.fastq.gz" --r2-suffix "_R2.fastq.gz" --paired \
    --output samples.txt \
    --interactive

# Then run the pipeline using the samples file
humann3-tools --run-preprocessing \
    --samples-file samples.txt \
    --paired \
    --kneaddata-dbs /path/to/human_db /path/to/vector_db \
    --humann3-nucleotide-db /path/to/chocophlan \
    --humann3-protein-db /path/to/uniref \
    --output-dir results/ \
    --run-diff-abundance \
    --group-col "Group"
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
   - Consider running the pipeline in stages.  
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
