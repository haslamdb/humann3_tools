HUMAnN3 Tools
A comprehensive Python package for preprocessing, processing, and analyzing metagenomic sequence data using KneadData, HUMAnN3, and downstream analysis tools.

Overview
HUMAnN3 Tools provides an end-to-end solution for metagenomic analysis:

Preprocessing: Run quality control and host depletion with KneadData
Functional Profiling: Process raw sequences through HUMAnN3
HUMAnN3 Processing: Normalize, join, and split HUMAnN3 output files
Downstream Analysis: Perform statistical tests, PCA, and visualization
Differential Abundance Analysis: Apply methods like ALDEx2, ANCOM, and ANCOM-BC
Installation
Install from GitHub:


Copy
git clone https://github.com/yourusername/humann3_tools.git
cd humann3_tools
pip install -e .
Prerequisites
KneadData and HUMAnN3 must be installed and in your PATH
Reference databases for KneadData and HUMAnN3 (can be specified at runtime)
Dependencies
The package automatically installs Python dependencies including:

pandas, numpy, scipy
scikit-bio, scikit-learn, scikit-posthocs
statsmodels, matplotlib, seaborn
matplotlib-venn (optional, for visualization of method comparisons)
Command-Line Usage
End-to-End Pipeline (Raw Sequences to Analysis)

Copy
humann3-tools --run-preprocessing --input-fastq reads_1.fastq reads_2.fastq --paired \
    --kneaddata-db /path/to/kneaddata_db \
    --humann3-nucleotide-db /path/to/chocophlan \
    --humann3-protein-db /path/to/uniref \
    --sample-key /path/to/metadata.csv \
    --output-dir /path/to/output \
    --group-col "Group" \
    --threads 8
Parallel Processing for Large Datasets

Copy
humann3-tools --run-preprocessing --input-fastq reads_*.fastq --paired \
    --kneaddata-db /path/to/kneaddata_db \
    --humann3-nucleotide-db /path/to/chocophlan \
    --humann3-protein-db /path/to/uniref \
    --sample-key /path/to/metadata.csv \
    --output-dir /path/to/output \
    --group-col "Group" \
    --use-parallel \
    --threads-per-sample 4 \
    --max-parallel 8
Process Existing HUMAnN3 Output Files

Copy
humann3-tools --sample-key /path/to/metadata.csv \
    --pathway-dir /path/to/pathways \
    --gene-dir /path/to/genes \
    --output-dir /path/to/output \
    --group-col "Treatment"
Run Differential Abundance Analysis

Copy
humann3-tools --sample-key /path/to/metadata.csv \
    --pathway-dir /path/to/pathways \
    --gene-dir /path/to/genes \
    --output-dir /path/to/output \
    --run-diff-abundance \
    --diff-methods aldex2,ancom,ancom-bc \
    --group-col "Group"
Python API Usage
Run End-to-End Preprocessing and Analysis
python

Copy
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
Process Existing HUMAnN3 Output Files
python

Copy
from humann3_tools import run_full_pipeline

pathway_file, gene_file, success = run_full_pipeline(
    sample_key="metadata.csv",
    pathway_dir="pathway_files/",
    gene_dir="gene_files/",
    output_dir="results/",
    group_col="Treatment"
)
Run Only Differential Abundance Analysis
python

Copy
from humann3_tools import run_pathway_differential_abundance

results = run_pathway_differential_abundance(
    pathway_file="path/to/pathway_abundance.tsv",
    sample_key="metadata.csv",
    output_dir="results/",
    group_col="Group",
    methods=["aldex2", "ancom"],
    include_unmapped=False
)
Output Directory Structure

Copy
output_dir/
├── PreprocessedData/           # Raw data preprocessing results
│   ├── kneaddata_output/       # KneadData results
│   └── humann3_output/         # Direct HUMAnN3 output
├── pathways/                   # HUMAnN3 pathway processing
│   └── ProcessedFiles/
│       ├── Normalized/
│       ├── ProcessedFiles_pathabundance-cpm.tsv
│       ├── ProcessedFiles_pathabundance-cpm_stratified.tsv
│       ├── ProcessedFiles_pathabundance-cpm_unstratified.tsv
│       └── pathway_abundance.tsv
├── genes/                      # HUMAnN3 gene processing
│   └── ProcessedFiles/
│       ├── Normalized/
│       ├── ProcessedFiles_genefamilies-cpm.tsv
│       ├── ProcessedFiles_genefamilies-cpm_stratified.tsv
│       └── ProcessedFiles_genefamilies-cpm_unstratified.tsv
├── DownstreamAnalysis/         # Statistical analysis & visualization
│   ├── gene_families_bar.svg
│   ├── gene_families_pca.svg
│   ├── pathways_bar.svg
│   ├── pathways_pca.svg
│   ├── kruskal_wallis_pathways.csv
│   └── dunn_posthoc_tests/
└── DifferentialAbundance/      # Differential abundance results
    ├── Pathways/
    │   ├── aldex2_results.csv
    │   ├── ancom_results.csv
    │   ├── ancom_bc_results.csv
    │   ├── venn_diagram.png
    │   └── method_comparison.txt
    └── Genes/
Sample Key Format
The sample key CSV file should contain:

A column with sample identifiers matching the file names in the input directories
Additional columns for grouping and metadata information
Example:


Copy
SampleName,Group,Treatment,Site
sample1,Control,Placebo,Gut
sample2,Treatment,Drug,Gut
sample3,Control,Placebo,Skin
Advanced Features
Resource Management
Optimize memory and CPU usage:


Copy
humann3-tools --run-preprocessing --input-fastq reads_*.fastq \
    --max-memory 32000 \        # Limit memory usage to 32GB
    --threads-per-sample 4 \    # Threads per sample
    --max-parallel 6            # Max samples to process in parallel
Non-interactive Mode
For automated pipelines:


Copy
humann3-tools --sample-key samples.csv --pathway-dir pathways/ --gene-dir genes/ \
    --output-dir results/ --no-interactive
Handling Unmapped Reads
Control how unmapped reads are treated in differential abundance analysis:


Copy
humann3-tools --sample-key samples.csv --pathway-dir pathways/ --gene-dir genes/ \
    --output-dir results/ --run-diff-abundance --exclude-unmapped
Troubleshooting
Common Issues
Database Paths: Ensure KneadData and HUMAnN3 databases are correctly specified
Sample Key Format: Check that sample identifiers match file naming patterns
Memory Usage: For large datasets, use --max-memory and adjust parallel processing parameters
Missing Results: View log files with --log-level DEBUG --log-file run.log
Getting Help
For detailed information about parameters:


Copy
humann3-tools --help
Citation
If you use HUMAnN3 Tools in your research, please cite:

The original HUMAnN3 paper: Franzosa EA, et al. (2018). Species-level functional profiling of metagenomes and metatranscriptomes. Nature Methods, 15(11), 962-968.
KneadData: The Huttenhower Lab (https://github.com/biobakery/kneaddata)
This tool: Haslam, D. (2025). HUMAnN3 Tools: A comprehensive framework for metagenomic analysis.
License
This project is licensed under the MIT License - see the LICENSE file for details.




Retry

