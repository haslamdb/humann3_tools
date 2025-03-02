#!/bin/bash
set -e

# Create and activate environment
conda create -n humann3_tools python=3.12 -y
eval "$(conda shell.bash hook)"
conda activate humann3_tools

# Add channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels biobakery

# Install core dependencies
conda install -y humann=3.8 kneaddata
conda install -y pandas numpy scipy scikit-bio scikit-learn statsmodels matplotlib seaborn
conda install -y -c conda-forge scikit-posthocs matplotlib-venn tqdm psutil

# Install the package
# navigate to the main package directory if necessary
# cd humann3_analysis
pip install --ignore-installed --no-cache-dir -e .

# Verify installation
echo "Installation complete. Testing command availability:"
which humann3-tools
humann3-tools --help

# Optional: Download test databases
echo "Do you want to download test databases? (y/n)"
read answer
if [ "$answer" = "y" ]; then
  echo "Enter directory for databases:"
  read db_dir
  mkdir -p $db_dir
  humann_databases --download chocophlan test $db_dir
  humann_databases --download uniref test $db_dir
  kneaddata_database --download test_database_human $db_dir
  echo "Test databases downloaded to $db_dir"
fi

echo "Setup complete!"