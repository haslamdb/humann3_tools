#!/bin/bash
#
# conda_setup.sh - Set up Conda environment for HUMAnN3 Tools
#
# This script creates a conda environment with all necessary dependencies
# for running the HUMAnN3 Tools package on Unix/Linux/macOS systems.
#
# Usage:
#   ./conda_setup.sh [--name ENV_NAME] [--python PYTHON_VERSION] 
#                    [--biobakery-channel CHANNEL] [--force]
#
# Examples:
#   # Create default environment (humann3-tools)
#   ./conda_setup.sh
#   
#   # Create environment with specific name and Python version
#   ./conda_setup.sh --name my-humann-env --python 3.9
#   
#   # Force recreation of an existing environment
#   ./conda_setup.sh --force

set -e

# Default values
ENV_NAME="humann3-tools"
PYTHON_VERSION="3.9"
BIOBAKERY_CHANNEL="biobakery"
FORCE=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --name)
      ENV_NAME="$2"
      shift 2
      ;;
    --python)
      PYTHON_VERSION="$2"
      shift 2
      ;;
    --biobakery-channel)
      BIOBAKERY_CHANNEL="$2"
      shift 2
      ;;
    --force)
      FORCE=true
      shift
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: ./conda_setup.sh [--name ENV_NAME] [--python PYTHON_VERSION] [--biobakery-channel CHANNEL] [--force]"
      exit 1
      ;;
  esac
done

# Print banner
echo "================================================================================"
echo "HUMAnN3 Tools Conda Environment Setup"
echo "================================================================================"
echo "Environment name: $ENV_NAME"
echo "Python version: $PYTHON_VERSION"
echo "Biobakery channel: $BIOBAKERY_CHANNEL"
echo "Force recreation: $FORCE"
echo "================================================================================"

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found in PATH. Please install conda first."
    exit 1
fi

echo "Using conda installation at: $(which conda)"

# Create or use existing environment
if conda env list | grep -q "^$ENV_NAME "; then
    if [ "$FORCE" = true ]; then
        echo "Environment '$ENV_NAME' already exists. Removing it as requested..."
        conda env remove --name $ENV_NAME
        conda create --name $ENV_NAME python=$PYTHON_VERSION -y
    else
        echo "Environment '$ENV_NAME' already exists. Use --force to recreate it."
        echo "Proceeding with dependency installation..."
    fi
else
    echo "Creating conda environment: $ENV_NAME with Python $PYTHON_VERSION"
    conda create --name $ENV_NAME python=$PYTHON_VERSION -y
fi

# Install core dependencies from conda-forge
echo "Installing core dependencies..."
conda install --name $ENV_NAME -y -c conda-forge \
    pandas \
    numpy \
    matplotlib \
    seaborn \
    scipy \
    scikit-learn \
    statsmodels \
    tqdm \
    psutil \
    biopython || {
        echo "Warning: Some dependencies failed to install together. Trying individually..."
        for pkg in pandas numpy matplotlib seaborn scipy scikit-learn statsmodels tqdm psutil biopython; do
            conda install --name $ENV_NAME -y -c conda-forge $pkg || echo "Warning: Failed to install $pkg"
        done
    }

# Install biobakery tools
echo "Installing biobakery tools from $BIOBAKERY_CHANNEL channel..."
conda install --name $ENV_NAME -y -c $BIOBAKERY_CHANNEL \
    humann=3.6 \
    kneaddata \
    metaphlan=4.0 || {
        echo "Warning: Some biobakery tools failed to install together. Trying individually..."
        for pkg in "humann=3.6" kneaddata "metaphlan=4.0"; do
            conda install --name $ENV_NAME -y -c $BIOBAKERY_CHANNEL $pkg || echo "Warning: Failed to install $pkg"
        done
    }

# Install pip dependencies
echo "Installing additional dependencies with pip..."
conda run --name $ENV_NAME pip install \
    scikit-posthocs \
    scikit-bio || {
        echo "Warning: Some pip dependencies failed to install together. Trying individually..."
        for pkg in scikit-posthocs scikit-bio; do
            conda run --name $ENV_NAME pip install $pkg || echo "Warning: Failed to install $pkg"
        done
    }

# Install current package in development mode
echo "Installing humann3_tools package..."
conda run --name $ENV_NAME pip install -e . || {
    echo "Warning: Could not install humann3_tools in development mode."
    echo "Please install it manually after activation:"
    echo "    conda activate $ENV_NAME"
    echo "    pip install -e ."
}

# Create activation script
cat > activate_env.sh << EOF
#!/bin/bash
# Activate the $ENV_NAME conda environment
conda activate $ENV_NAME
EOF

chmod +x activate_env.sh
echo "Created activation script: activate_env.sh"

# Print success message and next steps
echo
echo "================================================================================"
echo "HUMAnN3 Tools conda environment '$ENV_NAME' has been set up successfully!"
echo "================================================================================"
echo
echo "To activate the environment, run:"
echo "    conda activate $ENV_NAME"
echo
echo "To verify the installation, run:"
echo "    humann3-tools --help"
echo
echo "To get started with a metadata-driven workflow, run:"
echo "    humann3-tools --run-preprocessing --use-metadata \\"
echo "        --sample-key metadata.csv \\"
echo "        --seq-dir /path/to/sequence/files \\"
echo "        --paired \\"
echo "        --r1-suffix '_R1.fastq.gz' --r2-suffix '_R2.fastq.gz' \\"
echo "        --kneaddata-dbs /path/to/kneaddata_db \\"
echo "        --humann3-nucleotide-db /path/to/chocophlan \\"
echo "        --humann3-protein-db /path/to/uniref \\"
echo "        --output-dir ./results \\"
echo "        --group-col 'Treatment'"
echo
echo "Refer to the usage guide for more examples and options."
echo "================================================================================"
