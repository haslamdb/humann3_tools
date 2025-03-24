# HUMAnN3 Tools - Installation Guide

This guide provides detailed instructions for installing the HUMAnN3 Tools package using various methods. Choose the installation method that works best for your system and requirements.

## Prerequisites

Before installing HUMAnN3 Tools, ensure you have the following prerequisites:

### For All Methods:
- Git (to clone the repository)
- Internet connection (to download dependencies)

### For Conda Installation:
- Conda (Anaconda or Miniconda)
- Python 3.8, 3.9, or 3.10 (3.9 recommended)

### For Pip Installation:
- Python 3.8, 3.9, or 3.10 (3.9 recommended)
- pip (Python package installer)

## Installation Methods

HUMAnN3 Tools offers three primary installation methods:

1. **Conda Setup Script (Python)**: Recommended for Windows users or those who prefer Python scripts
2. **Conda Setup Shell Script**: Recommended for Unix/Linux/macOS users
3. **Pip Installation**: For users who prefer to manage dependencies manually or don't use conda

## 1. Conda Setup Script (Python)

This method uses a Python script to create a conda environment with all required dependencies.

```bash
# Clone the repository
git clone https://github.com/dhaslam/humann3_tools.git

# Navigate to the package directory
cd humann3_tools

# Run the setup script
python conda_setup.py
```

### Options

The conda setup script supports several command-line options:

```bash
python conda_setup.py [--name ENV_NAME] [--python PYTHON_VERSION] [--biobakery-channel CHANNEL] [--force]
```

- `--name`: Name of the conda environment (default: humann3-tools)
- `--python`: Python version to use (default: 3.9)
- `--biobakery-channel`: Conda channel for biobakery tools (default: biobakery)
- `--force`: Force recreation of an existing environment

### Examples

```bash
# Create a custom environment with Python 3.8
python conda_setup.py --name my-humann-env --python 3.8

# Recreate an existing environment
python conda_setup.py --force
```

## 2. Conda Setup Shell Script

For Unix/Linux/macOS users, a shell script is provided for convenience.

```bash
# Clone the repository
git clone https://github.com/dhaslam/humann3_tools.git

# Navigate to the package directory
cd humann3_tools

# Make the script executable
chmod +x conda_setup.sh

# Run the setup script
./conda_setup.sh
```

### Options

The shell script supports the same options as the Python script:

```bash
./conda_setup.sh [--name ENV_NAME] [--python PYTHON_VERSION] [--biobakery-channel CHANNEL] [--force]
```

### Examples

```bash
# Create a custom environment with Python 3.8
./conda_setup.sh --name my-humann-env --python 3.8

# Recreate an existing environment
./conda_setup.sh --force
```

## 3. Pip Installation

For users who prefer pip or want to manage dependencies manually:

```bash
# Clone the repository
git clone https://github.com/dhaslam/humann3_tools.git

# Navigate to the package directory
cd humann3_tools

# Install the package and dependencies
pip install -e .

# Or install with development dependencies
pip install -e ".[dev]"
```

### Note on Biobakery Tools

When using pip installation, you'll need to manually install HUMAnN3, KneadData, and MetaPhlAn3. These are typically installed via conda:

```bash
# Install biobakery tools with conda
conda install -c biobakery humann=3.6 kneaddata metaphlan=4.0
```

## Verifying the Installation

After installation, verify that HUMAnN3 Tools is correctly installed:

```bash
# Activate the conda environment (if using conda)
conda activate humann3-tools  # or your custom environment name

# Test the command-line tool
humann3-tools --help
```

You should see the help message for the `humann3-tools` command.

## Troubleshooting

### Common Installation Issues

#### Conda Environment Creation Fails

If you encounter errors during conda environment creation:

1. Check that conda is properly installed and in your PATH
2. Try specifying an older Python version: `--python 3.8`
3. On Windows, run the script from Anaconda Prompt
4. Check that you have administrative privileges if needed

#### Dependency Installation Failures

If some dependencies fail to install:

1. Make sure you have internet connectivity
2. Try updating conda: `conda update -n base conda`
3. Try installing the problematic package manually:
   ```bash
   conda install -n humann3-tools -c conda-forge package_name
   ```

#### Biobakery Tools Installation Issues

If HUMAnN3, KneadData, or MetaPhlAn fail to install:

1. Check the biobakery channel: `conda config --show channels`
2. Add the biobakery channel: `conda config --add channels biobakery`
3. Try installing each tool individually:
   ```bash
   conda install -n humann3-tools -c biobakery humann=3.6
   conda install -n humann3-tools -c biobakery kneaddata
   conda install -n humann3-tools -c biobakery metaphlan=4.0
   ```

### Getting Help

If you encounter issues not covered in this guide:

1. Check the [GitHub repository issues](https://github.com/dhaslam/humann3_tools/issues)
2. Run commands with verbose logging: `--log-level DEBUG`
3. Open a new issue on GitHub with details about your environment and the error messages
