# Create the main package folder
import os
from setuptools import setup, find_packages

def create_package_structure():
    """Create the initial package directory structure"""
    
    directories = [
        "humann3_tools",
        "humann3_tools/humann3_tools",
        "humann3_tools/humann3_tools/utils",
        "humann3_tools/humann3_tools/humann3",
        "humann3_tools/humann3_tools/analysis"
    ]
    
    for directory in directories:
        os.makedirs(directory, exist_ok=True)
        # Create __init__.py files in each directory
        if directory.startswith("humann3_tools/humann3_tools"):
            with open(os.path.join(directory, "__init__.py"), 'w') as f:
                f.write("# Package initialization\n")
    
    # Create the root level files
    with open("humann3_tools/setup.py", 'w') as f:
        f.write("""from setuptools import setup, find_packages

setup(
    name="humann3_tools",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "numpy",
        "scipy",
        "scikit-bio",
        "statsmodels",
        "matplotlib",
        "seaborn",
        "scikit-posthocs",
        "scikit-learn",
        # New dependencies for differential abundance
        "matplotlib-venn"  # Optional but useful for method comparisons
    ],
    entry_points={
        'console_scripts': [
            'humann3-tools=humann3_tools.cli:main',
        ],
    },
    author="Your Name",
    author_email="your.email@example.com",
    description="Tools for HUMAnN3 output data processing and analysis",
        )
            """)
        
        
with open("humann3_tools/README.md", 'w') as f:
    f.write("""# HUMAnN3 Tools

A Python package for processing and analyzing HUMAnN3 output data.

## Installation

```bash
pip install -e .

            """)

