#!/usr/bin/env python3
"""
Setup script for humann3_tools package.
"""

from setuptools import setup, find_packages

# Read the long description from README.md
try:
    with open("README.md", "r") as fh:
        long_description = fh.read()
except FileNotFoundError:
    long_description = "HUMAnN3 Tools - A comprehensive toolkit for metagenomic analysis"

# Define package metadata
setup(
    name="humann3_tools",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A comprehensive toolkit for metagenomic analysis with HUMAnN3",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/humann3_tools",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.12",
    install_requires=[
        "numpy>=1.16.0",
        "pandas>=1.0.0",
        "matplotlib>=3.1.0",
        "seaborn>=0.10.0",
        "scikit-learn>=0.22.0",
        "scipy>=1.4.0",
        "statsmodels>=0.11.0",
        "scikit-posthocs>=0.6.0",
    ],
    entry_points={
        "console_scripts": [
            "humann3-tools=humann3_tools.cli.main_cli:main",
            "humann3-kneaddata=humann3_tools.cli.kneaddata_cli:main",
            "humann3-humann3=humann3_tools.cli.humann3_cli:main",
            "humann3-join=humann3_tools.cli.join_cli:main",
            "humann3-stats=humann3_tools.cli.stats_cli:main",
            "humann3-diff=humann3_tools.cli.diff_cli:main",
            "humann3-viz=humann3_tools.cli.viz_cli:main",
        ],
    },
)