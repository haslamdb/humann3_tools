#!/usr/bin/env python3
# setup.py for humann3_tools

from setuptools import setup, find_packages

setup(
    name="humann3_tools",
    version="0.2.0",
    description="Tools for processing and analyzing HUMAnN3 output data",
    author="Bioinformatics Team",
    packages=find_packages(),
    install_requires=[
        "pandas>=1.0.0",
        "numpy>=1.18.0",
        "matplotlib>=3.1.0",
        "seaborn>=0.10.0",
        "scikit-learn>=0.22.0",
        "scipy>=1.4.0",
        "scikit-posthocs>=0.6.0",
        "statsmodels>=0.11.0",
        "tqdm>=4.42.0",
        "psutil>=5.7.0",
        "logging>=0.4.9"
    ],
    entry_points={
        'console_scripts': [
            'humann3-join=humann3_tools.humann3.join_cli:main',
            'humann3-diff=humann3_tools.analysis.differential_abundance_cli:main',
            'humann3-stats=humann3_tools.analysis.statistical_cli:main',
            'humann3-viz=humann3_tools.analysis.visualize_cli:main',
        ],
    },
    python_requires='>=3.6',
)
