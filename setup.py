from setuptools import setup, find_packages

setup(
    name="vfam_trees",
    version="1.1.2",
    packages=find_packages(),
    install_requires=[
        "biopython>=1.81",
        "click>=8.1",
        "pyyaml>=6.0",
        "snakemake>=7.0",
        "requests>=2.31",
        "matplotlib>=3.9",
    ],
    entry_points={
        "console_scripts": [
            "vfam_trees=vfam_trees.cli:main",
        ],
    },
    python_requires=">=3.9",
)
