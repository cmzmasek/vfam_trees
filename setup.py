import re
from pathlib import Path

from setuptools import setup, find_packages

# Single source of truth: vfam_trees/__init__.py
_init = Path(__file__).parent / "vfam_trees" / "__init__.py"
_m = re.search(r'^__version__\s*=\s*"([^"]+)"', _init.read_text(), re.M)
if not _m:
    raise RuntimeError(f"Could not find __version__ in {_init}")
VERSION = _m.group(1)

setup(
    name="vfam_trees",
    version=VERSION,
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
