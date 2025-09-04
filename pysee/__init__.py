"""
PySEE - Interactive, Reproducible Bioinformatics Visualization for Python

PySEE brings iSEE-style linked dashboards to the Python bioinformatics ecosystem.
It provides interactive visualization panels for AnnData objects with linked selections
and reproducible code export capabilities.
"""

__version__ = "0.1.0"
__author__ = "PySEE Contributors"

# Core imports
from .core.data import AnnDataWrapper
from .core.dashboard import PySEE

# Panel imports
from .panels.base import BasePanel
from .panels.umap import UMAPPanel
from .panels.violin import ViolinPanel

__all__ = [
    "PySEE",
    "AnnDataWrapper",
    "BasePanel",
    "UMAPPanel",
    "ViolinPanel",
]
