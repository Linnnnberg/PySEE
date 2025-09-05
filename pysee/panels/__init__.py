"""
PySEE visualization panels.

This module contains all the visualization panel classes for PySEE.
"""

from .base import BasePanel
from .umap import UMAPPanel
from .violin import ViolinPanel
from .heatmap import HeatmapPanel
from .qc import QCPanel
from .dotplot import DotPlotPanel

__all__ = ["BasePanel", "UMAPPanel", "ViolinPanel", "HeatmapPanel", "QCPanel", "DotPlotPanel"]
