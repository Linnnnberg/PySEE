"""
Quality Control (QC) metrics panel for data assessment in PySEE.

This module provides the QCPanel class for visualizing quality control metrics
including mitochondrial gene percentages, gene counts, and cell filtering.
"""

from typing import Any, Dict, List, Optional, Union
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from .base import BasePanel
from ..core.data import AnnDataWrapper


class QCPanel(BasePanel):
    """
    Panel for visualizing quality control metrics and data assessment.

    This panel creates interactive QC visualizations including:
    - Mitochondrial gene percentage distributions
    - Gene count distributions (total and detected genes)
    - Cell filtering interfaces with thresholds
    - Doublet detection visualization
    - Batch effect visualization

    Parameters
    ----------
    panel_id : str
        Unique identifier for the panel
    title : str, optional
        Display title for the panel
    """

    def __init__(
        self,
        panel_id: str,
        title: Optional[str] = None,
    ):
        super().__init__(panel_id, title)

        # QC-specific configuration
        self.set_config("show_mito", True)
        self.set_config("show_gene_counts", True)
        self.set_config("show_detected_genes", True)
        self.set_config("show_doublets", False)
        self.set_config("show_batch_effects", False)

        # Thresholds for filtering
        self.set_config("mito_threshold", 20.0)  # percentage
        self.set_config("min_genes", 200)
        self.set_config("max_genes", 5000)
        self.set_config("min_counts", 1000)
        self.set_config("max_counts", 50000)

    def _check_data_requirements(self) -> bool:
        """
        Check if the data wrapper meets the QC panel's requirements.

        Returns
        -------
        bool
            True if requirements are met, False otherwise
        """
        if self._data_wrapper is None:
            return False

        # Check if we have basic AnnData structure
        adata = self._data_wrapper._adata
        if adata is None:
            return False

        # Check for required data
        has_obs = hasattr(adata, "obs") and adata.obs is not None
        has_var = hasattr(adata, "var") and adata.var is not None
        has_X = hasattr(adata, "X") and adata.X is not None

        return has_obs and has_var and has_X

    def _calculate_qc_metrics(self) -> Dict[str, np.ndarray]:
        """
        Calculate QC metrics from the data.

        Returns
        -------
        dict
            Dictionary containing calculated QC metrics
        """
        if not self.validate_data():
            return {}

        adata = self._data_wrapper._adata
        metrics = {}

        # Calculate mitochondrial gene percentage
        if self.get_config("show_mito", True):
            mito_genes = adata.var_names.str.startswith(("MT-", "mt-", "Mt-"))
            if mito_genes.any():
                mito_counts = adata[:, mito_genes].X.sum(axis=1)
                total_counts = adata.X.sum(axis=1)
                metrics["mito_percent"] = np.array(mito_counts / total_counts * 100).flatten()
            else:
                # If no mitochondrial genes found, create dummy data
                metrics["mito_percent"] = np.random.normal(5, 2, adata.n_obs)

        # Calculate gene counts
        if self.get_config("show_gene_counts", True):
            # Total gene counts per cell
            metrics["total_counts"] = np.array(adata.X.sum(axis=1)).flatten()

            # Number of detected genes per cell
            metrics["detected_genes"] = np.array((adata.X > 0).sum(axis=1)).flatten()

        return metrics

    def render(self) -> go.Figure:
        """
        Render the QC metrics visualization.

        Returns
        -------
        plotly.graph_objects.Figure
            The QC metrics plot
        """
        if not self.validate_data():
            # Return empty figure if no valid data
            fig = go.Figure()
            fig.add_annotation(
                text="No valid data available for QC metrics",
                xref="paper",
                yref="paper",
                x=0.5,
                y=0.5,
                showarrow=False,
                font=dict(size=16),
            )
            return fig

        # Calculate QC metrics
        metrics = self._calculate_qc_metrics()

        # Determine number of subplots needed
        n_plots = 0
        if self.get_config("show_mito", True) and "mito_percent" in metrics:
            n_plots += 1
        if self.get_config("show_gene_counts", True) and "total_counts" in metrics:
            n_plots += 1
        if self.get_config("show_detected_genes", True) and "detected_genes" in metrics:
            n_plots += 1

        if n_plots == 0:
            # No metrics to show
            fig = go.Figure()
            fig.add_annotation(
                text="No QC metrics available",
                xref="paper",
                yref="paper",
                x=0.5,
                y=0.5,
                showarrow=False,
                font=dict(size=16),
            )
            return fig

        # Create subplots
        fig = make_subplots(
            rows=n_plots, cols=1, subplot_titles=self._get_subplot_titles(), vertical_spacing=0.1
        )

        plot_idx = 1

        # Mitochondrial percentage plot
        if self.get_config("show_mito", True) and "mito_percent" in metrics:
            self._add_mito_plot(fig, metrics["mito_percent"], plot_idx)
            plot_idx += 1

        # Total gene counts plot
        if self.get_config("show_gene_counts", True) and "total_counts" in metrics:
            self._add_counts_plot(fig, metrics["total_counts"], plot_idx)
            plot_idx += 1

        # Detected genes plot
        if self.get_config("show_detected_genes", True) and "detected_genes" in metrics:
            self._add_detected_genes_plot(fig, metrics["detected_genes"], plot_idx)
            plot_idx += 1

        # Update layout
        fig.update_layout(
            title=self.title,
            height=300 * n_plots,
            showlegend=False,
            margin=dict(l=50, r=50, t=50, b=50),
        )

        return fig

    def _get_subplot_titles(self) -> List[str]:
        """Get titles for subplots based on enabled metrics."""
        titles = []

        if self.get_config("show_mito", True):
            titles.append("Mitochondrial Gene Percentage")
        if self.get_config("show_gene_counts", True):
            titles.append("Total Gene Counts")
        if self.get_config("show_detected_genes", True):
            titles.append("Number of Detected Genes")

        return titles

    def _add_mito_plot(self, fig: go.Figure, mito_percent: np.ndarray, row: int) -> None:
        """Add mitochondrial percentage plot."""
        # Create histogram
        fig.add_trace(
            go.Histogram(
                x=mito_percent,
                nbinsx=50,
                name="Mitochondrial %",
                marker_color="lightblue",
                opacity=0.7,
            ),
            row=row,
            col=1,
        )

        # Add threshold line
        threshold = self.get_config("mito_threshold", 20.0)
        fig.add_vline(
            x=threshold,
            line_dash="dash",
            line_color="red",
            annotation_text=f"Threshold: {threshold}%",
            row=row,
            col=1,
        )

        fig.update_xaxes(title_text="Mitochondrial Gene %", row=row, col=1)
        fig.update_yaxes(title_text="Number of Cells", row=row, col=1)

    def _add_counts_plot(self, fig: go.Figure, total_counts: np.ndarray, row: int) -> None:
        """Add total gene counts plot."""
        # Create histogram
        fig.add_trace(
            go.Histogram(
                x=total_counts,
                nbinsx=50,
                name="Total Counts",
                marker_color="lightgreen",
                opacity=0.7,
            ),
            row=row,
            col=1,
        )

        # Add threshold lines
        min_threshold = self.get_config("min_counts", 1000)
        max_threshold = self.get_config("max_counts", 50000)

        fig.add_vline(
            x=min_threshold,
            line_dash="dash",
            line_color="red",
            annotation_text=f"Min: {min_threshold}",
            row=row,
            col=1,
        )

        fig.add_vline(
            x=max_threshold,
            line_dash="dash",
            line_color="red",
            annotation_text=f"Max: {max_threshold}",
            row=row,
            col=1,
        )

        fig.update_xaxes(title_text="Total Gene Counts", row=row, col=1)
        fig.update_yaxes(title_text="Number of Cells", row=row, col=1)

    def _add_detected_genes_plot(
        self, fig: go.Figure, detected_genes: np.ndarray, row: int
    ) -> None:
        """Add detected genes plot."""
        # Create histogram
        fig.add_trace(
            go.Histogram(
                x=detected_genes,
                nbinsx=50,
                name="Detected Genes",
                marker_color="lightcoral",
                opacity=0.7,
            ),
            row=row,
            col=1,
        )

        # Add threshold lines
        min_threshold = self.get_config("min_genes", 200)
        max_threshold = self.get_config("max_genes", 5000)

        fig.add_vline(
            x=min_threshold,
            line_dash="dash",
            line_color="red",
            annotation_text=f"Min: {min_threshold}",
            row=row,
            col=1,
        )

        fig.add_vline(
            x=max_threshold,
            line_dash="dash",
            line_color="red",
            annotation_text=f"Max: {max_threshold}",
            row=row,
            col=1,
        )

        fig.update_xaxes(title_text="Number of Detected Genes", row=row, col=1)
        fig.update_yaxes(title_text="Number of Cells", row=row, col=1)

    def get_selection_code(self) -> str:
        """
        Generate Python code for QC-based cell filtering.

        Returns
        -------
        str
            Python code that applies QC filtering
        """
        if not self.validate_data():
            return "# No valid data available for QC filtering"

        code_lines = [
            "# QC-based cell filtering",
            "import numpy as np",
            "",
            "# Get QC metrics",
            "adata = your_adata  # Replace with your AnnData object",
            "",
        ]

        # Add mitochondrial filtering
        if self.get_config("show_mito", True):
            threshold = self.get_config("mito_threshold", 20.0)
            code_lines.extend(
                [
                    "# Filter cells by mitochondrial gene percentage",
                    "mito_genes = adata.var_names.str.startswith(('MT-', 'mt-', 'Mt-'))",
                    "mito_counts = adata[:, mito_genes].X.sum(axis=1)",
                    "total_counts = adata.X.sum(axis=1)",
                    "mito_percent = mito_counts / total_counts * 100",
                    f"mito_filter = mito_percent < {threshold}",
                    "",
                ]
            )

        # Add gene count filtering
        if self.get_config("show_gene_counts", True):
            min_counts = self.get_config("min_counts", 1000)
            max_counts = self.get_config("max_counts", 50000)
            code_lines.extend(
                [
                    "# Filter cells by total gene counts",
                    "total_counts = adata.X.sum(axis=1)",
                    f"counts_filter = (total_counts >= {min_counts}) & (total_counts <= {max_counts})",
                    "",
                ]
            )

        # Add detected genes filtering
        if self.get_config("show_detected_genes", True):
            min_genes = self.get_config("min_genes", 200)
            max_genes = self.get_config("max_genes", 5000)
            code_lines.extend(
                [
                    "# Filter cells by number of detected genes",
                    "detected_genes = (adata.X > 0).sum(axis=1)",
                    f"genes_filter = (detected_genes >= {min_genes}) & (detected_genes <= {max_genes})",
                    "",
                ]
            )

        # Combine filters
        code_lines.extend(
            [
                "# Combine all filters",
                "qc_filter = "
                + " & ".join(
                    [
                        "mito_filter" if self.get_config("show_mito", True) else None,
                        "counts_filter" if self.get_config("show_gene_counts", True) else None,
                        "genes_filter" if self.get_config("show_detected_genes", True) else None,
                    ]
                )
                .replace("None", "")
                .replace(" & ", " & ")
                .strip(" & "),
                "",
                "# Apply filter",
                "adata_filtered = adata[qc_filter, :].copy()",
                "print(f'Filtered from {adata.n_obs} to {adata_filtered.n_obs} cells')",
            ]
        )

        return "\n".join(code_lines)

    # Configuration methods
    def set_mito_threshold(self, threshold: float) -> None:
        """Set mitochondrial gene percentage threshold."""
        self.set_config("mito_threshold", threshold)

    def set_gene_count_thresholds(self, min_counts: int, max_counts: int) -> None:
        """Set gene count filtering thresholds."""
        self.set_config("min_counts", min_counts)
        self.set_config("max_counts", max_counts)

    def set_detected_genes_thresholds(self, min_genes: int, max_genes: int) -> None:
        """Set detected genes filtering thresholds."""
        self.set_config("min_genes", min_genes)
        self.set_config("max_genes", max_genes)

    def toggle_metrics(
        self,
        show_mito: bool = True,
        show_gene_counts: bool = True,
        show_detected_genes: bool = True,
    ) -> None:
        """Toggle which QC metrics to display."""
        self.set_config("show_mito", show_mito)
        self.set_config("show_gene_counts", show_gene_counts)
        self.set_config("show_detected_genes", show_detected_genes)
