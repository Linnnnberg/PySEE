"""
Heatmap visualization panel for PySEE.

This module provides the HeatmapPanel class for visualizing gene expression matrices
with hierarchical clustering, dendrograms, and interactive selection.
"""

from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import pdist

from ..core.data import AnnDataWrapper
from .base import BasePanel


class HeatmapPanel(BasePanel):
    """
    Panel for visualizing gene expression heatmaps with clustering.

    This panel creates interactive heatmaps of gene expression data with support for
    hierarchical clustering, dendrograms, gene/cell filtering, and linked interactions.

    Parameters
    ----------
    panel_id : str
        Unique identifier for the panel
    genes : List[str], optional
        List of genes to display. If None, uses top variable genes
    cells : List[str], optional
        List of cell IDs to display. If None, uses all cells
    group_by : str, optional
        Column name in adata.obs to group cells by
    title : str, optional
        Display title for the panel
    """

    def __init__(
        self,
        panel_id: str,
        genes: Optional[List[str]] = None,
        cells: Optional[List[str]] = None,
        group_by: Optional[str] = None,
        title: Optional[str] = None,
    ):
        super().__init__(panel_id, title)

        self.set_config("genes", genes)
        self.set_config("cells", cells)
        self.set_config("group_by", group_by)
        self.set_config("n_top_genes", 50)
        self.set_config("cluster_genes", True)
        self.set_config("cluster_cells", True)
        self.set_config("clustering_method", "ward")
        self.set_config("color_scale", "RdBu_r")
        self.set_config("show_dendrograms", True)
        self.set_config("show_colorbar", True)
        self.set_config("gene_labels", True)
        self.set_config("cell_labels", True)
        self.set_config("max_genes", 100)
        self.set_config("max_cells", 500)

    def _check_data_requirements(self) -> bool:
        """
        Check if the data wrapper meets the panel's requirements.

        Returns
        -------
        bool
            True if requirements are met, False otherwise
        """
        if self._data_wrapper is None:
            return False

        # Check if we have expression data
        try:
            expression_matrix = self._data_wrapper.get_expression_data()
            return expression_matrix.size > 0
        except Exception:
            return False

    def _get_expression_data(self) -> Tuple[np.ndarray, List[str], List[str]]:
        """
        Get the expression data for the heatmap.

        Returns
        -------
        expression_matrix : np.ndarray
            Expression matrix (genes x cells)
        gene_names : List[str]
            List of gene names
        cell_names : List[str]
            List of cell names
        """
        if self._data_wrapper is None:
            raise ValueError("No data wrapper set for this panel")

        # Get genes to display
        genes = self.get_config("genes")
        if genes is None:
            # Use top variable genes
            n_top = self.get_config("n_top_genes")
            genes = self._data_wrapper.get_top_variable_genes(n_top)

        # Get cells to display
        cells = self.get_config("cells")
        if cells is None:
            # Use all cells or a subset
            max_cells = self.get_config("max_cells")
            all_cells = self._data_wrapper._adata.obs_names.tolist()
            if len(all_cells) > max_cells:
                # Randomly sample cells
                np.random.seed(42)  # For reproducibility
                cells = np.random.choice(all_cells, max_cells, replace=False).tolist()
            else:
                cells = all_cells

        # Limit genes if too many
        max_genes = self.get_config("max_genes")
        if len(genes) > max_genes:
            genes = genes[:max_genes]

        # Get expression matrix
        expression_matrix = self._data_wrapper.get_gene_expression(genes, cells)

        return expression_matrix, genes, cells

    def _perform_clustering(
        self, expression_matrix: np.ndarray, cluster_genes: bool = True, cluster_cells: bool = True
    ) -> Tuple[
        Optional[np.ndarray], Optional[np.ndarray], Optional[List[int]], Optional[List[int]]
    ]:
        """
        Perform hierarchical clustering on the expression matrix.

        Parameters
        ----------
        expression_matrix : np.ndarray
            Expression matrix (genes x cells)
        cluster_genes : bool
            Whether to cluster genes
        cluster_cells : bool
            Whether to cluster cells

        Returns
        -------
        gene_linkage : np.ndarray, optional
            Linkage matrix for genes
        cell_linkage : np.ndarray, optional
            Linkage matrix for cells
        gene_order : List[int], optional
            Reordered gene indices
        cell_order : List[int], optional
            Reordered cell indices
        """
        method = self.get_config("clustering_method")

        gene_linkage = None
        cell_linkage = None
        gene_order = None
        cell_order = None

        if cluster_genes and expression_matrix.shape[0] > 1:
            try:
                # Cluster genes (rows)
                gene_distances = pdist(expression_matrix, metric="correlation")
                # Check for infinite or NaN values
                if np.isfinite(gene_distances).all():
                    gene_linkage = linkage(gene_distances, method=method)
                    gene_clusters = fcluster(gene_linkage, t=0.7, criterion="distance")
                    gene_order = np.argsort(gene_clusters)
                else:
                    # Fallback: use variance-based ordering
                    gene_variances = np.var(expression_matrix, axis=1)
                    gene_order = np.argsort(gene_variances)[::-1]
            except Exception:
                # Fallback: use variance-based ordering
                gene_variances = np.var(expression_matrix, axis=1)
                gene_order = np.argsort(gene_variances)[::-1]

        if cluster_cells and expression_matrix.shape[1] > 1:
            try:
                # Cluster cells (columns)
                cell_distances = pdist(expression_matrix.T, metric="correlation")
                # Check for infinite or NaN values
                if np.isfinite(cell_distances).all():
                    cell_linkage = linkage(cell_distances, method=method)
                    cell_clusters = fcluster(cell_linkage, t=0.7, criterion="distance")
                    cell_order = np.argsort(cell_clusters)
                else:
                    # Fallback: use variance-based ordering
                    cell_variances = np.var(expression_matrix, axis=0)
                    cell_order = np.argsort(cell_variances)[::-1]
            except Exception:
                # Fallback: use variance-based ordering
                cell_variances = np.var(expression_matrix, axis=0)
                cell_order = np.argsort(cell_variances)[::-1]

        return gene_linkage, cell_linkage, gene_order, cell_order  # type: ignore[return-value]

    def _create_dendrogram_trace(
        self, linkage_matrix: np.ndarray, orientation: str = "top", side: str = "top"
    ) -> go.Scatter:
        """
        Create a dendrogram trace for the heatmap.

        Parameters
        ----------
        linkage_matrix : np.ndarray
            Linkage matrix from hierarchical clustering
        orientation : str
            Orientation of the dendrogram ("top", "bottom", "left", "right")
        side : str
            Which side of the heatmap to place the dendrogram

        Returns
        -------
        go.Scatter
            Dendrogram trace
        """
        # Create dendrogram coordinates
        dendro = dendrogram(linkage_matrix, no_plot=True)

        # Extract coordinates
        icoord = np.array(dendro["icoord"])
        dcoord = np.array(dendro["dcoord"])

        # Flatten coordinates for plotting
        x_coords = []
        y_coords = []

        for i, d in zip(icoord, dcoord):
            x_coords.extend([i[0], i[1], i[2], i[3], None])
            y_coords.extend([d[0], d[1], d[2], d[3], None])

        return go.Scatter(
            x=x_coords,
            y=y_coords,
            mode="lines",
            line=dict(color="black", width=1),
            showlegend=False,
            hoverinfo="skip",
        )

    def render(self) -> go.Figure:
        """
        Render the heatmap visualization.

        Returns
        -------
        go.Figure
            Plotly figure containing the heatmap
        """
        if self._data_wrapper is None:
            raise ValueError("No data wrapper set for this panel")

        # Get expression data
        expression_matrix, gene_names, cell_names = self._get_expression_data()

        if expression_matrix.size == 0:
            # Return empty figure if no data
            fig = go.Figure()
            fig.add_annotation(
                text="No data available for heatmap",
                xref="paper",
                yref="paper",
                x=0.5,
                y=0.5,
                showarrow=False,
                font=dict(size=16),
            )
            return fig

        # Perform clustering
        cluster_genes = self.get_config("cluster_genes")
        cluster_cells = self.get_config("cluster_cells")
        gene_linkage, cell_linkage, gene_order, cell_order = self._perform_clustering(
            expression_matrix, cluster_genes, cluster_cells
        )

        # Reorder data if clustering was performed
        if gene_order is not None:
            expression_matrix = expression_matrix[gene_order]
            gene_names = [gene_names[i] for i in gene_order]

        if cell_order is not None:
            expression_matrix = expression_matrix[:, cell_order]
            cell_names = [cell_names[i] for i in cell_order]

        # Create subplots for dendrograms
        show_dendrograms = self.get_config("show_dendrograms")
        show_gene_dendro = show_dendrograms and cluster_genes and gene_linkage is not None
        show_cell_dendro = show_dendrograms and cluster_cells and cell_linkage is not None

        # Calculate subplot layout
        rows = 2 if show_gene_dendro else 1
        cols = 2 if show_cell_dendro else 1

        # Create subplots
        fig = make_subplots(
            rows=rows,
            cols=cols,
            subplot_titles=([self.title] if self.title else None),
            horizontal_spacing=0.1,
            vertical_spacing=0.1,
            specs=[[{"secondary_y": False}] * cols] * rows,
        )

        # Main heatmap subplot
        heatmap_row = 2 if show_gene_dendro else 1
        heatmap_col = 2 if show_cell_dendro else 1

        # Create heatmap trace
        heatmap_trace = go.Heatmap(
            z=expression_matrix,
            x=cell_names if self.get_config("cell_labels") else None,
            y=gene_names if self.get_config("gene_labels") else None,
            colorscale=self.get_config("color_scale"),
            showscale=self.get_config("show_colorbar"),
            hoverongaps=False,
            hovertemplate="<b>%{y}</b><br>"
            + "Cell: %{x}<br>"
            + "Expression: %{z:.3f}<br>"
            + "<extra></extra>",
        )

        fig.add_trace(heatmap_trace, row=heatmap_row, col=heatmap_col)

        # Add dendrograms if requested
        if show_gene_dendro and gene_linkage is not None:
            dendro_trace = self._create_dendrogram_trace(gene_linkage, "top")
            fig.add_trace(dendro_trace, row=1, col=heatmap_col)

        if show_cell_dendro and cell_linkage is not None:
            dendro_trace = self._create_dendrogram_trace(cell_linkage, "left")
            fig.add_trace(dendro_trace, row=heatmap_row, col=1)

        # Update layout
        fig.update_layout(
            title=self.title or "Gene Expression Heatmap", height=600, showlegend=False
        )

        # Update axes
        if show_gene_dendro:
            fig.update_xaxes(showticklabels=False, row=1, col=heatmap_col)
            fig.update_yaxes(showticklabels=False, row=1, col=heatmap_col)

        if show_cell_dendro:
            fig.update_xaxes(showticklabels=False, row=heatmap_row, col=1)
            fig.update_yaxes(showticklabels=False, row=heatmap_row, col=1)

        # Update main heatmap axes
        fig.update_xaxes(title="Cells", tickangle=45, row=heatmap_row, col=heatmap_col)
        fig.update_yaxes(title="Genes", row=heatmap_row, col=heatmap_col)

        return fig

    def get_selection_code(self) -> str:
        """
        Generate Python code for the current selection.

        Returns
        -------
        str
            Python code snippet
        """
        if self._selection is None:
            return "# No selection in heatmap panel"

        genes = self.get_config("genes")
        cells = self.get_config("cells")

        code_lines = ["# Heatmap panel selection"]

        if genes:
            code_lines.append(f"selected_genes = {genes}")

        if cells:
            code_lines.append(f"selected_cells = {cells}")

        code_lines.append("# Use these for downstream analysis")

        return "\n".join(code_lines)

    def set_genes(self, genes: List[str]) -> None:
        """Set the genes to display in the heatmap."""
        self.set_config("genes", genes)

    def set_cells(self, cells: List[str]) -> None:
        """Set the cells to display in the heatmap."""
        self.set_config("cells", cells)

    def set_clustering_method(self, method: str) -> None:
        """Set the clustering method (ward, complete, average, single)."""
        valid_methods = ["ward", "complete", "average", "single"]
        if method not in valid_methods:
            raise ValueError(f"Clustering method must be one of {valid_methods}")
        self.set_config("clustering_method", method)

    def set_color_scale(self, scale: str) -> None:
        """Set the color scale for the heatmap."""
        self.set_config("color_scale", scale)

    def set_max_genes(self, max_genes: int) -> None:
        """Set the maximum number of genes to display."""
        self.set_config("max_genes", max_genes)

    def set_max_cells(self, max_cells: int) -> None:
        """Set the maximum number of cells to display."""
        self.set_config("max_cells", max_cells)

    def toggle_clustering(self, cluster_genes: bool = None, cluster_cells: bool = None) -> None:
        """Toggle clustering for genes and/or cells."""
        if cluster_genes is not None:
            self.set_config("cluster_genes", cluster_genes)
        if cluster_cells is not None:
            self.set_config("cluster_cells", cluster_cells)

    def toggle_dendrograms(self, show: bool = None) -> None:
        """Toggle display of dendrograms."""
        if show is not None:
            self.set_config("show_dendrograms", show)
