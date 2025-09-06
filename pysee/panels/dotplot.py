"""
Dot plot visualization panel for PySEE.

This module provides the DotPlotPanel class for visualizing marker gene expression
with dot plots showing expression levels and percentage of cells expressing each gene.
"""

from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ..core.data import AnnDataWrapper
from .base import BasePanel


class DotPlotPanel(BasePanel):
    """
    Panel for visualizing marker gene expression with dot plots.

    This panel creates interactive dot plots showing gene expression levels and
    percentage of cells expressing each gene across different groups or clusters.

    Parameters
    ----------
    panel_id : str
        Unique identifier for the panel
    genes : List[str], optional
        List of genes to display. If None, uses top variable genes
    group_by : str, optional
        Column name in adata.obs to group cells by
    title : str, optional
        Display title for the panel
    """

    def __init__(
        self,
        panel_id: str,
        genes: Optional[List[str]] = None,
        group_by: Optional[str] = None,
        title: Optional[str] = None,
    ):
        super().__init__(panel_id, title)

        self.set_config("genes", genes)
        self.set_config("group_by", group_by)
        self.set_config("n_top_genes", 20)
        self.set_config("min_pct", 0.05)  # Minimum percentage of cells expressing gene
        self.set_config("max_pct", 1.0)  # Maximum percentage of cells expressing gene
        self.set_config("min_expr", 0.0)  # Minimum expression level
        self.set_config("max_expr", 5.0)  # Maximum expression level
        self.set_config("dot_size_range", [5, 50])  # Size range for dots
        self.set_config("color_scale", "viridis")
        self.set_config("show_legend", True)
        self.set_config("show_axes", True)
        self.set_config("sort_genes", True)
        self.set_config("sort_groups", True)
        self.set_config("max_genes", 50)
        self.set_config("max_groups", 20)

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
            if expression_matrix is None or expression_matrix.shape[0] == 0:
                return False
        except Exception:
            return False

        # Check if we have genes
        genes = self.get_config("genes")
        if genes is None or len(genes) == 0:
            return False

        # Check if genes exist in the data
        try:
            available_genes = self._data_wrapper.adata.var_names
            if not all(gene in available_genes for gene in genes):
                return False
        except Exception:
            return False

        # Check if group_by column exists (if specified)
        group_by = self.get_config("group_by")
        if group_by is not None:
            try:
                obs_columns = self._data_wrapper.get_obs_columns()
                if group_by not in obs_columns:
                    return False
            except Exception:
                return False

        return True

    def render(self) -> go.Figure:
        """
        Render the dot plot visualization.

        Returns
        -------
        plotly.graph_objects.Figure
            The dot plot figure
        """
        if not self.validate_data():
            return self._create_error_figure("Panel data requirements not met")

        try:
            # Get configuration
            genes = self.get_config("genes")
            group_by = self.get_config("group_by")
            dot_size_range = self.get_config("dot_size_range")
            color_scale = self.get_config("color_scale")
            show_legend = self.get_config("show_legend")
            show_axes = self.get_config("show_axes")
            sort_genes = self.get_config("sort_genes")
            sort_groups = self.get_config("sort_groups")
            max_genes = self.get_config("max_genes")
            max_groups = self.get_config("max_groups")

            # Get data
            expression_matrix = self._data_wrapper.get_expression_data()
            gene_names = self._data_wrapper.adata.var_names
            obs_data = self._data_wrapper.adata.obs

            # Filter genes if too many
            if len(genes) > max_genes:
                genes = genes[:max_groups]

            # Get group information
            if group_by is not None and group_by in obs_data.columns:
                groups = obs_data[group_by].values
                unique_groups = np.unique(groups)
                if len(unique_groups) > max_groups:
                    unique_groups = unique_groups[:max_groups]
            else:
                groups = np.array(["All"] * expression_matrix.shape[0])
                unique_groups = ["All"]

            # Calculate statistics for each gene-group combination
            dot_data = self._calculate_dot_plot_data(
                expression_matrix, gene_names, genes, groups, unique_groups
            )

            # Create the dot plot
            fig = self._create_dot_plot_figure(
                dot_data,
                genes,
                unique_groups,
                dot_size_range,
                color_scale,
                show_legend,
                show_axes,
                sort_genes,
                sort_groups,
            )

            return fig

        except Exception as e:
            return self._create_error_figure(f"Error rendering dot plot: {str(e)}")

    def _calculate_dot_plot_data(
        self,
        expression_matrix: np.ndarray,
        gene_names: List[str],
        genes: List[str],
        groups: np.ndarray,
        unique_groups: np.ndarray,
    ) -> pd.DataFrame:
        """
        Calculate dot plot statistics for each gene-group combination.

        Parameters
        ----------
        expression_matrix : np.ndarray
            Expression matrix (cells x genes)
        gene_names : List[str]
            List of gene names
        genes : List[str]
            Genes to include in the plot
        groups : np.ndarray
            Group assignments for each cell
        unique_groups : np.ndarray
            Unique group names

        Returns
        -------
        pd.DataFrame
            DataFrame with columns: gene, group, mean_expr, pct_expr
        """
        # Get gene indices
        gene_indices = [gene_names.index(gene) for gene in genes if gene in gene_names]

        # Initialize results
        results = []

        for i, gene_idx in enumerate(gene_indices):
            gene = genes[i]
            gene_expr = expression_matrix[:, gene_idx]

            for group in unique_groups:
                # Get cells in this group
                group_mask = groups == group
                group_expr = gene_expr[group_mask]

                # Calculate statistics
                mean_expr = np.mean(group_expr)
                pct_expr = np.mean(group_expr > 0)  # Percentage of cells expressing

                results.append(
                    {"gene": gene, "group": group, "mean_expr": mean_expr, "pct_expr": pct_expr}
                )

        return pd.DataFrame(results)

    def _create_dot_plot_figure(
        self,
        dot_data: pd.DataFrame,
        genes: List[str],
        groups: np.ndarray,
        dot_size_range: List[float],
        color_scale: str,
        show_legend: bool,
        show_axes: bool,
        sort_genes: bool,
        sort_groups: bool,
    ) -> go.Figure:
        """
        Create the actual dot plot figure.

        Parameters
        ----------
        dot_data : pd.DataFrame
            Data for the dot plot
        genes : List[str]
            List of genes
        groups : np.ndarray
            List of groups
        dot_size_range : List[float]
            Size range for dots
        color_scale : str
            Color scale name
        show_legend : bool
            Whether to show legend
        show_axes : bool
            Whether to show axes
        sort_genes : bool
            Whether to sort genes
        sort_groups : bool
            Whether to sort groups

        Returns
        -------
        go.Figure
            The dot plot figure
        """
        # Sort data if requested
        if sort_genes:
            genes = sorted(genes)
        if sort_groups:
            groups = np.array(sorted(groups))

        # Create figure
        fig = go.Figure()

        # Prepare data for single scatter plot
        x_coords = []
        y_coords = []
        sizes = []
        colors = []
        hover_texts = []

        for _, row in dot_data.iterrows():
            gene = row["gene"]
            group = row["group"]
            mean_expr = row["mean_expr"]
            pct_expr = row["pct_expr"]

            # Calculate dot size based on percentage
            dot_size = dot_size_range[0] + (dot_size_range[1] - dot_size_range[0]) * pct_expr

            # Collect data
            x_coords.append(group)
            y_coords.append(gene)
            sizes.append(dot_size)
            colors.append(mean_expr)
            hover_texts.append(
                f"<b>{gene}</b><br>"
                + f"Group: {group}<br>"
                + f"Mean Expression: {mean_expr:.3f}<br>"
                + f"Pct Expressing: {pct_expr:.1%}<br>"
                + "<extra></extra>"
            )

        # Add single scatter plot with all dots
        fig.add_trace(
            go.Scatter(
                x=x_coords,
                y=y_coords,
                mode="markers",
                marker=dict(
                    size=sizes,
                    color=colors,
                    colorscale=color_scale,
                    showscale=False,  # We'll handle colorbar in layout
                    coloraxis="coloraxis",  # Use the coloraxis from layout
                    line=dict(width=1, color="white"),
                    sizemode="diameter",
                    sizemin=dot_size_range[0],
                ),
                name="Gene Expression",
                showlegend=False,
                hovertemplate="%{customdata}<extra></extra>",
                customdata=hover_texts,
            )
        )

        # Update layout
        fig.update_layout(
            title=self.title,
            xaxis_title="Group",
            yaxis_title="Gene",
            showlegend=show_legend,
            hovermode="closest",
            width=900,  # Increased width for better spacing
            height=600,
            margin=dict(l=80, r=120, t=80, b=80),  # Added margins for color bar space
            coloraxis=dict(
                colorbar=dict(
                    len=0.8,  # Make color bar shorter
                    thickness=20,  # Make color bar thicker
                    x=1.02,  # Position color bar further right
                    xpad=20,  # Add padding between plot and color bar
                    title=dict(text="Mean Expression", side="right", font=dict(size=12)),
                    tickfont=dict(size=10),  # Smaller font for tick labels
                )
            ),
        )

        # Update axes
        if show_axes:
            fig.update_xaxes(categoryorder="array", categoryarray=groups, tickangle=45)
            fig.update_yaxes(categoryorder="array", categoryarray=genes)

        return fig

    def _create_error_figure(self, message: str) -> go.Figure:
        """Create an error figure with the given message."""
        fig = go.Figure()
        fig.add_annotation(
            text=message,
            xref="paper",
            yref="paper",
            x=0.5,
            y=0.5,
            showarrow=False,
            font=dict(size=16, color="red"),
        )
        fig.update_layout(
            title=self.title,
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            width=800,
            height=600,
        )
        return fig

    def get_selection_code(self) -> str:
        """
        Generate Python code for the current selection.

        Returns
        -------
        str
            Python code that reproduces the current selection
        """
        if self._selection is None:
            return "# No selection made in dot plot panel"

        # Get configuration
        genes = self.get_config("genes")
        group_by = self.get_config("group_by")

        code_lines = [
            "# Dot plot panel selection code",
            f"# Selected genes: {genes}",
            f"# Grouped by: {group_by}",
            "",
            "# Get selected cells",
            "selected_cells = adata.obs_names[selected_mask]",
            "",
            "# Filter data to selected cells",
            "selected_adata = adata[selected_cells]",
            "",
            "# Create dot plot for selected cells",
            "import pysee",
            f"dot_panel = pysee.DotPlotPanel('dotplot', genes={genes}, group_by='{group_by}')",
            "dot_panel.data_wrapper = pysee.AnnDataWrapper(selected_adata)",
            "fig = dot_panel.render()",
            "fig.show()",
        ]

        return "\n".join(code_lines)

    def _on_data_changed(self) -> None:
        """Called when the data wrapper changes."""
        # Reset any cached data
        pass

    def _on_selection_changed(self) -> None:
        """Called when the selection changes."""
        # Update visualization if needed
        pass

    def _on_config_changed(self) -> None:
        """Called when the configuration changes."""
        # Update visualization if needed
        pass
