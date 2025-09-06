"""
UMAP embedding visualization panel for PySEE.

This module provides the UMAPPanel class for visualizing dimensionality reduction
embeddings like UMAP, t-SNE, and PCA.
"""

from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ..core.data import AnnDataWrapper
from .base import BasePanel
from .plotly_selection import PlotlySelectionTools
from .advanced_selection import AdvancedSelectionManager


class UMAPPanel(BasePanel):
    """
    Panel for visualizing UMAP and other dimensionality reduction embeddings.

    This panel creates interactive scatter plots of embeddings with support for
    color mapping, point selection, and linked interactions.

    Parameters
    ----------
    panel_id : str
        Unique identifier for the panel
    embedding : str, default 'X_umap'
        Key for the embedding in adata.obsm
    color : str, optional
        Column name in adata.obs to use for coloring points
    title : str, optional
        Display title for the panel
    """

    def __init__(
        self,
        panel_id: str,
        embedding: str = "X_umap",
        color: Optional[str] = None,
        title: Optional[str] = None,
    ):
        super().__init__(panel_id, title)

        self.set_config("embedding", embedding)
        self.set_config("color", color)
        self.set_config("point_size", 3)
        self.set_config("opacity", 0.7)
        self.set_config("show_legend", True)
        self.set_config("width", 600)
        self.set_config("height", 400)

    def _check_data_requirements(self) -> bool:
        """Check if the data wrapper has the required embedding."""
        if self._data_wrapper is None:
            return False

        embedding_key = self.get_config("embedding")
        return embedding_key in self._data_wrapper.get_obsm_keys()

    def render(self) -> go.Figure:
        """
        Render the UMAP visualization.

        Returns
        -------
        plotly.graph_objects.Figure
            Interactive scatter plot of the embedding
        """
        if not self.validate_data():
            raise ValueError("Panel data requirements not met")

        # Get embedding data
        embedding_key = self.get_config("embedding")
        coords = self._data_wrapper.get_embedding_data(embedding_key)

        if coords.shape[1] < 2:
            raise ValueError(f"Embedding '{embedding_key}' must have at least 2 dimensions")

        # Prepare data for plotting
        plot_data = pd.DataFrame(
            {
                "x": coords[:, 0],
                "y": coords[:, 1],
            }
        )

        # Add color information if specified
        color_col = self.get_config("color")
        if color_col and color_col in self._data_wrapper.get_obs_columns():
            plot_data["color"] = self._data_wrapper.adata.obs[color_col].values

        # Add selection information
        if self._selection is not None:
            plot_data["selected"] = self._selection
        else:
            plot_data["selected"] = False

        # Create the plot
        fig = go.Figure()

        if color_col and color_col in self._data_wrapper.get_obs_columns():
            # Color by specified column
            color_values = plot_data["color"]

            # Check if column is categorical or should be treated as categorical
            is_categorical = (
                self._data_wrapper.adata.obs[color_col].dtype == "category"
                or self._data_wrapper.adata.obs[color_col].dtype == "object"
            )

            if is_categorical:
                # Categorical coloring
                if self._data_wrapper.adata.obs[color_col].dtype == "category":
                    categories = self._data_wrapper.adata.obs[color_col].cat.categories
                else:
                    categories = self._data_wrapper.adata.obs[color_col].unique()
                colors = px.colors.qualitative.Set1[: len(categories)]

                for i, category in enumerate(categories):
                    mask = color_values == category
                    fig.add_trace(
                        go.Scatter(
                            x=plot_data.loc[mask, "x"],
                            y=plot_data.loc[mask, "y"],
                            mode="markers",
                            marker=dict(
                                size=self.get_config("point_size"),
                                color=colors[i % len(colors)],
                                opacity=self.get_config("opacity"),
                            ),
                            name=str(category),
                            showlegend=self.get_config("show_legend"),
                            hovertemplate=f"<b>{category}</b><br>"
                            + "X: %{x:.2f}<br>"
                            + "Y: %{y:.2f}<br>"
                            + "<extra></extra>",
                        )
                    )
            else:
                # Continuous coloring
                fig.add_trace(
                    go.Scatter(
                        x=plot_data["x"],
                        y=plot_data["y"],
                        mode="markers",
                        marker=dict(
                            size=self.get_config("point_size"),
                            color=color_values,
                            colorscale="Viridis",
                            opacity=self.get_config("opacity"),
                            colorbar=dict(title=color_col),
                        ),
                        name=color_col,
                        showlegend=False,
                        hovertemplate=f"<b>{color_col}</b>: %{{marker.color:.2f}}<br>"
                        + "X: %{x:.2f}<br>"
                        + "Y: %{y:.2f}<br>"
                        + "<extra></extra>",
                    )
                )
        else:
            # No coloring, single color
            fig.add_trace(
                go.Scatter(
                    x=plot_data["x"],
                    y=plot_data["y"],
                    mode="markers",
                    marker=dict(
                        size=self.get_config("point_size"),
                        color="blue",
                        opacity=self.get_config("opacity"),
                    ),
                    name="Cells",
                    showlegend=False,
                    hovertemplate="X: %{x:.2f}<br>" + "Y: %{y:.2f}<br>" + "<extra></extra>",
                )
            )

        # Highlight selected points
        if self._selection is not None and np.any(self._selection):
            selected_data = plot_data[self._selection]
            fig.add_trace(
                go.Scatter(
                    x=selected_data["x"],
                    y=selected_data["y"],
                    mode="markers",
                    marker=dict(
                        size=self.get_config("point_size") + 2,
                        color="red",
                        opacity=1.0,
                        symbol="diamond",
                    ),
                    name="Selected",
                    showlegend=True,
                    hovertemplate="<b>Selected</b><br>"
                    + "X: %{x:.2f}<br>"
                    + "Y: %{y:.2f}<br>"
                    + "<extra></extra>",
                )
            )

        # Update layout
        fig.update_layout(
            title=self.title,
            xaxis_title=f"{embedding_key} 1",
            yaxis_title=f"{embedding_key} 2",
            width=self.get_config("width"),
            height=self.get_config("height"),
            showlegend=self.get_config("show_legend"),
            hovermode="closest",
        )

        # Add selection functionality
        fig.update_traces(
            selector=dict(type="scatter"),
            selected=dict(marker=dict(color="red", size=8)),
            unselected=dict(marker=dict(opacity=0.3)),
        )

        # Add advanced selection tools if enabled
        if self._advanced_selection_enabled and self._selection_manager is not None:
            # Add selection UI elements
            PlotlySelectionTools.add_selection_buttons(fig)
            PlotlySelectionTools.add_selection_mode_buttons(fig)
            PlotlySelectionTools.add_undo_redo_buttons(fig)
            PlotlySelectionTools.add_interactive_descriptions(fig)
            
            # Configure selection events
            PlotlySelectionTools.configure_selection_events(fig, self._on_plotly_selection)

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
            return "# No selection made"

        # Count selected cells
        n_selected = np.sum(self._selection)
        n_total = len(self._selection)

        code_lines = [
            f"# UMAP Panel Selection: {n_selected}/{n_total} cells selected",
            f"# Panel ID: {self.panel_id}",
            f"# Embedding: {self.get_config('embedding')}",
        ]

        if self.get_config("color"):
            code_lines.append(f"# Color by: {self.get_config('color')}")

        code_lines.extend(
            [
                "",
                "# Selection mask",
                f"selection_mask = np.array({self._selection.tolist()})",
                "",
                "# Get selected cells",
                "selected_cells = adata[selection_mask]",
                "",
                "# Get selected cell indices",
                "selected_indices = np.where(selection_mask)[0]",
                "",
                "# Get selected cell names (if available)",
                "if hasattr(adata, 'obs_names'):",
                "    selected_cell_names = adata.obs_names[selection_mask]",
            ]
        )

        return "\n".join(code_lines)

    def _on_data_changed(self) -> None:
        """Called when the data wrapper changes."""
        # Reset selection when data changes
        self._selection = None

    def _on_selection_changed(self) -> None:
        """Called when the selection changes."""
        # This could be used to notify linked panels
        pass

    def _on_config_changed(self) -> None:
        """Called when the configuration changes."""
        # This could be used to update the visualization
        pass

    def _on_plotly_selection(self, selection_data: dict) -> None:
        """Handle Plotly selection events for advanced selection tools."""
        if not self._advanced_selection_enabled or self._selection_manager is None:
            return

        try:
            # Get coordinates for selection
            coords = self._data_wrapper.get_embedding_data(self.get_config("embedding"))
            
            # Handle different selection types
            selection_type = selection_data.get("type", "point")
            
            if selection_type == "point":
                # Point selection
                point_indices = selection_data.get("points", [])
                indices = [p.get("pointIndex", 0) for p in point_indices if "pointIndex" in p]
                if indices:
                    self._selection_manager.create_point_selection(indices)
                    
            elif selection_type == "rectangular":
                # Rectangular selection
                bounds = selection_data.get("range", {})
                if "x" in bounds and "y" in bounds:
                    x_range = bounds["x"]
                    y_range = bounds["y"]
                    bounds_tuple = (x_range[0], y_range[0], x_range[1], y_range[1])
                    self._selection_manager.create_rectangular_selection(bounds_tuple, coords)
                    
            elif selection_type == "lasso":
                # Lasso selection
                path = selection_data.get("path", [])
                if path:
                    path_coords = [(p.get("x", 0), p.get("y", 0)) for p in path]
                    self._selection_manager.create_lasso_selection(path_coords, coords)
                    
            # Update the panel's selection
            if self._selection_manager.get_current_selection() is not None:
                self._selection = self._selection_manager.get_current_selection()
                self._on_selection_changed()
                
        except Exception as e:
            print(f"Error handling selection: {e}")

    def enable_advanced_selection(self) -> None:
        """Enable advanced selection tools for this panel."""
        self._advanced_selection_enabled = True
        if self._data_wrapper is not None:
            self._selection_manager = AdvancedSelectionManager(self._data_wrapper)

    def disable_advanced_selection(self) -> None:
        """Disable advanced selection tools for this panel."""
        self._advanced_selection_enabled = False
        self._selection_manager = None

    def set_embedding(self, embedding: str) -> None:
        """
        Set the embedding to visualize.

        Parameters
        ----------
        embedding : str
            Key for the embedding in adata.obsm
        """
        self.set_config("embedding", embedding)

    def set_color(self, color: Optional[str]) -> None:
        """
        Set the column to use for coloring points.

        Parameters
        ----------
        color : str or None
            Column name in adata.obs to use for coloring, or None for no coloring
        """
        self.set_config("color", color)

    def set_point_size(self, size: float) -> None:
        """
        Set the size of the points in the plot.

        Parameters
        ----------
        size : float
            Point size
        """
        self.set_config("point_size", size)

    def set_opacity(self, opacity: float) -> None:
        """
        Set the opacity of the points.

        Parameters
        ----------
        opacity : float
            Opacity value between 0 and 1
        """
        if not 0 <= opacity <= 1:
            raise ValueError("Opacity must be between 0 and 1")
        self.set_config("opacity", opacity)
