"""
Violin plot panel for gene expression visualization in PySEE.

This module provides the ViolinPanel class for visualizing gene expression
distributions with support for grouping and statistical comparisons.
"""

from typing import Any, Dict, List, Optional, Union
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from .base import BasePanel
from ..core.data import AnnDataWrapper


class ViolinPanel(BasePanel):
    """
    Panel for visualizing gene expression distributions using violin plots.
    
    This panel creates interactive violin plots for gene expression data with
    support for grouping, statistical comparisons, and linked selections.
    
    Parameters
    ----------
    panel_id : str
        Unique identifier for the panel
    gene : str, optional
        Gene name to visualize
    group_by : str, optional
        Column name in adata.obs to use for grouping
    title : str, optional
        Display title for the panel
    """
    
    def __init__(
        self,
        panel_id: str,
        gene: Optional[str] = None,
        group_by: Optional[str] = None,
        title: Optional[str] = None
    ):
        super().__init__(panel_id, title)
        
        self.set_config('gene', gene)
        self.set_config('group_by', group_by)
        self.set_config('plot_type', 'violin')  # 'violin', 'box', 'strip'
        self.set_config('show_points', True)
        self.set_config('show_box', True)
        self.set_config('width', 600)
        self.set_config('height', 400)
        self.set_config('color_by_group', True)
    
    def _check_data_requirements(self) -> bool:
        """Check if the data wrapper has the required data."""
        if self._data_wrapper is None:
            return False
        
        # Check if gene is specified and exists
        gene = self.get_config('gene')
        if gene and gene not in self._data_wrapper.adata.var_names:
            return False
        
        # Check if group_by column exists
        group_by = self.get_config('group_by')
        if group_by and group_by not in self._data_wrapper.get_obs_columns():
            return False
        
        return True
    
    def render(self) -> go.Figure:
        """
        Render the violin plot visualization.
        
        Returns
        -------
        plotly.graph_objects.Figure
            Interactive violin plot
        """
        if not self.validate_data():
            raise ValueError("Panel data requirements not met")
        
        gene = self.get_config('gene')
        if not gene:
            raise ValueError("No gene specified for visualization")
        
        # Get expression data
        expression_data = self._data_wrapper.get_expression_data(gene)
        if expression_data.ndim > 1:
            expression_data = expression_data.flatten()
        
        # Prepare data for plotting
        plot_data = pd.DataFrame({
            'expression': expression_data,
        })
        
        # Add grouping information
        group_by = self.get_config('group_by')
        if group_by and group_by in self._data_wrapper.get_obs_columns():
            plot_data['group'] = self._data_wrapper.adata.obs[group_by].values
        else:
            plot_data['group'] = 'All'
        
        # Add selection information
        if self._selection is not None:
            plot_data['selected'] = self._selection
        else:
            plot_data['selected'] = False
        
        # Create the plot
        fig = go.Figure()
        
        plot_type = self.get_config('plot_type')
        show_points = self.get_config('show_points')
        show_box = self.get_config('show_box')
        
        if group_by and group_by in self._data_wrapper.get_obs_columns():
            # Grouped plot
            groups = plot_data['group'].unique()
            colors = px.colors.qualitative.Set1[:len(groups)]
            
            for i, group in enumerate(groups):
                group_data = plot_data[plot_data['group'] == group]
                
                if plot_type == 'violin':
                    fig.add_trace(go.Violin(
                        y=group_data['expression'],
                        name=str(group),
                        box_visible=show_box,
                        meanline_visible=True,
                        fillcolor=colors[i % len(colors)],
                        opacity=0.7,
                        line_color=colors[i % len(colors)],
                        hovertemplate=f'<b>{group}</b><br>' +
                                    f'{gene}: %{{y:.2f}}<br>' +
                                    '<extra></extra>',
                    ))
                elif plot_type == 'box':
                    fig.add_trace(go.Box(
                        y=group_data['expression'],
                        name=str(group),
                        fillcolor=colors[i % len(colors)],
                        opacity=0.7,
                        line_color=colors[i % len(colors)],
                        hovertemplate=f'<b>{group}</b><br>' +
                                    f'{gene}: %{{y:.2f}}<br>' +
                                    '<extra></extra>',
                    ))
                elif plot_type == 'strip':
                    fig.add_trace(go.Box(
                        y=group_data['expression'],
                        name=str(group),
                        boxpoints='all',
                        jitter=0.3,
                        pointpos=-1.8,
                        fillcolor=colors[i % len(colors)],
                        opacity=0.7,
                        line_color=colors[i % len(colors)],
                        hovertemplate=f'<b>{group}</b><br>' +
                                    f'{gene}: %{{y:.2f}}<br>' +
                                    '<extra></extra>',
                    ))
                
                # Add points if requested
                if show_points and plot_type != 'strip':
                    fig.add_trace(go.Scatter(
                        y=group_data['expression'],
                        mode='markers',
                        marker=dict(
                            color=colors[i % len(colors)],
                            size=4,
                            opacity=0.6,
                        ),
                        name=f'{group} points',
                        showlegend=False,
                        hovertemplate=f'<b>{group}</b><br>' +
                                    f'{gene}: %{{y:.2f}}<br>' +
                                    '<extra></extra>',
                    ))
        else:
            # Single group plot
            if plot_type == 'violin':
                fig.add_trace(go.Violin(
                    y=plot_data['expression'],
                    name='All',
                    box_visible=show_box,
                    meanline_visible=True,
                    fillcolor='lightblue',
                    opacity=0.7,
                    line_color='blue',
                    hovertemplate=f'{gene}: %{{y:.2f}}<br>' +
                                '<extra></extra>',
                ))
            elif plot_type == 'box':
                fig.add_trace(go.Box(
                    y=plot_data['expression'],
                    name='All',
                    fillcolor='lightblue',
                    opacity=0.7,
                    line_color='blue',
                    hovertemplate=f'{gene}: %{{y:.2f}}<br>' +
                                '<extra></extra>',
                ))
            elif plot_type == 'strip':
                fig.add_trace(go.Box(
                    y=plot_data['expression'],
                    name='All',
                    boxpoints='all',
                    jitter=0.3,
                    pointpos=-1.8,
                    fillcolor='lightblue',
                    opacity=0.7,
                    line_color='blue',
                    hovertemplate=f'{gene}: %{{y:.2f}}<br>' +
                                '<extra></extra>',
                ))
            
            # Add points if requested
            if show_points and plot_type != 'strip':
                fig.add_trace(go.Scatter(
                    y=plot_data['expression'],
                    mode='markers',
                    marker=dict(
                        color='blue',
                        size=4,
                        opacity=0.6,
                    ),
                    name='Points',
                    showlegend=False,
                    hovertemplate=f'{gene}: %{{y:.2f}}<br>' +
                                '<extra></extra>',
                ))
        
        # Highlight selected points
        if self._selection is not None and np.any(self._selection):
            selected_data = plot_data[self._selection]
            fig.add_trace(go.Scatter(
                y=selected_data['expression'],
                mode='markers',
                marker=dict(
                    color='red',
                    size=8,
                    symbol='diamond',
                    opacity=1.0,
                ),
                name='Selected',
                showlegend=True,
                hovertemplate=f'<b>Selected</b><br>' +
                            f'{gene}: %{{y:.2f}}<br>' +
                            '<extra></extra>',
            ))
        
        # Update layout
        fig.update_layout(
            title=f"{self.title} - {gene}",
            yaxis_title=f'{gene} Expression',
            xaxis_title=group_by if group_by else 'Group',
            width=self.get_config('width'),
            height=self.get_config('height'),
            showlegend=True,
            hovermode='closest',
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
            return "# No selection made"
        
        gene = self.get_config('gene')
        group_by = self.get_config('group_by')
        
        # Count selected cells
        n_selected = np.sum(self._selection)
        n_total = len(self._selection)
        
        code_lines = [
            f"# Violin Panel Selection: {n_selected}/{n_total} cells selected",
            f"# Panel ID: {self.panel_id}",
            f"# Gene: {gene}",
        ]
        
        if group_by:
            code_lines.append(f"# Group by: {group_by}")
        
        code_lines.extend([
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
            "# Get expression values for selected cells",
            f"selected_expression = adata[selection_mask, '{gene}'].X",
            "",
            "# Get selected cell names (if available)",
            "if hasattr(adata, 'obs_names'):",
            "    selected_cell_names = adata.obs_names[selection_mask]",
        ])
        
        if group_by:
            code_lines.extend([
                "",
                "# Get group information for selected cells",
                f"selected_groups = adata.obs.loc[selection_mask, '{group_by}']",
            ])
        
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
    
    def set_gene(self, gene: str) -> None:
        """
        Set the gene to visualize.
        
        Parameters
        ----------
        gene : str
            Gene name to visualize
        """
        self.set_config('gene', gene)
    
    def set_group_by(self, group_by: Optional[str]) -> None:
        """
        Set the column to use for grouping.
        
        Parameters
        ----------
        group_by : str or None
            Column name in adata.obs to use for grouping, or None for no grouping
        """
        self.set_config('group_by', group_by)
    
    def set_plot_type(self, plot_type: str) -> None:
        """
        Set the type of plot to render.
        
        Parameters
        ----------
        plot_type : str
            Type of plot: 'violin', 'box', or 'strip'
        """
        if plot_type not in ['violin', 'box', 'strip']:
            raise ValueError("plot_type must be 'violin', 'box', or 'strip'")
        self.set_config('plot_type', plot_type)
    
    def set_show_points(self, show_points: bool) -> None:
        """
        Set whether to show individual points.
        
        Parameters
        ----------
        show_points : bool
            Whether to show individual points
        """
        self.set_config('show_points', show_points)
    
    def set_show_box(self, show_box: bool) -> None:
        """
        Set whether to show box plot inside violin plot.
        
        Parameters
        ----------
        show_box : bool
            Whether to show box plot inside violin plot
        """
        self.set_config('show_box', show_box)
