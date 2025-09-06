"""
Plotly integration for advanced selection tools.

This module provides utilities for integrating advanced selection tools
with Plotly figures, including selection buttons, event handling, and
visual feedback.
"""

from typing import Dict, Any, List, Tuple, Callable, Optional
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

from .selection_types import SelectionType, SelectionMode


class PlotlySelectionTools:
    """Utilities for adding advanced selection tools to Plotly figures."""
    
    @staticmethod
    def add_selection_buttons(fig: go.Figure, 
                            position: str = "top",
                            show_labels: bool = True) -> go.Figure:
        """
        Add selection tool buttons to plotly figure.
        
        Parameters
        ----------
        fig : go.Figure
            Plotly figure to modify
        position : str
            Button position ("top", "bottom", "left", "right")
        show_labels : bool
            Whether to show button labels
            
        Returns
        -------
        go.Figure
            Modified figure with selection buttons
        """
        # Define selection tool buttons with descriptions
        buttons = [
            {
                "label": "👆 Point Select" if show_labels else "👆",
                "method": "restyle",
                "args": [{"dragmode": "select"}],
                "args2": [{"dragmode": "pan"}],
                "name": "point_select",
                "visible": True
            },
            {
                "label": "▭ Box Select" if show_labels else "▭",
                "method": "restyle", 
                "args": [{"dragmode": "select"}],
                "name": "box_select",
                "visible": True
            },
            {
                "label": "⚬ Lasso Select" if show_labels else "⚬",
                "method": "restyle",
                "args": [{"dragmode": "lasso"}],
                "name": "lasso_select", 
                "visible": True
            },
            {
                "label": "◇ Pan" if show_labels else "◇",
                "method": "restyle",
                "args": [{"dragmode": "pan"}],
                "name": "pan_tool",
                "visible": True
            },
            {
                "label": "🔍 Zoom" if show_labels else "🔍",
                "method": "restyle",
                "args": [{"dragmode": "zoom"}],
                "name": "zoom_tool",
                "visible": True
            }
        ]
        
        # Add buttons to figure
        fig.update_layout(
            updatemenus=[
                {
                    "type": "buttons",
                    "direction": "left" if position in ["top", "bottom"] else "down",
                    "x": 0.0 if position == "left" else (1.0 if position == "right" else 0.5),
                    "y": 1.02 if position == "top" else (-0.1 if position == "bottom" else 0.5),
                    "xanchor": "left" if position == "left" else ("right" if position == "right" else "center"),
                    "yanchor": "bottom" if position in ["top", "bottom"] else "middle",
                    "buttons": buttons,
                    "bgcolor": "rgba(255,255,255,0.8)",
                    "bordercolor": "rgba(0,0,0,0.2)",
                    "borderwidth": 1,
                    "font": {"size": 10}
                }
            ]
        )
        
        return fig
    
    @staticmethod
    def add_selection_mode_buttons(fig: go.Figure,
                                 position: str = "topright",
                                 show_labels: bool = True) -> go.Figure:
        """
        Add selection mode buttons (replace, add, subtract, intersect).
        
        Parameters
        ----------
        fig : go.Figure
            Plotly figure to modify
        position : str
            Button position
        show_labels : bool
            Whether to show button labels
            
        Returns
        -------
        go.Figure
            Modified figure with mode buttons
        """
        # Define selection mode buttons with descriptions
        mode_buttons = [
            {
                "label": "🔄 Replace" if show_labels else "🔄",
                "method": "restyle",
                "args": [{"selectedpoints": None}],  # Clear selection
                "args2": [{"selectedpoints": None}],
                "name": "replace_mode"
            },
            {
                "label": "➕ Add" if show_labels else "➕",
                "method": "restyle",
                "args": [{}],  # Placeholder - actual logic handled in callbacks
                "name": "add_mode"
            },
            {
                "label": "➖ Subtract" if show_labels else "➖", 
                "method": "restyle",
                "args": [{}],  # Placeholder
                "name": "subtract_mode"
            },
            {
                "label": "∩ Intersect" if show_labels else "∩",
                "method": "restyle", 
                "args": [{}],  # Placeholder
                "name": "intersect_mode"
            }
        ]
        
        # Get existing updatemenus or create new list
        existing_menus = list(fig.layout.updatemenus or [])
        
        # Add mode buttons
        new_menu = {
            "type": "buttons",
            "direction": "down",
            "x": 1.02,
            "y": 1.0,
            "xanchor": "left",
            "yanchor": "top", 
            "buttons": mode_buttons,
            "bgcolor": "rgba(240,240,240,0.8)",
            "bordercolor": "rgba(0,0,0,0.2)",
            "borderwidth": 1,
            "font": {"size": 9}
        }
        
        fig.update_layout(updatemenus=existing_menus + [new_menu])
        
        return fig
    
    @staticmethod
    def add_undo_redo_buttons(fig: go.Figure,
                            position: str = "bottomright") -> go.Figure:
        """
        Add undo/redo buttons.
        
        Parameters
        ----------
        fig : go.Figure
            Plotly figure to modify
        position : str
            Button position
            
        Returns
        -------
        go.Figure
            Modified figure with undo/redo buttons
        """
        # Define undo/redo buttons
        history_buttons = [
            {
                "label": "↶ Undo",
                "method": "restyle",
                "args": [{}]  # Placeholder - handled in callbacks
            },
            {
                "label": "↷ Redo",
                "method": "restyle",
                "args": [{}]  # Placeholder
            },
            {
                "label": "🗑️ Clear",
                "method": "restyle",
                "args": [{"selectedpoints": None}]
            }
        ]
        
        # Get existing updatemenus
        existing_menus = list(fig.layout.updatemenus or [])
        
        # Add history buttons
        new_menu = {
            "type": "buttons", 
            "direction": "left",
            "x": 1.0,
            "y": -0.1,
            "xanchor": "right",
            "yanchor": "top",
            "buttons": history_buttons,
            "bgcolor": "rgba(255,255,255,0.8)",
            "bordercolor": "rgba(0,0,0,0.2)",
            "borderwidth": 1,
            "font": {"size": 9}
        }
        
        fig.update_layout(updatemenus=existing_menus + [new_menu])
        
        return fig
    
    @staticmethod
    def configure_selection_events(fig: go.Figure,
                                 selection_callback: Optional[Callable] = None) -> go.Figure:
        """
        Configure plotly selection events and callbacks.
        
        Parameters
        ----------
        fig : go.Figure
            Plotly figure to modify
        selection_callback : Callable, optional
            Function to call when selection changes
            
        Returns
        -------
        go.Figure
            Modified figure with selection events configured
        """
        # Enable selection events
        fig.update_traces(
            selected=dict(marker=dict(opacity=1.0, size=8)),
            unselected=dict(marker=dict(opacity=0.3, size=4)),
            selectedpoints=[]
        )
        
        # Configure layout for selection
        fig.update_layout(
            clickmode='event+select',
            selectdirection='diagonal',
            dragmode='select'
        )
        
        return fig
    
    @staticmethod
    def add_selection_info_box(fig: go.Figure,
                             selection_stats: Dict[str, Any],
                             position: str = "topleft") -> go.Figure:
        """
        Add information box showing selection statistics.
        
        Parameters
        ----------
        fig : go.Figure
            Plotly figure to modify
        selection_stats : Dict[str, Any]
            Selection statistics to display
        position : str
            Position of info box
            
        Returns
        -------
        go.Figure
            Modified figure with info box
        """
        if not selection_stats or selection_stats.get("n_selected", 0) == 0:
            info_text = "No selection"
        else:
            n_selected = selection_stats.get("n_selected", 0)
            n_total = selection_stats.get("n_total", 0)
            percentage = selection_stats.get("selection_percentage", 0)
            
            info_text = f"Selected: {n_selected:,} / {n_total:,} ({percentage:.1f}%)"
        
        # Position mapping
        pos_map = {
            "topleft": {"x": 0.02, "y": 0.98, "xanchor": "left", "yanchor": "top"},
            "topright": {"x": 0.98, "y": 0.98, "xanchor": "right", "yanchor": "top"},
            "bottomleft": {"x": 0.02, "y": 0.02, "xanchor": "left", "yanchor": "bottom"},
            "bottomright": {"x": 0.98, "y": 0.02, "xanchor": "right", "yanchor": "bottom"}
        }
        
        pos = pos_map.get(position, pos_map["topleft"])
        
        # Add annotation for info box
        fig.add_annotation(
            text=info_text,
            x=pos["x"],
            y=pos["y"],
            xref="paper",
            yref="paper",
            xanchor=pos["xanchor"],
            yanchor=pos["yanchor"],
            bgcolor="rgba(255,255,255,0.9)",
            bordercolor="rgba(0,0,0,0.3)",
            borderwidth=1,
            font=dict(size=10, color="black")
        )
        
        return fig
    
    @staticmethod
    def highlight_selection(fig: go.Figure,
                          selection_mask: np.ndarray,
                          trace_index: int = 0) -> go.Figure:
        """
        Highlight selected points in the figure.
        
        Parameters
        ----------
        fig : go.Figure
            Plotly figure to modify
        selection_mask : np.ndarray
            Boolean mask of selected points
        trace_index : int
            Index of trace to update
            
        Returns
        -------
        go.Figure
            Modified figure with highlighted selection
        """
        if trace_index >= len(fig.data):
            return fig
        
        # Get selected indices
        selected_indices = np.where(selection_mask)[0].tolist()
        
        # Update trace with selection
        fig.data[trace_index].selectedpoints = selected_indices
        
        return fig
    
    @staticmethod
    def create_selection_overlay(coords: np.ndarray,
                               selection_path: List[Tuple[float, float]],
                               selection_type: SelectionType) -> go.Scatter:
        """
        Create overlay showing selection area.
        
        Parameters
        ----------
        coords : np.ndarray
            Coordinate array
        selection_path : List[Tuple[float, float]]
            Path defining selection area
        selection_type : SelectionType
            Type of selection
            
        Returns
        -------
        go.Scatter
            Scatter trace for selection overlay
        """
        if not selection_path:
            return go.Scatter(x=[], y=[], mode='lines', showlegend=False)
        
        x_coords, y_coords = zip(*selection_path)
        
        # Close the path for polygon/lasso
        if selection_type in [SelectionType.LASSO, SelectionType.POLYGON]:
            x_coords = list(x_coords) + [x_coords[0]]
            y_coords = list(y_coords) + [y_coords[0]]
        
        return go.Scatter(
            x=x_coords,
            y=y_coords,
            mode='lines',
            line=dict(color='red', width=2, dash='dash'),
            name=f'{selection_type.value.title()} Selection',
            showlegend=False,
            hoverinfo='skip'
        )
    
    @staticmethod
    def add_interactive_descriptions(fig: go.Figure, 
                                   title: str = "",
                                   descriptions: Optional[Dict[str, str]] = None) -> go.Figure:
        """
        Add interactive descriptions and tooltips to the figure.
        
        Parameters
        ----------
        fig : go.Figure
            Plotly figure to modify
        title : str
            Main title for the chart
        descriptions : Dict[str, str], optional
            Dictionary of descriptions for various elements
            
        Returns
        -------
        go.Figure
            Modified figure with descriptions
        """
        if descriptions is None:
            descriptions = {
                "selection_tools": "Use selection tools to interact with data points",
                "selection_modes": "Choose how new selections combine with existing ones",
                "navigation": "Pan and zoom to explore the data"
            }
        
        # Create comprehensive description text
        desc_text = f"""
        <b>{title}</b><br><br>
        
        <b>🔧 Selection Tools:</b><br>
        • <b>👆 Point Select</b>: Click individual points to select them<br>
        • <b>▭ Box Select</b>: Drag to draw a rectangle and select points inside<br>
        • <b>⚬ Lasso Select</b>: Draw a freeform shape to select points inside<br>
        • <b>◇ Pan</b>: Click and drag to move around the plot<br>
        • <b>🔍 Zoom</b>: Drag to zoom into a specific area<br><br>
        
        <b>🎯 Selection Modes:</b><br>
        • <b>🔄 Replace</b>: New selection replaces the current selection<br>
        • <b>➕ Add</b>: Add new points to the current selection (union)<br>
        • <b>➖ Subtract</b>: Remove points from the current selection<br>
        • <b>∩ Intersect</b>: Keep only points in both selections<br><br>
        
        <b>📚 History Controls:</b><br>
        • <b>↶ Undo</b>: Go back to previous selection state<br>
        • <b>↷ Redo</b>: Go forward to next selection state<br>
        • <b>🗑️ Clear</b>: Remove all selections<br><br>
        
        <b>💡 Tips:</b><br>
        • Hover over points to see details<br>
        • Use keyboard shortcuts: Ctrl+Z (undo), Ctrl+Y (redo)<br>
        • Double-click to reset zoom<br>
        • Hold Shift while selecting to add to selection
        """
        
        # Add detailed help (initially hidden) - positioned below the help button
        fig.add_annotation(
            text=desc_text,
            x=1.02,  # Same x as the help button
            y=0.95,  # Below the help button
            xref="paper",
            yref="paper",
            xanchor="left",
            yanchor="top",
            bgcolor="rgba(255,255,255,0.98)",
            bordercolor="rgba(0,0,0,0.4)",
            borderwidth=2,
            font=dict(size=10, color="black", family="Arial"),
            align="left",
            showarrow=False,
            width=400,
            visible=False,  # Initially hidden
            name="help_overlay"
        )
        
        # Get help annotation index
        help_index = len(fig.layout.annotations) - 1
        
        # Add single toggle button outside the chart area
        help_button = {
            "type": "buttons",
            "direction": "left",
            "x": 1.02,  # Outside the chart on the right
            "y": 1.0,
            "xanchor": "left",
            "yanchor": "top",
            "buttons": [
                {
                    "label": "❓ Show Help",
                    "method": "relayout",
                    "args": [{f"annotations[{help_index}].visible": True}],
                    "args2": [{f"annotations[{help_index}].visible": False}]
                }
            ],
            "bgcolor": "rgba(135,206,250,0.9)",
            "bordercolor": "rgba(0,0,0,0.3)",
            "borderwidth": 1,
            "font": {"size": 11, "color": "black"},
            "pad": {"r": 10, "t": 10}
        }
        
        # Get existing updatemenus
        existing_menus = list(fig.layout.updatemenus or [])
        fig.update_layout(updatemenus=existing_menus + [help_button])
        
        return fig
    
    @staticmethod
    def add_keyboard_shortcuts_info(fig: go.Figure) -> go.Figure:
        """
        Add keyboard shortcuts information to the figure.
        
        Parameters
        ----------
        fig : go.Figure
            Plotly figure to modify
            
        Returns
        -------
        go.Figure
            Modified figure with keyboard shortcuts info
        """
        shortcuts_text = """
        <b>⌨️ Keyboard Shortcuts:</b><br>
        • <b>Ctrl + Z</b>: Undo selection<br>
        • <b>Ctrl + Y</b>: Redo selection<br>
        • <b>Ctrl + A</b>: Select all points<br>
        • <b>Ctrl + D</b>: Clear selection<br>
        • <b>Shift + Click</b>: Add to selection<br>
        • <b>Double Click</b>: Reset zoom<br>
        • <b>Mouse Wheel</b>: Zoom in/out
        """
        
        fig.add_annotation(
            text=shortcuts_text,
            x=0.98,
            y=0.02,
            xref="paper",
            yref="paper",
            xanchor="right",
            yanchor="bottom",
            bgcolor="rgba(240,248,255,0.95)",
            bordercolor="rgba(0,0,0,0.2)",
            borderwidth=1,
            font=dict(size=9, color="black", family="Arial"),
            align="left",
            showarrow=False,
            visible=True
        )
        
        return fig
