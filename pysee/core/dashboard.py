"""
Main PySEE dashboard engine.

This module provides the PySEE class that manages panels, interactions,
and code export functionality.
"""

from typing import Any, Dict, List, Optional, Union
import numpy as np
import pandas as pd
from .data import AnnDataWrapper
from ..panels.base import BasePanel


class PySEE:
    """
    Main PySEE dashboard class for managing interactive visualizations.

    This class orchestrates multiple visualization panels, handles panel linking,
    manages selections, and provides code export functionality.

    Parameters
    ----------
    adata : anndata.AnnData
        The AnnData object to visualize
    title : str, optional
        Title for the dashboard
    """

    def __init__(self, adata: Any, title: Optional[str] = None) -> None:
        self._data_wrapper = AnnDataWrapper(adata)
        self._panels: Dict[str, BasePanel] = {}
        self._panel_order: List[str] = []
        self._title = title or "PySEE Dashboard"
        self._global_selection: Optional[np.ndarray] = None

    @property
    def data_wrapper(self) -> AnnDataWrapper:
        """Get the data wrapper."""
        return self._data_wrapper

    @property
    def panels(self) -> Dict[str, BasePanel]:
        """Get dictionary of all panels."""
        return self._panels.copy()

    @property
    def panel_order(self) -> List[str]:
        """Get ordered list of panel IDs."""
        return self._panel_order.copy()

    @property
    def title(self) -> str:
        """Get the dashboard title."""
        return self._title

    @title.setter
    def title(self, title: str) -> None:
        """Set the dashboard title."""
        self._title = title

    def add_panel(self, panel_id: str, panel: BasePanel) -> None:
        """
        Add a panel to the dashboard.

        Parameters
        ----------
        panel_id : str
            Unique identifier for the panel
        panel : BasePanel
            The panel to add
        """
        if panel_id in self._panels:
            raise ValueError(f"Panel with ID '{panel_id}' already exists")

        # Set the data wrapper for the panel
        panel.data_wrapper = self._data_wrapper

        # Add to dashboard
        self._panels[panel_id] = panel
        self._panel_order.append(panel_id)

    def remove_panel(self, panel_id: str) -> None:
        """
        Remove a panel from the dashboard.

        Parameters
        ----------
        panel_id : str
            ID of the panel to remove
        """
        if panel_id not in self._panels:
            raise ValueError(f"Panel with ID '{panel_id}' not found")

        # Remove from all linked panels
        for other_panel in self._panels.values():
            other_panel.unlink_from(panel_id)

        # Remove from dashboard
        del self._panels[panel_id]
        self._panel_order.remove(panel_id)

    def get_panel(self, panel_id: str) -> BasePanel:
        """
        Get a panel by ID.

        Parameters
        ----------
        panel_id : str
            ID of the panel to get

        Returns
        -------
        BasePanel
            The requested panel
        """
        if panel_id not in self._panels:
            raise ValueError(f"Panel with ID '{panel_id}' not found")
        return self._panels[panel_id]

    def link(self, source: str, target: str) -> None:
        """
        Link two panels so that selections in the source panel affect the target panel.

        Parameters
        ----------
        source : str
            ID of the source panel (where selections are made)
        target : str
            ID of the target panel (that will be affected by selections)
        """
        if source not in self._panels:
            raise ValueError(f"Source panel '{source}' not found")
        if target not in self._panels:
            raise ValueError(f"Target panel '{target}' not found")

        self._panels[target].link_to(source)

    def unlink(self, source: str, target: str) -> None:
        """
        Unlink two panels.

        Parameters
        ----------
        source : str
            ID of the source panel
        target : str
            ID of the target panel
        """
        if source not in self._panels:
            raise ValueError(f"Source panel '{source}' not found")
        if target not in self._panels:
            raise ValueError(f"Target panel '{target}' not found")

        self._panels[target].unlink_from(source)

    def set_global_selection(self, selection: Optional[np.ndarray]) -> None:
        """
        Set a global selection that affects all panels.

        Parameters
        ----------
        selection : np.ndarray or None
            Boolean array indicating selected cells, or None to clear selection
        """
        if selection is not None:
            if len(selection) != self._data_wrapper.adata.n_obs:
                raise ValueError(
                    f"Selection length ({len(selection)}) doesn't match number of observations ({self._data_wrapper.adata.n_obs})"
                )
            selection = np.array(selection, dtype=bool)

        self._global_selection = selection

        # Update all panels with the global selection
        for panel in self._panels.values():
            panel.selection = selection

    def get_global_selection(self) -> Optional[np.ndarray]:
        """
        Get the current global selection.

        Returns
        -------
        np.ndarray or None
            Boolean array indicating selected cells, or None if no selection
        """
        return self._global_selection.copy() if self._global_selection is not None else None

    def clear_selection(self) -> None:
        """Clear the global selection."""
        self.set_global_selection(None)

    def render_panel(self, panel_id: str) -> Any:
        """
        Render a specific panel.

        Parameters
        ----------
        panel_id : str
            ID of the panel to render

        Returns
        -------
        Any
            The rendered visualization
        """
        panel = self.get_panel(panel_id)
        return panel.render()

    def render_all_panels(self) -> Dict[str, Any]:
        """
        Render all panels in the dashboard.

        Returns
        -------
        dict
            Dictionary mapping panel IDs to their rendered visualizations
        """
        rendered_panels = {}
        for panel_id in self._panel_order:
            try:
                rendered_panels[panel_id] = self.render_panel(panel_id)
            except Exception as e:
                print(f"Error rendering panel '{panel_id}': {e}")
                rendered_panels[panel_id] = None

        return rendered_panels

    def export_code(self) -> str:
        """
        Export Python code that reproduces the current state of the dashboard.

        Returns
        -------
        str
            Python code that reproduces the current dashboard state
        """
        code_lines = [
            "# PySEE Dashboard Code Export",
            f"# Dashboard Title: {self._title}",
            f"# Generated on: {pd.Timestamp.now()}",
            "",
            "import numpy as np",
            "import pandas as pd",
            "import scanpy as sc",
            "from pysee import PySEE, UMAPPanel, ViolinPanel",
            "",
            "# Load your data (replace with your actual data loading code)",
            "# adata = sc.read('your_data.h5ad')",
            "",
            "# Create dashboard",
            f"app = PySEE(adata, title='{self._title}')",
            "",
        ]

        # Add panel creation code
        for panel_id in self._panel_order:
            panel = self._panels[panel_id]
            panel_info = panel.get_panel_info()

            if panel_info["panel_type"] == "UMAPPanel":
                embedding = panel.get_config("embedding")
                color = panel.get_config("color")
                code_lines.extend(
                    [
                        f"# Add UMAP panel: {panel_id}",
                        f"app.add_panel('{panel_id}', UMAPPanel(",
                        f"    panel_id='{panel_id}',",
                        f"    embedding='{embedding}',",
                        f"    color={'None' if color is None else repr(color)},",
                        f"    title='{panel.title}'",
                        "))",
                        "",
                    ]
                )
            elif panel_info["panel_type"] == "ViolinPanel":
                gene = panel.get_config("gene")
                group_by = panel.get_config("group_by")
                code_lines.extend(
                    [
                        f"# Add Violin panel: {panel_id}",
                        f"app.add_panel('{panel_id}', ViolinPanel(",
                        f"    panel_id='{panel_id}',",
                        f"    gene={'None' if gene is None else repr(gene)},",
                        f"    group_by={'None' if group_by is None else repr(group_by)},",
                        f"    title='{panel.title}'",
                        "))",
                        "",
                    ]
                )

        # Add linking code
        for panel_id in self._panel_order:
            panel = self._panels[panel_id]
            linked_panels = panel.linked_panels
            for linked_panel in linked_panels:
                code_lines.append(f"app.link('{linked_panel}', '{panel_id}')")

        if any(panel.linked_panels for panel in self._panels.values()):
            code_lines.append("")

        # Add selection code if there's a global selection
        if self._global_selection is not None:
            code_lines.extend(
                [
                    "# Set global selection",
                    f"selection_mask = np.array({self._global_selection.tolist()})",
                    "app.set_global_selection(selection_mask)",
                    "",
                ]
            )

        # Add rendering code
        code_lines.extend(
            [
                "# Render panels",
                "rendered_panels = app.render_all_panels()",
                "",
                "# Display panels (in Jupyter notebook)",
                "# for panel_id, fig in rendered_panels.items():",
                "#     if fig is not None:",
                "#         fig.show()",
            ]
        )

        return "\n".join(code_lines)

    def get_dashboard_info(self) -> Dict[str, Any]:
        """
        Get information about the dashboard.

        Returns
        -------
        dict
            Dictionary containing dashboard information
        """
        panel_info = {}
        for panel_id, panel in self._panels.items():
            panel_info[panel_id] = panel.get_panel_info()

        return {
            "title": self._title,
            "n_panels": len(self._panels),
            "panel_order": self._panel_order,
            "panels": panel_info,
            "has_global_selection": self._global_selection is not None,
            "n_selected_cells": (np.sum(self._global_selection) if self._global_selection is not None else 0),
            "data_summary": self._data_wrapper.get_summary_stats(),
        }

    def show(self) -> None:
        """
        Display the dashboard (placeholder for Jupyter integration).

        This method will be enhanced in future versions to provide
        proper Jupyter notebook display functionality.
        """
        print(f"PySEE Dashboard: {self._title}")
        print(f"Number of panels: {len(self._panels)}")
        print(f"Data: {self._data_wrapper.adata.n_obs} cells, {self._data_wrapper.adata.n_vars} genes")

        if self._global_selection is not None:
            n_selected = np.sum(self._global_selection)
            print(f"Selection: {n_selected} cells selected")

        print("\nPanels:")
        for panel_id in self._panel_order:
            panel = self._panels[panel_id]
            print(f"  - {panel_id}: {panel.__class__.__name__} - {panel.title}")

    def __repr__(self) -> str:
        """String representation of the dashboard."""
        return f"PySEE(title='{self._title}', n_panels={len(self._panels)})"
