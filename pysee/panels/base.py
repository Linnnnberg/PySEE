"""
Base panel class for PySEE visualization panels.

This module defines the abstract base class that all PySEE panels must inherit from.
"""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd

from ..core.data import AnnDataWrapper
from ..utils.export import PublicationExporter
from .advanced_selection import AdvancedSelectionManager
from .selection_types import SelectionHistory, SelectionMode, SelectionType


class BasePanel(ABC):
    """
    Abstract base class for all PySEE visualization panels.

    All panels must inherit from this class and implement the required methods.
    """

    def __init__(self, panel_id: str, title: Optional[str] = None):
        """
        Initialize the base panel.

        Parameters
        ----------
        panel_id : str
            Unique identifier for the panel
        title : str, optional
            Display title for the panel
        """
        self.panel_id = panel_id
        self.title = title or panel_id
        self._data_wrapper: Optional[AnnDataWrapper] = None
        self._selection: Optional[np.ndarray] = None
        self._linked_panels: List[str] = []
        self._config: Dict[str, Any] = {}
        self._selection_manager: Optional[AdvancedSelectionManager] = None
        self._advanced_selection_enabled = False

    @property
    def data_wrapper(self) -> Optional[AnnDataWrapper]:
        """Get the data wrapper associated with this panel."""
        return self._data_wrapper

    @data_wrapper.setter
    def data_wrapper(self, wrapper: AnnDataWrapper) -> None:
        """Set the data wrapper for this panel."""
        self._data_wrapper = wrapper
        # Initialize selection manager if advanced selection is enabled
        if self._advanced_selection_enabled and wrapper is not None:
            self._selection_manager = AdvancedSelectionManager(wrapper)
        self._on_data_changed()

    @property
    def selection(self) -> Optional[np.ndarray]:
        """Get the current selection mask."""
        return self._selection

    @selection.setter
    def selection(self, mask: Optional[np.ndarray]) -> None:
        """Set the selection mask."""
        self._selection = mask
        self._on_selection_changed()

    @property
    def linked_panels(self) -> List[str]:
        """Get list of linked panel IDs."""
        return self._linked_panels.copy()

    def link_to(self, panel_id: str) -> None:
        """Link this panel to another panel."""
        if panel_id not in self._linked_panels:
            self._linked_panels.append(panel_id)

    def unlink_from(self, panel_id: str) -> None:
        """Unlink this panel from another panel."""
        if panel_id in self._linked_panels:
            self._linked_panels.remove(panel_id)

    def get_config(self, key: str, default: Any = None) -> Any:
        """Get a configuration value."""
        return self._config.get(key, default)

    def set_config(self, key: str, value: Any) -> None:
        """Set a configuration value."""
        self._config[key] = value
        self._on_config_changed()

    def update_config(self, config: Dict[str, Any]) -> None:
        """Update multiple configuration values."""
        self._config.update(config)
        self._on_config_changed()

    @abstractmethod
    def render(self) -> Any:
        """
        Render the panel visualization.

        Returns
        -------
        Any
            The rendered visualization (e.g., plotly figure, bokeh plot, etc.)
        """
        pass

    @abstractmethod
    def get_selection_code(self) -> str:
        """
        Generate Python code for the current selection.

        Returns
        -------
        str
            Python code that reproduces the current selection
        """
        pass

    def _on_data_changed(self) -> None:
        """Called when the data wrapper changes. Override in subclasses."""
        pass

    def _on_selection_changed(self) -> None:
        """Called when the selection changes. Override in subclasses."""
        pass

    def _on_config_changed(self) -> None:
        """Called when the configuration changes. Override in subclasses."""
        pass

    def enable_advanced_selection(self) -> None:
        """Enable advanced selection tools for this panel."""
        self._advanced_selection_enabled = True
        if self._data_wrapper is not None:
            self._selection_manager = AdvancedSelectionManager(self._data_wrapper)
            # Set up callback to sync with panel selection
            self._selection_manager.add_callback(self._on_advanced_selection_changed)

    def disable_advanced_selection(self) -> None:
        """Disable advanced selection tools for this panel."""
        self._advanced_selection_enabled = False
        self._selection_manager = None

    def is_advanced_selection_enabled(self) -> bool:
        """Check if advanced selection is enabled."""
        return self._advanced_selection_enabled and self._selection_manager is not None

    def set_selection_type(self, selection_type: SelectionType) -> None:
        """Set the selection tool type."""
        if self._selection_manager is not None:
            self._selection_manager.set_selection_type(selection_type)

    def set_selection_mode(self, mode: SelectionMode) -> None:
        """Set the selection mode."""
        if self._selection_manager is not None:
            self._selection_manager.set_selection_mode(mode)

    def undo_selection(self) -> bool:
        """Undo last selection. Returns True if successful."""
        if self._selection_manager is not None:
            return self._selection_manager.undo_selection()
        return False

    def redo_selection(self) -> bool:
        """Redo next selection. Returns True if successful."""
        if self._selection_manager is not None:
            return self._selection_manager.redo_selection()
        return False

    def clear_selection(self) -> None:
        """Clear current selection."""
        if self._selection_manager is not None:
            self._selection_manager.clear_selection()
        else:
            self._selection = None

    def select_all(self) -> None:
        """Select all points."""
        if self._selection_manager is not None:
            self._selection_manager.select_all()

    def invert_selection(self) -> None:
        """Invert current selection."""
        if self._selection_manager is not None:
            self._selection_manager.invert_selection()

    def get_selection_statistics(self) -> Dict[str, Any]:
        """Get selection statistics."""
        if self._selection_manager is not None:
            return self._selection_manager.get_selection_statistics()
        elif self._selection is not None:
            from .selection_types import SelectionStatistics

            return SelectionStatistics.calculate_stats(self._selection, len(self._selection))
        else:
            return {"error": "No selection available"}

    def _on_advanced_selection_changed(self, selection: np.ndarray) -> None:
        """Callback for when advanced selection changes."""
        self._selection = selection
        self._on_selection_changed()

    def validate_data(self) -> bool:
        """
        Validate that the panel has the required data.

        Returns
        -------
        bool
            True if data is valid, False otherwise
        """
        if self._data_wrapper is None:
            return False

        # Check if data wrapper has required data
        return self._check_data_requirements()

    @abstractmethod
    def _check_data_requirements(self) -> bool:
        """
        Check if the data wrapper meets the panel's requirements.

        Returns
        -------
        bool
            True if requirements are met, False otherwise
        """
        pass

    def export(
        self,
        output_path: Union[str, Path],
        format: str = "png",
        template: str = "custom",
        width: Optional[int] = None,
        height: Optional[int] = None,
        dpi: Optional[int] = None,
        title: Optional[str] = None,
        caption: Optional[str] = None,
        config_overrides: Optional[Dict[str, Any]] = None,
    ) -> str:
        """
        Export the panel visualization in publication-ready format.

        Parameters
        ----------
        output_path : str or Path
            Output file path
        format : str, default 'png'
            Export format ('png', 'svg', 'pdf', 'html')
        template : str, default 'custom'
            Journal template ('nature', 'science', 'cell', 'custom')
        width : int, optional
            Width in pixels (overrides template)
        height : int, optional
            Height in pixels (overrides template)
        dpi : int, optional
            DPI for raster formats (overrides template)
        title : str, optional
            Figure title (uses panel title if None)
        caption : str, optional
            Figure caption
        config_overrides : dict, optional
            Configuration overrides for styling

        Returns
        -------
        str
            Path to the exported file
        """
        if not self.validate_data():
            raise ValueError("Panel data requirements not met")

        # Render the figure
        figure = self.render()

        # Use panel title if no title provided
        if title is None:
            title = self.title

        # Export using PublicationExporter
        exporter = PublicationExporter(template)
        return exporter.export_panel(
            figure, output_path, format, width, height, dpi, title, caption, config_overrides
        )

    def export_as_bytes(
        self,
        format: str = "png",
        template: str = "custom",
        width: Optional[int] = None,
        height: Optional[int] = None,
        dpi: Optional[int] = None,
    ) -> bytes:
        """
        Export the panel visualization as bytes.

        Parameters
        ----------
        format : str, default 'png'
            Export format ('png', 'svg')
        template : str, default 'custom'
            Journal template
        width : int, optional
            Width in pixels
        height : int, optional
            Height in pixels
        dpi : int, optional
            DPI for raster formats

        Returns
        -------
        bytes
            Figure as bytes
        """
        if not self.validate_data():
            raise ValueError("Panel data requirements not met")

        # Render the figure
        figure = self.render()

        # Export using PublicationExporter
        exporter = PublicationExporter(template)
        return exporter.get_figure_as_bytes(figure, format, width, height, dpi)

    def export_as_base64(self, format: str = "png", template: str = "custom", **kwargs: Any) -> str:
        """
        Export the panel visualization as base64 string.

        Parameters
        ----------
        format : str, default 'png'
            Export format
        template : str, default 'custom'
            Journal template
        **kwargs
            Additional arguments passed to get_figure_as_bytes

        Returns
        -------
        str
            Base64 encoded figure
        """
        if not self.validate_data():
            raise ValueError("Panel data requirements not met")

        # Render the figure
        figure = self.render()

        # Export using PublicationExporter
        exporter = PublicationExporter(template)
        return exporter.get_figure_as_base64(figure, format, **kwargs)

    def configure_for_export(self, config_overrides: Dict[str, Any]) -> "BasePanel":
        """
        Configure panel for export with temporary overrides.

        This method creates a copy of the panel with modified configuration
        for export purposes without affecting the original panel.

        Parameters
        ----------
        config_overrides : dict
            Configuration overrides to apply

        Returns
        -------
        BasePanel
            New panel instance with overridden configuration
        """
        # Create a copy of the panel
        import copy

        panel_copy = copy.deepcopy(self)

        # Apply configuration overrides
        panel_copy.update_config(config_overrides)

        return panel_copy

    def get_export_config(self) -> Dict[str, Any]:
        """
        Get current export configuration.

        Returns
        -------
        dict
            Current export configuration
        """
        return {
            "panel_config": self._config.copy(),
            "available_templates": PublicationExporter().get_available_templates(),
            "recommended_formats": ["png", "svg", "pdf", "html"],
            "default_dpi": 300,
            "default_width_mm": 180,
            "default_height_mm": 120,
            "current_template_config": PublicationExporter().config,
        }

    def get_publication_info(self) -> Dict[str, Any]:
        """
        Get publication-ready information about this panel.

        Returns
        -------
        dict
            Dictionary containing publication information
        """
        return {
            "panel_id": self.panel_id,
            "title": self.title,
            "panel_type": self.__class__.__name__,
            "has_data": self._data_wrapper is not None,
            "has_selection": self._selection is not None,
            "available_templates": PublicationExporter().get_available_templates(),
            "recommended_formats": ["png", "svg", "pdf", "html"],
            "default_dpi": 300,
            "default_width_mm": 180,
            "default_height_mm": 120,
        }

    def get_panel_info(self) -> Dict[str, Any]:
        """
        Get information about this panel.

        Returns
        -------
        dict
            Dictionary containing panel information
        """
        return {
            "panel_id": self.panel_id,
            "title": self.title,
            "panel_type": self.__class__.__name__,
            "linked_panels": self._linked_panels,
            "config": self._config,
            "has_data": self._data_wrapper is not None,
            "has_selection": self._selection is not None,
        }

    def __repr__(self) -> str:
        """String representation of the panel."""
        return f"{self.__class__.__name__}(id='{self.panel_id}', title='{self.title}')"
