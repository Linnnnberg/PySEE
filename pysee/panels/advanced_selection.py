"""
Advanced selection manager for PySEE panels.

This module provides the AdvancedSelectionManager class that handles
complex selection operations, history management, and coordinate-based
selection calculations.
"""

from datetime import datetime
from typing import Any, Callable, Dict, List, Optional, Tuple

import matplotlib.path as mpath
import numpy as np
from shapely.geometry import Point, Polygon
from shapely.geometry.polygon import LinearRing

from ..core.data import AnnDataWrapper
from .selection_types import (
    SelectionHistory,
    SelectionMode,
    SelectionState,
    SelectionStatistics,
    SelectionType,
)


class AdvancedSelectionManager:
    """
    Manages advanced selection operations for PySEE panels.

    This class handles different types of selections (rectangular, lasso, polygon),
    selection modes (add, subtract, intersect), and maintains selection history
    with undo/redo capabilities.
    """

    def __init__(self, data_wrapper: AnnDataWrapper):
        """
        Initialize the advanced selection manager.

        Parameters
        ----------
        data_wrapper : AnnDataWrapper
            The data wrapper containing the dataset
        """
        self._data_wrapper = data_wrapper
        self._history = SelectionHistory()
        self._current_selection: Optional[np.ndarray] = None
        self._selection_type = SelectionType.POINT
        self._selection_mode = SelectionMode.REPLACE
        self._callbacks: List[Callable[[np.ndarray], None]] = []

    @property
    def current_selection(self) -> Optional[np.ndarray]:
        """Get the current selection mask."""
        return self._current_selection

    @property
    def selection_type(self) -> SelectionType:
        """Get the current selection type."""
        return self._selection_type

    @property
    def selection_mode(self) -> SelectionMode:
        """Get the current selection mode."""
        return self._selection_mode

    def set_selection_type(self, selection_type: SelectionType) -> None:
        """
        Set the type of selection tool.

        Parameters
        ----------
        selection_type : SelectionType
            The selection tool type to use
        """
        self._selection_type = selection_type

    def set_selection_mode(self, mode: SelectionMode) -> None:
        """
        Set how new selections interact with existing selections.

        Parameters
        ----------
        mode : SelectionMode
            The selection mode to use
        """
        self._selection_mode = mode

    def add_callback(self, callback: Callable[[np.ndarray], None]) -> None:
        """
        Add a callback function to be called when selection changes.

        Parameters
        ----------
        callback : Callable
            Function to call with new selection mask
        """
        self._callbacks.append(callback)

    def remove_callback(self, callback: Callable[[np.ndarray], None]) -> None:
        """
        Remove a callback function.

        Parameters
        ----------
        callback : Callable
            Function to remove
        """
        if callback in self._callbacks:
            self._callbacks.remove(callback)

    def _notify_callbacks(self, selection: np.ndarray) -> None:
        """Notify all registered callbacks of selection change."""
        for callback in self._callbacks:
            try:
                callback(selection)
            except Exception as e:
                print(f"Warning: Selection callback failed: {e}")

    def create_point_selection(self, indices: List[int]) -> np.ndarray:
        """
        Create selection from point indices.

        Parameters
        ----------
        indices : List[int]
            List of point indices to select

        Returns
        -------
        np.ndarray
            Boolean selection mask
        """
        n_points = self._data_wrapper.adata.n_obs
        selection = np.zeros(n_points, dtype=bool)

        # Ensure indices are within bounds
        valid_indices = [i for i in indices if 0 <= i < n_points]
        selection[valid_indices] = True

        return selection

    def create_rectangular_selection(
        self, bounds: Tuple[float, float, float, float], coords: np.ndarray
    ) -> np.ndarray:
        """
        Create selection from rectangular bounds.

        Parameters
        ----------
        bounds : Tuple[float, float, float, float]
            Rectangle bounds as (xmin, ymin, xmax, ymax)
        coords : np.ndarray
            2D coordinates array (n_points, 2)

        Returns
        -------
        np.ndarray
            Boolean selection mask
        """
        xmin, ymin, xmax, ymax = bounds

        # Ensure proper ordering
        if xmin > xmax:
            xmin, xmax = xmax, xmin
        if ymin > ymax:
            ymin, ymax = ymax, ymin

        # Create selection mask
        x_coords = coords[:, 0]
        y_coords = coords[:, 1]

        selection = (
            (x_coords >= xmin) & (x_coords <= xmax) & (y_coords >= ymin) & (y_coords <= ymax)
        )

        return selection

    def create_lasso_selection(
        self, path: List[Tuple[float, float]], coords: np.ndarray
    ) -> np.ndarray:
        """
        Create selection from lasso path.

        Parameters
        ----------
        path : List[Tuple[float, float]]
            List of (x, y) coordinates defining the lasso path
        coords : np.ndarray
            2D coordinates array (n_points, 2)

        Returns
        -------
        np.ndarray
            Boolean selection mask
        """
        if len(path) < 3:
            # Need at least 3 points to form a polygon
            return np.zeros(len(coords), dtype=bool)

        # Use matplotlib path for point-in-polygon test
        lasso_path = mpath.Path(path)
        selection = lasso_path.contains_points(coords)

        return selection

    def create_polygon_selection(
        self, vertices: List[Tuple[float, float]], coords: np.ndarray
    ) -> np.ndarray:
        """
        Create selection from polygon vertices.

        Parameters
        ----------
        vertices : List[Tuple[float, float]]
            List of (x, y) coordinates defining polygon vertices
        coords : np.ndarray
            2D coordinates array (n_points, 2)

        Returns
        -------
        np.ndarray
            Boolean selection mask
        """
        if len(vertices) < 3:
            # Need at least 3 vertices to form a polygon
            return np.zeros(len(coords), dtype=bool)

        try:
            # Use shapely for robust polygon operations
            polygon = Polygon(vertices)
            if not polygon.is_valid:
                # Try to fix invalid polygon
                polygon = polygon.buffer(0)

            # Check which points are inside the polygon
            selection = np.array([polygon.contains(Point(x, y)) for x, y in coords])

        except Exception:
            # Fallback to matplotlib path if shapely fails
            poly_path = mpath.Path(vertices)
            selection = poly_path.contains_points(coords)

        return selection

    def create_circle_selection(
        self, center: Tuple[float, float], radius: float, coords: np.ndarray
    ) -> np.ndarray:
        """
        Create selection from circle.

        Parameters
        ----------
        center : Tuple[float, float]
            Center coordinates (x, y)
        radius : float
            Circle radius
        coords : np.ndarray
            2D coordinates array (n_points, 2)

        Returns
        -------
        np.ndarray
            Boolean selection mask
        """
        cx, cy = center

        # Calculate distances from center
        distances = np.sqrt((coords[:, 0] - cx) ** 2 + (coords[:, 1] - cy) ** 2)
        selection = distances <= radius

        return selection

    def combine_selections(
        self, new_selection: np.ndarray, existing_selection: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Combine new selection with existing selection based on current mode.

        Parameters
        ----------
        new_selection : np.ndarray
            New selection mask
        existing_selection : np.ndarray, optional
            Existing selection mask (uses current selection if None)

        Returns
        -------
        np.ndarray
            Combined selection mask
        """
        if existing_selection is None:
            existing_selection = self._current_selection

        if existing_selection is None or self._selection_mode == SelectionMode.REPLACE:
            return new_selection.copy()

        if self._selection_mode == SelectionMode.ADD:
            result = existing_selection | new_selection
            return np.array(result, dtype=bool)
        elif self._selection_mode == SelectionMode.SUBTRACT:
            result = existing_selection & ~new_selection
            return np.array(result, dtype=bool)
        elif self._selection_mode == SelectionMode.INTERSECT:
            result = existing_selection & new_selection
            return np.array(result, dtype=bool)
        else:
            return new_selection.copy()

    def apply_selection(
        self,
        selection: np.ndarray,
        selection_type: Optional[SelectionType] = None,
        add_to_history: bool = True,
    ) -> None:
        """
        Apply a new selection and update history.

        Parameters
        ----------
        selection : np.ndarray
            Selection mask to apply
        selection_type : SelectionType, optional
            Type of selection (uses current if None)
        add_to_history : bool
            Whether to add this selection to history
        """
        if selection_type is None:
            selection_type = self._selection_type

        # Combine with existing selection based on mode
        final_selection = self.combine_selections(selection)

        # Update current selection
        self._current_selection = final_selection

        # Add to history
        if add_to_history:
            selection_state = SelectionState(
                selection=final_selection.copy(),
                selection_type=selection_type,
                selection_mode=self._selection_mode,
                timestamp=datetime.now(),
                metadata={
                    "n_selected": int(np.sum(final_selection)),
                    "n_total": len(final_selection),
                },
            )
            self._history.add_selection(selection_state)

        # Notify callbacks
        self._notify_callbacks(final_selection)

    def undo_selection(self) -> bool:
        """
        Undo to previous selection.

        Returns
        -------
        bool
            True if undo was successful
        """
        previous_state = self._history.undo()
        if previous_state is not None:
            self._current_selection = previous_state.selection.copy()
            self._notify_callbacks(self._current_selection)
            return True
        return False

    def redo_selection(self) -> bool:
        """
        Redo to next selection.

        Returns
        -------
        bool
            True if redo was successful
        """
        next_state = self._history.redo()
        if next_state is not None:
            self._current_selection = next_state.selection.copy()
            self._notify_callbacks(self._current_selection)
            return True
        return False

    def clear_selection(self, add_to_history: bool = True) -> None:
        """
        Clear current selection.

        Parameters
        ----------
        add_to_history : bool
            Whether to add clear operation to history
        """
        if self._current_selection is not None:
            empty_selection = np.zeros_like(self._current_selection, dtype=bool)
            self.apply_selection(empty_selection, add_to_history=add_to_history)

    def select_all(self, add_to_history: bool = True) -> None:
        """
        Select all points.

        Parameters
        ----------
        add_to_history : bool
            Whether to add select all operation to history
        """
        n_points = self._data_wrapper.adata.n_obs
        full_selection = np.ones(n_points, dtype=bool)
        self.apply_selection(full_selection, add_to_history=add_to_history)

    def invert_selection(self, add_to_history: bool = True) -> None:
        """
        Invert current selection.

        Parameters
        ----------
        add_to_history : bool
            Whether to add invert operation to history
        """
        if self._current_selection is not None:
            inverted_selection = ~self._current_selection
            self.apply_selection(inverted_selection, add_to_history=add_to_history)

    def get_selection_statistics(self) -> Dict[str, Any]:
        """
        Get statistics about current selection.

        Returns
        -------
        Dict[str, Any]
            Selection statistics
        """
        if self._current_selection is None:
            return {"error": "No selection available"}

        stats = SelectionStatistics.calculate_stats(
            self._current_selection, self._data_wrapper.adata.n_obs
        )

        # Add history information
        history_summary = self._history.get_history_summary()
        stats.update(
            {
                "history": history_summary,
                "current_mode": self._selection_mode.value,
                "current_type": self._selection_type.value,
            }
        )

        return stats

    def get_selected_indices(self) -> np.ndarray:
        """
        Get indices of selected points.

        Returns
        -------
        np.ndarray
            Array of selected indices
        """
        if self._current_selection is None:
            return np.array([], dtype=int)

        return np.where(self._current_selection)[0]

    def get_selected_data(self) -> Optional[Any]:
        """
        Get subset of data for selected points.

        Returns
        -------
        AnnDataWrapper or None
            Data subset for selected points
        """
        if self._current_selection is None or not np.any(self._current_selection):
            return None

        return self._data_wrapper.get_cell_subset(self._current_selection)
