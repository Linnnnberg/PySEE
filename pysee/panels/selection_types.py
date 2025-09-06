"""
Advanced selection types and modes for PySEE panels.

This module defines the selection types, modes, and core infrastructure
for advanced selection tools in PySEE visualization panels.
"""

from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np


class SelectionType(Enum):
    """Types of selection tools available."""

    POINT = "point"
    RECTANGULAR = "rectangular"
    LASSO = "lasso"
    POLYGON = "polygon"
    CIRCLE = "circle"


class SelectionMode(Enum):
    """How new selections interact with existing selections."""

    REPLACE = "replace"  # Replace current selection
    ADD = "add"  # Add to current selection (union)
    SUBTRACT = "subtract"  # Remove from current selection (difference)
    INTERSECT = "intersect"  # Keep only intersection


@dataclass
class SelectionState:
    """Represents a selection state with metadata."""

    selection: np.ndarray
    selection_type: SelectionType
    selection_mode: SelectionMode
    timestamp: datetime
    metadata: Dict[str, Any]

    def __post_init__(self) -> None:
        """Validate selection state after initialization."""
        if not isinstance(self.selection, np.ndarray):
            self.selection = np.array(self.selection, dtype=bool)

        if self.selection.dtype != bool:
            self.selection = self.selection.astype(bool)


class SelectionHistory:
    """Manages selection history with undo/redo capabilities."""

    def __init__(self, max_history: int = 50):
        """
        Initialize selection history.

        Parameters
        ----------
        max_history : int
            Maximum number of selection states to keep in history
        """
        self._history: List[SelectionState] = []
        self._current_index = -1
        self._max_history = max_history

    def add_selection(self, selection_state: SelectionState) -> None:
        """
        Add a new selection state to history.

        Parameters
        ----------
        selection_state : SelectionState
            The selection state to add
        """
        # Remove any future history if we're not at the end
        if self._current_index < len(self._history) - 1:
            self._history = self._history[: self._current_index + 1]

        # Add new state
        self._history.append(selection_state)
        self._current_index = len(self._history) - 1

        # Trim history if too long
        if len(self._history) > self._max_history:
            self._history = self._history[-self._max_history :]
            self._current_index = len(self._history) - 1

    def undo(self) -> Optional[SelectionState]:
        """
        Undo to previous selection state.

        Returns
        -------
        SelectionState or None
            Previous selection state, or None if no undo available
        """
        if self.can_undo():
            self._current_index -= 1
            return self._history[self._current_index]
        return None

    def redo(self) -> Optional[SelectionState]:
        """
        Redo to next selection state.

        Returns
        -------
        SelectionState or None
            Next selection state, or None if no redo available
        """
        if self.can_redo():
            self._current_index += 1
            return self._history[self._current_index]
        return None

    def can_undo(self) -> bool:
        """Check if undo is possible."""
        return self._current_index > 0

    def can_redo(self) -> bool:
        """Check if redo is possible."""
        return self._current_index < len(self._history) - 1

    def get_current_state(self) -> Optional[SelectionState]:
        """Get current selection state."""
        if 0 <= self._current_index < len(self._history):
            return self._history[self._current_index]
        return None

    def clear(self) -> None:
        """Clear all history."""
        self._history.clear()
        self._current_index = -1

    def get_history_summary(self) -> Dict[str, Any]:
        """Get summary of selection history."""
        return {
            "total_states": len(self._history),
            "current_index": self._current_index,
            "can_undo": self.can_undo(),
            "can_redo": self.can_redo(),
            "max_history": self._max_history,
        }


class SelectionStatistics:
    """Calculate statistics for selections."""

    @staticmethod
    def calculate_stats(selection: np.ndarray, total_points: int) -> Dict[str, Any]:
        """
        Calculate selection statistics.

        Parameters
        ----------
        selection : np.ndarray
            Boolean selection array
        total_points : int
            Total number of points in dataset

        Returns
        -------
        Dict[str, Any]
            Selection statistics
        """
        n_selected = np.sum(selection)
        n_total = len(selection)

        return {
            "n_selected": int(n_selected),
            "n_total": int(n_total),
            "n_unselected": int(n_total - n_selected),
            "selection_percentage": float(n_selected / n_total * 100) if n_total > 0 else 0.0,
            "selection_ratio": float(n_selected / n_total) if n_total > 0 else 0.0,
            "has_selection": bool(n_selected > 0),
            "is_empty": bool(n_selected == 0),
            "is_full": bool(n_selected == n_total),
            "selection_indices": (
                np.where(selection)[0].tolist() if n_selected < 1000 else "too_many_to_list"
            ),
        }

    @staticmethod
    def compare_selections(selection1: np.ndarray, selection2: np.ndarray) -> Dict[str, Any]:
        """
        Compare two selections.

        Parameters
        ----------
        selection1, selection2 : np.ndarray
            Boolean selection arrays to compare

        Returns
        -------
        Dict[str, Any]
            Comparison statistics
        """
        if len(selection1) != len(selection2):
            raise ValueError("Selections must have the same length")

        intersection = selection1 & selection2
        union = selection1 | selection2
        difference1 = selection1 & ~selection2
        difference2 = selection2 & ~selection1

        return {
            "intersection_size": int(np.sum(intersection)),
            "union_size": int(np.sum(union)),
            "difference1_size": int(np.sum(difference1)),  # in sel1 but not sel2
            "difference2_size": int(np.sum(difference2)),  # in sel2 but not sel1
            "jaccard_index": (
                float(np.sum(intersection) / np.sum(union)) if np.sum(union) > 0 else 0.0
            ),
            "overlap_percentage": (
                float(np.sum(intersection) / min(np.sum(selection1), np.sum(selection2)) * 100)
                if min(np.sum(selection1), np.sum(selection2)) > 0
                else 0.0
            ),
        }
