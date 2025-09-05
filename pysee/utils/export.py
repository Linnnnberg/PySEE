"""
Publication-quality export utilities for PySEE visualizations.

This module provides functionality for exporting PySEE visualizations in
publication-ready formats with high resolution and journal-specific styling.
"""

import os
import io
from typing import Any, Dict, List, Optional, Union, Tuple
from pathlib import Path
import plotly.graph_objects as go
from plotly.io import to_image, to_html, write_image, write_html
import base64


class PublicationExporter:
    """
    Export PySEE visualizations in publication-ready formats.

    Supports high-resolution image export, journal-specific styling,
    and batch export functionality.
    """

    # Journal-specific styling templates
    JOURNAL_TEMPLATES = {
        "nature": {
            "width": 180,  # mm
            "height": 120,  # mm
            "dpi": 300,
            "font_family": "Arial",
            "font_size": 8,
            "line_width": 0.5,
            "marker_size": 3,
            "color_scheme": "viridis",
        },
        "science": {
            "width": 180,  # mm
            "height": 120,  # mm
            "dpi": 300,
            "font_family": "Arial",
            "font_size": 7,
            "line_width": 0.5,
            "marker_size": 2.5,
            "color_scheme": "viridis",
        },
        "cell": {
            "width": 180,  # mm
            "height": 120,  # mm
            "dpi": 300,
            "font_family": "Arial",
            "font_size": 8,
            "line_width": 0.5,
            "marker_size": 3,
            "color_scheme": "viridis",
        },
        "custom": {
            "width": 180,  # mm
            "height": 120,  # mm
            "dpi": 300,
            "font_family": "Arial",
            "font_size": 8,
            "line_width": 0.5,
            "marker_size": 3,
            "color_scheme": "viridis",
        },
    }

    def __init__(self, template: str = "custom"):
        """
        Initialize the publication exporter.

        Parameters
        ----------
        template : str, default 'custom'
            Journal template to use ('nature', 'science', 'cell', 'custom')
        """
        self.template = template
        self.config = self.JOURNAL_TEMPLATES[template].copy()

    def export_panel(
        self,
        figure: go.Figure,
        output_path: Union[str, Path],
        format: str = "png",
        width: Optional[int] = None,
        height: Optional[int] = None,
        dpi: Optional[int] = None,
        title: Optional[str] = None,
        caption: Optional[str] = None,
        config_overrides: Optional[Dict[str, Any]] = None,
    ) -> str:
        """
        Export a single panel to publication-ready format.

        Parameters
        ----------
        figure : plotly.graph_objects.Figure
            Plotly figure to export
        output_path : str or Path
            Output file path
        format : str, default 'png'
            Export format ('png', 'svg', 'pdf', 'html')
        width : int, optional
            Width in pixels (overrides template)
        height : int, optional
            Height in pixels (overrides template)
        dpi : int, optional
            DPI for raster formats (overrides template)
        title : str, optional
            Figure title to add
        caption : str, optional
            Figure caption to add
        config_overrides : dict, optional
            Configuration overrides for styling

        Returns
        -------
        str
            Path to the exported file
        """
        # Convert mm to pixels if needed
        if width is None:
            width = int(self.config["width"] * self.config["dpi"] / 25.4)
        if height is None:
            height = int(self.config["height"] * self.config["dpi"] / 25.4)
        if dpi is None:
            dpi = self.config["dpi"]

        # Apply publication styling with config overrides
        styled_figure = self._apply_publication_styling(figure, title, caption, config_overrides)

        # Ensure output directory exists
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Export based on format
        if format.lower() == "png":
            write_image(
                styled_figure,
                str(output_path),
                width=width,
                height=height,
                scale=dpi / 72,  # Plotly uses 72 DPI as base
                format="png",
            )
        elif format.lower() == "svg":
            write_image(styled_figure, str(output_path), width=width, height=height, format="svg")
        elif format.lower() == "pdf":
            write_image(styled_figure, str(output_path), width=width, height=height, format="pdf")
        elif format.lower() == "html":
            write_html(styled_figure, str(output_path), config={"displayModeBar": False})
        else:
            raise ValueError(f"Unsupported format: {format}")

        return str(output_path)

    def export_batch(
        self,
        panels: List[Tuple[go.Figure, str]],  # (figure, panel_name)
        output_dir: Union[str, Path],
        format: str = "png",
        prefix: str = "panel",
        **kwargs,
    ) -> List[str]:
        """
        Export multiple panels in batch.

        Parameters
        ----------
        panels : List[Tuple[plotly.graph_objects.Figure, str]]
            List of (figure, panel_name) tuples
        output_dir : str or Path
            Output directory
        format : str, default 'png'
            Export format
        prefix : str, default 'panel'
            Prefix for output filenames
        **kwargs
            Additional arguments passed to export_panel

        Returns
        -------
        List[str]
            List of exported file paths
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        exported_files = []
        for i, (figure, panel_name) in enumerate(panels):
            filename = f"{prefix}_{i+1:02d}_{panel_name}.{format}"
            output_path = output_dir / filename

            exported_path = self.export_panel(figure, output_path, format, **kwargs)
            exported_files.append(exported_path)

        return exported_files

    def get_figure_as_bytes(
        self,
        figure: go.Figure,
        format: str = "png",
        width: Optional[int] = None,
        height: Optional[int] = None,
        dpi: Optional[int] = None,
    ) -> bytes:
        """
        Get figure as bytes for in-memory operations.

        Parameters
        ----------
        figure : plotly.graph_objects.Figure
            Plotly figure
        format : str, default 'png'
            Export format
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
        # Apply styling
        styled_figure = self._apply_publication_styling(figure)

        # Convert mm to pixels if needed
        if width is None:
            width = int(self.config["width"] * self.config["dpi"] / 25.4)
        if height is None:
            height = int(self.config["height"] * self.config["dpi"] / 25.4)
        if dpi is None:
            dpi = self.config["dpi"]

        # Convert to bytes
        if format.lower() == "png":
            return to_image(styled_figure, format="png", width=width, height=height, scale=dpi / 72)
        elif format.lower() == "svg":
            return to_image(styled_figure, format="svg", width=width, height=height)
        else:
            raise ValueError(f"Unsupported format for bytes: {format}")

    def get_figure_as_base64(self, figure: go.Figure, format: str = "png", **kwargs) -> str:
        """
        Get figure as base64 string for embedding in HTML.

        Parameters
        ----------
        figure : plotly.graph_objects.Figure
            Plotly figure
        format : str, default 'png'
            Export format
        **kwargs
            Additional arguments passed to get_figure_as_bytes

        Returns
        -------
        str
            Base64 encoded figure
        """
        figure_bytes = self.get_figure_as_bytes(figure, format, **kwargs)
        return base64.b64encode(figure_bytes).decode("utf-8")

    def _apply_publication_styling(
        self,
        figure: go.Figure,
        title: Optional[str] = None,
        caption: Optional[str] = None,
        config_overrides: Optional[Dict[str, Any]] = None,
    ) -> go.Figure:
        """
        Apply publication-quality styling to a figure.

        Parameters
        ----------
        figure : plotly.graph_objects.Figure
            Input figure
        title : str, optional
            Figure title
        caption : str, optional
            Figure caption
        config_overrides : dict, optional
            Configuration overrides for styling

        Returns
        -------
        plotly.graph_objects.Figure
            Styled figure
        """
        # Create a copy to avoid modifying the original
        import copy

        styled_figure = copy.deepcopy(figure)

        # Apply config overrides if provided
        if config_overrides:
            self.config.update(config_overrides)

        # Update layout with publication styling
        styled_figure.update_layout(
            font=dict(
                family=self.config["font_family"], size=self.config["font_size"], color="black"
            ),
            plot_bgcolor="white",
            paper_bgcolor="white",
            margin=dict(l=60, r=60, t=60, b=60),
            showlegend=True,
            legend=dict(
                font=dict(size=self.config["font_size"]),
                bgcolor="rgba(255,255,255,0.8)",
                bordercolor="black",
                borderwidth=0.5,
            ),
        )

        # Update axes styling
        styled_figure.update_xaxes(
            showline=True,
            linewidth=self.config["line_width"],
            linecolor="black",
            mirror=True,
            gridcolor="lightgray",
            gridwidth=0.5,
        )

        styled_figure.update_yaxes(
            showline=True,
            linewidth=self.config["line_width"],
            linecolor="black",
            mirror=True,
            gridcolor="lightgray",
            gridwidth=0.5,
        )

        # Add title if provided
        if title:
            styled_figure.update_layout(
                title=dict(
                    text=title,
                    font=dict(
                        family=self.config["font_family"],
                        size=self.config["font_size"] + 2,
                        color="black",
                    ),
                    x=0.5,
                    xanchor="center",
                )
            )

        # Add caption if provided
        if caption:
            styled_figure.add_annotation(
                text=caption,
                xref="paper",
                yref="paper",
                x=0.5,
                y=-0.1,
                showarrow=False,
                font=dict(
                    family=self.config["font_family"],
                    size=self.config["font_size"] - 1,
                    color="black",
                ),
                xanchor="center",
            )

        return styled_figure

    def update_template(self, template: str) -> None:
        """
        Update the journal template.

        Parameters
        ----------
        template : str
            Template name ('nature', 'science', 'cell', 'custom')
        """
        if template not in self.JOURNAL_TEMPLATES:
            raise ValueError(f"Unknown template: {template}")

        self.template = template
        self.config = self.JOURNAL_TEMPLATES[template].copy()

    def get_available_templates(self) -> List[str]:
        """
        Get list of available journal templates.

        Returns
        -------
        List[str]
            List of template names
        """
        return list(self.JOURNAL_TEMPLATES.keys())

    def get_template_config(self, template: str) -> Dict[str, Any]:
        """
        Get configuration for a specific template.

        Parameters
        ----------
        template : str
            Template name

        Returns
        -------
        Dict[str, Any]
            Template configuration
        """
        if template not in self.JOURNAL_TEMPLATES:
            raise ValueError(f"Unknown template: {template}")

        return self.JOURNAL_TEMPLATES[template].copy()


def export_panel(
    figure: go.Figure,
    output_path: Union[str, Path],
    template: str = "custom",
    format: str = "png",
    **kwargs,
) -> str:
    """
    Convenience function for exporting a single panel.

    Parameters
    ----------
    figure : plotly.graph_objects.Figure
        Plotly figure to export
    output_path : str or Path
        Output file path
    template : str, default 'custom'
        Journal template to use
    format : str, default 'png'
        Export format
    **kwargs
        Additional arguments passed to export_panel

    Returns
    -------
    str
        Path to the exported file
    """
    exporter = PublicationExporter(template)
    return exporter.export_panel(figure, output_path, format, **kwargs)


def export_batch(
    panels: List[Tuple[go.Figure, str]],
    output_dir: Union[str, Path],
    template: str = "custom",
    format: str = "png",
    **kwargs,
) -> List[str]:
    """
    Convenience function for batch export.

    Parameters
    ----------
    panels : List[Tuple[plotly.graph_objects.Figure, str]]
        List of (figure, panel_name) tuples
    output_dir : str or Path
        Output directory
    template : str, default 'custom'
        Journal template to use
    format : str, default 'png'
        Export format
    **kwargs
        Additional arguments passed to export_batch

    Returns
    -------
    List[str]
        List of exported file paths
    """
    exporter = PublicationExporter(template)
    return exporter.export_batch(panels, output_dir, format, **kwargs)
