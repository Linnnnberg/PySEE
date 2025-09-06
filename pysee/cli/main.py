"""
Command-line interface for PySEE.

This module provides a simple CLI for running PySEE dashboards.
"""

import argparse
import sys
from pathlib import Path

import scanpy as sc

from ..core.dashboard import PySEE
from ..panels.umap import UMAPPanel
from ..panels.violin import ViolinPanel


def main() -> None:
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="PySEE - Interactive, Reproducible Bioinformatics Visualization for Python"
    )

    parser.add_argument("data_file", help="Path to the data file (h5ad, h5, or csv)")

    parser.add_argument(
        "--title",
        default="PySEE Dashboard",
        help="Title for the dashboard (default: 'PySEE Dashboard')",
    )

    parser.add_argument(
        "--umap-embedding", default="X_umap", help="Embedding key for UMAP plot (default: 'X_umap')"
    )

    parser.add_argument("--umap-color", help="Column name to use for coloring UMAP points")

    parser.add_argument("--violin-gene", help="Gene name for violin plot")

    parser.add_argument("--violin-group", help="Column name to use for grouping violin plot")

    parser.add_argument(
        "--export-code", action="store_true", help="Export Python code instead of running dashboard"
    )

    args = parser.parse_args()

    # Load data
    data_path = Path(args.data_file)
    if not data_path.exists():
        print(f"Error: Data file '{data_path}' not found", file=sys.stderr)
        sys.exit(1)

    try:
        if data_path.suffix == ".h5ad":
            adata = sc.read_h5ad(data_path)
        elif data_path.suffix == ".h5":
            adata = sc.read_h5ad(data_path)
        else:
            print(f"Error: Unsupported file format '{data_path.suffix}'", file=sys.stderr)
            sys.exit(1)
    except Exception as e:
        print(f"Error loading data: {e}", file=sys.stderr)
        sys.exit(1)

    # Create dashboard
    app = PySEE(adata, title=args.title)

    # Add UMAP panel
    app.add_panel(
        "umap",
        UMAPPanel(
            panel_id="umap", embedding=args.umap_embedding, color=args.umap_color, title="UMAP Plot"
        ),
    )

    # Add violin panel if gene is specified
    if args.violin_gene:
        app.add_panel(
            "violin",
            ViolinPanel(
                panel_id="violin",
                gene=args.violin_gene,
                group_by=args.violin_group,
                title="Gene Expression",
            ),
        )

        # Link panels
        app.link("umap", "violin")

    # Export code or show dashboard
    if args.export_code:
        print(app.export_code())
    else:
        app.show()


if __name__ == "__main__":
    main()
