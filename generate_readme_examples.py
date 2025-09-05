#!/usr/bin/env python3
"""
Generate example images for README documentation.

This script creates sample visualizations to showcase PySEE's
publication export capabilities.
"""

import os
import sys
from pathlib import Path

# Add the current directory to Python path
sys.path.insert(0, str(Path(__file__).parent))

import scanpy as sc
import pandas as pd
import numpy as np
from pysee import PySEE
from pysee.panels import UMAPPanel, ViolinPanel, HeatmapPanel, QCPanel, DotPlotPanel


def create_sample_data():
    """Create sample single-cell data for examples."""
    print("Creating sample data...")

    # Load the PBMC dataset
    adata = sc.datasets.pbmc3k()

    # Basic preprocessing
    adata.var_names_make_unique()
    adata.obs["n_counts"] = adata.X.sum(axis=1)
    adata.obs["n_genes"] = (adata.X > 0).sum(axis=1)

    # Calculate QC metrics
    mt_genes = adata.var_names.str.startswith("MT-")
    adata.obs["pct_counts_mt"] = (
        np.array(adata[:, mt_genes].X.sum(axis=1)).flatten() / adata.obs["n_counts"] * 100
    )

    # Filter cells
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Normalize and log transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Find highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]

    # Scale data
    sc.pp.scale(adata, max_value=10)

    # PCA
    sc.tl.pca(adata, svd_solver="arpack")

    # UMAP
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)

    # Clustering
    sc.tl.leiden(adata, resolution=0.5)

    # Find marker genes
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")

    print(f"Created dataset with {adata.n_obs} cells and {adata.n_vars} genes")
    return adata


def generate_examples():
    """Generate example images for README."""
    print("Generating README examples...")

    # Create sample data
    adata = create_sample_data()

    # Create output directory
    output_dir = Path("examples")
    output_dir.mkdir(exist_ok=True)

    # Create PySEE dashboard
    app = PySEE(adata, title="PySEE Example Dashboard")

    # Add UMAP panel
    umap_panel = UMAPPanel(
        panel_id="umap", embedding="X_umap", color="leiden", title="UMAP Visualization"
    )
    umap_panel.set_config("point_size", 4)
    umap_panel.set_config("opacity", 0.8)
    app.add_panel("umap", umap_panel)

    # Add Violin panel - use a gene that exists in the dataset
    available_genes = adata.var_names.tolist()[:10]  # Get first 10 genes
    violin_gene = available_genes[0]  # Use first available gene

    violin_panel = ViolinPanel(
        panel_id="violin", gene=violin_gene, group_by="leiden", title=f"{violin_gene} Expression"
    )
    violin_panel.set_config("plot_type", "violin")
    violin_panel.set_config("show_points", True)
    app.add_panel("violin", violin_panel)

    # Add QC panel
    qc_panel = QCPanel(panel_id="qc", title="Quality Control")
    app.add_panel("qc", qc_panel)

    # Ensure all panels have data
    for panel in app.panels.values():
        panel.data_wrapper = app.data_wrapper

    # Generate example images
    print("\nGenerating example images...")

    try:
        # Export UMAP panel in different journal styles
        print("  - UMAP Nature style...")
        umap_panel.export(
            output_dir / "umap_nature.png",
            format="png",
            template="nature",
            title="UMAP Visualization",
            caption="Single-cell UMAP embedding colored by Leiden clusters",
        )

        print("  - UMAP Science style...")
        umap_panel.export(
            output_dir / "umap_science.png", format="png", template="science", title="UMAP Plot"
        )

        print("  - UMAP Cell style...")
        umap_panel.export(
            output_dir / "umap_cell.png", format="png", template="cell", title="UMAP Embedding"
        )

        # Export Violin panel
        print("  - Violin plot...")
        violin_panel.export(
            output_dir / "violin_plot.png",
            format="png",
            template="nature",
            title=f"{violin_gene} Expression",
            caption=f"Violin plot showing {violin_gene} expression across cell clusters",
        )

        # Export QC panel
        print("  - QC panel...")
        qc_panel.export(
            output_dir / "qc_panel.png",
            format="png",
            template="nature",
            title="Quality Control Metrics",
            caption="Quality control metrics for single-cell data",
        )

        # Export as SVG (vector format)
        print("  - UMAP SVG...")
        umap_panel.export(output_dir / "umap_vector.svg", format="svg", template="nature")

        # Export as PDF
        print("  - UMAP PDF...")
        umap_panel.export(
            output_dir / "umap_publication.pdf",
            format="pdf",
            template="nature",
            title="UMAP Visualization for Publication",
        )

        # Export dashboard summary
        print("  - Dashboard summary...")
        app.export_dashboard_summary(output_dir / "dashboard_summary.html", template="custom")

        print(f"\n✅ Successfully generated examples in '{output_dir}' directory")

        # List generated files
        print("\n📋 Generated files:")
        for file_path in sorted(output_dir.rglob("*")):
            if file_path.is_file():
                size = file_path.stat().st_size
                print(f"   - {file_path.name} ({size:,} bytes)")

        return output_dir

    except Exception as e:
        print(f"❌ Error generating examples: {e}")
        import traceback

        traceback.print_exc()
        return None


def main():
    """Main function to generate README examples."""
    print("PySEE README Example Generator")
    print("==============================")

    try:
        output_dir = generate_examples()

        if output_dir:
            print(f"\n🎉 Examples generated successfully!")
            print(f"📁 Check the '{output_dir}' directory for example images")
            print("\nThese images can be used in the README to showcase PySEE's capabilities.")
        else:
            print("\n❌ Failed to generate examples")
            return 1

    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback

        traceback.print_exc()
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
