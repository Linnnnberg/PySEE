"""
Demo script to showcase QC Panel and other PySEE visualizations.

This script creates interactive HTML files showing the QC panel and other panels.
"""

import numpy as np
import pandas as pd
import anndata as ad
from pysee import PySEE
from pysee.panels.qc import QCPanel
from pysee.panels.umap import UMAPPanel
from pysee.panels.violin import ViolinPanel
from pysee.panels.heatmap import HeatmapPanel


def create_realistic_sample_data():
    """Create realistic sample data with proper QC metrics."""
    print("Creating realistic sample data...")

    n_cells = 2000
    n_genes = 3000

    # Create expression matrix with realistic structure
    np.random.seed(42)

    # Base expression with some structure
    expression_matrix = np.random.negative_binomial(3, 0.4, size=(n_cells, n_genes)).astype(float)

    # Add some highly variable genes
    n_hvg = 200
    hvg_indices = np.random.choice(n_genes, n_hvg, replace=False)
    expression_matrix[:, hvg_indices] *= np.random.lognormal(0, 1, n_hvg)

    # Create gene names
    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]

    # Add mitochondrial genes (15% of total)
    n_mito_genes = int(0.15 * n_genes)
    mito_gene_indices = np.random.choice(n_genes, n_mito_genes, replace=False)
    for i, idx in enumerate(mito_gene_indices):
        gene_names[idx] = f"MT-Gene_{i:03d}"

    # Make mitochondrial genes more highly expressed
    expression_matrix[:, mito_gene_indices] *= 2.5

    # Create cell names
    cell_names = [f"Cell_{i:04d}" for i in range(n_cells)]

    # Create cell metadata with realistic QC metrics
    obs_data = {
        "cell_type": np.random.choice(
            ["T_cell", "B_cell", "NK_cell", "Monocyte", "Dendritic"],
            n_cells,
            p=[0.4, 0.25, 0.15, 0.15, 0.05],
        ),
        "batch": np.random.choice(["Batch_1", "Batch_2", "Batch_3", "Batch_4"], n_cells),
        "total_counts": expression_matrix.sum(axis=1),
        "detected_genes": (expression_matrix > 0).sum(axis=1),
    }

    # Add realistic mitochondrial percentage
    mito_genes = [name.startswith("MT-") for name in gene_names]
    mito_counts = expression_matrix[:, mito_genes].sum(axis=1)
    total_counts = expression_matrix.sum(axis=1)
    obs_data["mito_percent"] = mito_counts / total_counts * 100

    # Add some cells with high mitochondrial content (dying cells)
    high_mito_cells = np.random.choice(n_cells, int(0.1 * n_cells), replace=False)
    obs_data["mito_percent"][high_mito_cells] += np.random.exponential(15, len(high_mito_cells))

    # Create AnnData
    adata = ad.AnnData(
        X=expression_matrix.astype(float),
        obs=pd.DataFrame(obs_data, index=cell_names),
        var=pd.DataFrame(index=gene_names),
    )

    # Add UMAP coordinates
    np.random.seed(42)
    adata.obsm["X_umap"] = np.random.randn(adata.n_obs, 2)

    print(f"✅ Created realistic sample data: {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"   - Mitochondrial genes: {sum(mito_genes)}")
    print(f"   - Mean mito %: {obs_data['mito_percent'].mean():.1f}%")
    print(f"   - Mean total counts: {obs_data['total_counts'].mean():.0f}")
    print(f"   - Mean detected genes: {obs_data['detected_genes'].mean():.0f}")
    print(f"   - Cell types: {list(adata.obs['cell_type'].unique())}")

    return adata


def demo_qc_panel():
    """Demo the QC panel with different configurations."""
    print("\n🔥 QC Panel Demo")
    print("=" * 50)

    # Create data
    adata = create_realistic_sample_data()

    # Initialize PySEE
    app = PySEE(adata)

    # Create QC panel with default settings
    qc_panel = QCPanel("qc", title="Quality Control Metrics")
    app.add_panel("qc", qc_panel)

    # Render and save
    qc_fig = app.render_panel("qc")
    qc_fig.write_html("qc_panel_demo.html")
    print("✅ QC Panel saved to: qc_panel_demo.html")

    # Demo with custom thresholds
    print("\nCreating QC panel with custom thresholds...")
    qc_panel.set_mito_threshold(25.0)
    qc_panel.set_gene_count_thresholds(500, 15000)
    qc_panel.set_detected_genes_thresholds(100, 2500)

    qc_fig_custom = app.render_panel("qc")
    qc_fig_custom.write_html("qc_panel_custom_thresholds.html")
    print("✅ QC Panel with custom thresholds saved to: qc_panel_custom_thresholds.html")

    return qc_fig


def demo_multi_panel_dashboard():
    """Demo a complete multi-panel dashboard."""
    print("\n🔥 Multi-Panel Dashboard Demo")
    print("=" * 50)

    # Create data
    adata = create_realistic_sample_data()

    # Initialize PySEE
    app = PySEE(adata, title="PySEE Multi-Panel Dashboard")

    # Create all panels
    umap_panel = UMAPPanel("umap", title="UMAP Visualization", color="cell_type")
    violin_panel = ViolinPanel(
        "violin", gene="Gene_0001", group_by="cell_type", title="Gene Expression by Cell Type"
    )
    heatmap_panel = HeatmapPanel("heatmap", title="Gene Expression Heatmap")
    qc_panel = QCPanel("qc", title="Quality Control Metrics")

    # Add all panels
    app.add_panel("umap", umap_panel)
    app.add_panel("violin", violin_panel)
    app.add_panel("heatmap", heatmap_panel)
    app.add_panel("qc", qc_panel)

    # Render each panel
    print("Rendering UMAP panel...")
    umap_fig = app.render_panel("umap")
    umap_fig.write_html("umap_panel_demo.html")
    print("✅ UMAP Panel saved to: umap_panel_demo.html")

    print("Rendering Violin panel...")
    violin_fig = app.render_panel("violin")
    violin_fig.write_html("violin_panel_demo.html")
    print("✅ Violin Panel saved to: violin_panel_demo.html")

    print("Rendering Heatmap panel...")
    heatmap_fig = app.render_panel("heatmap")
    heatmap_fig.write_html("heatmap_panel_demo.html")
    print("✅ Heatmap Panel saved to: heatmap_panel_demo.html")

    print("Rendering QC panel...")
    qc_fig = app.render_panel("qc")
    qc_fig.write_html("qc_panel_multi_demo.html")
    print("✅ QC Panel (multi-demo) saved to: qc_panel_multi_demo.html")

    return umap_fig, violin_fig, heatmap_fig, qc_fig


def demo_qc_filtering():
    """Demo QC-based filtering workflow."""
    print("\n🔥 QC Filtering Workflow Demo")
    print("=" * 50)

    # Create data
    adata = create_realistic_sample_data()

    # Initialize PySEE
    app = PySEE(adata)

    # Create QC panel with filtering thresholds
    qc_panel = QCPanel("qc", title="QC Filtering Workflow")
    qc_panel.set_mito_threshold(20.0)
    qc_panel.set_gene_count_thresholds(1000, 20000)
    qc_panel.set_detected_genes_thresholds(200, 3000)

    app.add_panel("qc", qc_panel)

    # Render QC panel
    qc_fig = app.render_panel("qc")
    qc_fig.write_html("qc_filtering_workflow.html")
    print("✅ QC Filtering Workflow saved to: qc_filtering_workflow.html")

    # Show the generated filtering code
    print("\nGenerated QC filtering code:")
    print("-" * 40)
    code = qc_panel.get_selection_code()
    print(code)

    return qc_fig


def main():
    """Run all QC panel demos."""
    print("🎨 PySEE QC Panel Chart Demo")
    print("=" * 60)

    try:
        # Demo 1: Basic QC Panel
        demo_qc_panel()

        # Demo 2: Multi-Panel Dashboard
        demo_multi_panel_dashboard()

        # Demo 3: QC Filtering Workflow
        demo_qc_filtering()

        print("\n" + "=" * 60)
        print("📊 Demo Results:")
        print("✅ QC Panel Demo: qc_panel_demo.html")
        print("✅ QC Panel (Custom): qc_panel_custom_thresholds.html")
        print("✅ UMAP Panel: umap_panel_demo.html")
        print("✅ Violin Panel: violin_panel_demo.html")
        print("✅ Heatmap Panel: heatmap_panel_demo.html")
        print("✅ QC Panel (Multi): qc_panel_multi_demo.html")
        print("✅ QC Filtering: qc_filtering_workflow.html")
        print("\n🎉 All demos completed! Open the HTML files in your browser to view the charts.")

    except Exception as e:
        print(f"\n❌ Demo failed with error: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
