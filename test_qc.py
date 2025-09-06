"""
Test script for the QC Metrics Panel implementation.

This script tests the QCPanel functionality with sample data.
"""

import anndata as ad
import numpy as np
import pandas as pd

from pysee import PySEE
from pysee.panels.qc import QCPanel


def create_sample_data():
    """Create sample AnnData for testing QC panel."""
    print("Creating sample data...")

    # Create sample data with realistic QC metrics
    n_cells = 1000
    n_genes = 2000

    # Create expression matrix with some structure
    np.random.seed(42)
    expression_matrix = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes))

    # Add some mitochondrial genes (10% of total)
    n_mito_genes = int(0.1 * n_genes)
    mito_gene_indices = np.random.choice(n_genes, n_mito_genes, replace=False)

    # Make mitochondrial genes more highly expressed
    expression_matrix[:, mito_gene_indices] *= 3

    # Create gene names
    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
    # Make some genes mitochondrial
    for i, idx in enumerate(mito_gene_indices):
        gene_names[idx] = f"MT-Gene_{i:03d}"

    # Create cell names
    cell_names = [f"Cell_{i:04d}" for i in range(n_cells)]

    # Create cell metadata with some QC metrics
    obs_data = {
        "cell_type": np.random.choice(["T_cell", "B_cell", "NK_cell", "Monocyte"], n_cells),
        "batch": np.random.choice(["Batch_1", "Batch_2", "Batch_3"], n_cells),
        "total_counts": expression_matrix.sum(axis=1),
        "detected_genes": (expression_matrix > 0).sum(axis=1),
    }

    # Add mitochondrial percentage
    mito_genes = [name.startswith("MT-") for name in gene_names]
    mito_counts = expression_matrix[:, mito_genes].sum(axis=1)
    total_counts = expression_matrix.sum(axis=1)
    obs_data["mito_percent"] = mito_counts / total_counts * 100

    # Create AnnData
    adata = ad.AnnData(
        X=expression_matrix,
        obs=pd.DataFrame(obs_data, index=cell_names),
        var=pd.DataFrame(index=gene_names),
    )

    print(f"✅ Created sample data: {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"   - Mitochondrial genes: {sum(mito_genes)}")
    print(f"   - Mean mito %: {obs_data['mito_percent'].mean():.1f}%")
    print(f"   - Mean total counts: {obs_data['total_counts'].mean():.0f}")
    print(f"   - Mean detected genes: {obs_data['detected_genes'].mean():.0f}")

    return adata


def test_qc_panel():
    """Test basic QC panel functionality."""
    print("\n🔥 Testing QC Panel Implementation")
    print("=" * 50)

    # Create sample data
    adata = create_sample_data()

    # Initialize PySEE
    print("\nInitializing PySEE dashboard...")
    app = PySEE(adata)

    # Create QC panel
    print("Creating QCPanel...")
    qc_panel = QCPanel("qc", title="Quality Control Metrics")

    # Add panel to dashboard
    app.add_panel("qc", qc_panel)

    # Render the panel
    print("Rendering QC panel...")
    qc_fig = app.render_panel("qc")
    print("✅ QC panel rendered successfully!")
    print(f"   Figure type: {type(qc_fig)}")

    # Test configuration methods
    print("\nTesting configuration methods...")
    qc_panel.set_mito_threshold(15.0)
    qc_panel.set_gene_count_thresholds(500, 10000)
    qc_panel.set_detected_genes_thresholds(100, 3000)
    qc_panel.toggle_metrics(show_mito=True, show_gene_counts=True, show_detected_genes=True)
    print("✅ Configuration methods work!")

    # Test code export
    print("\nTesting code export...")
    code = qc_panel.get_selection_code()
    print("✅ Code export works!")
    print("   Generated code:")
    print(code[:200] + "..." if len(code) > 200 else code)

    return qc_fig


def test_qc_with_other_panels():
    """Test QC panel integration with other panels."""
    print("\n" + "=" * 50)
    print("Testing QC panel integration with other panels...")

    # Create sample data
    adata = create_sample_data()

    # Add some mock UMAP coordinates for testing
    import numpy as np

    np.random.seed(42)
    adata.obsm["X_umap"] = np.random.randn(adata.n_obs, 2)

    # Initialize PySEE
    app = PySEE(adata)

    # Create multiple panels
    from pysee.panels.heatmap import HeatmapPanel
    from pysee.panels.umap import UMAPPanel
    from pysee.panels.violin import ViolinPanel

    umap_panel = UMAPPanel("umap", title="UMAP Visualization")
    violin_panel = ViolinPanel("violin", gene="Gene_0001", title="Gene Expression")
    heatmap_panel = HeatmapPanel("heatmap", title="Gene Expression Heatmap")
    qc_panel = QCPanel("qc", title="Quality Control")

    # Add all panels
    app.add_panel("umap", umap_panel)
    app.add_panel("violin", violin_panel)
    app.add_panel("heatmap", heatmap_panel)
    app.add_panel("qc", qc_panel)

    # Render all panels
    print("Rendering all panels...")
    umap_fig = app.render_panel("umap")
    violin_fig = app.render_panel("violin")
    heatmap_fig = app.render_panel("heatmap")
    qc_fig = app.render_panel("qc")

    print("✅ All panels rendered successfully!")
    print(f"   UMAP: {type(umap_fig)}")
    print(f"   Violin: {type(violin_fig)}")
    print(f"   Heatmap: {type(heatmap_fig)}")
    print(f"   QC: {type(qc_fig)}")

    # Test panel integration
    print("\nTesting panel integration...")
    print(f"   Dashboard has {len(app.panels)} panels")
    print(f"   Panel IDs: {list(app.panels.keys())}")
    print("✅ Panels integrated successfully!")

    return umap_fig, violin_fig, heatmap_fig, qc_fig


def main():
    """Run all QC panel tests."""
    print("🧪 QC Panel Test Suite")
    print("=" * 50)

    try:
        # Test basic QC panel
        test_qc_panel()

        # Test integration with other panels
        test_qc_with_other_panels()

        print("\n" + "=" * 50)
        print("📊 Test Results:")
        print("   Basic QCPanel: ✅ PASS")
        print("   Multi-panel Integration: ✅ PASS")
        print("\n🎉 All tests passed! QCPanel is ready for use.")

    except Exception as e:
        print(f"\n❌ Test failed with error: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
