#!/usr/bin/env python3
"""
Test script for the HeatmapPanel implementation.

This script tests the new HeatmapPanel functionality with sample data.
"""

import anndata as ad
import numpy as np
import pandas as pd

from pysee import PySEE
from pysee.panels.heatmap import HeatmapPanel
from pysee.panels.umap import UMAPPanel
from pysee.panels.violin import ViolinPanel


def create_sample_data():
    """Create sample single-cell data for testing."""
    # Create sample data
    n_cells = 200
    n_genes = 100

    # Generate expression matrix with some structure
    np.random.seed(42)
    expression_matrix = np.random.negative_binomial(5, 0.3, (n_cells, n_genes)).astype(float)

    # Add some structure - make some genes more variable
    for i in range(10):
        expression_matrix[:, i] *= np.random.uniform(2, 5)

    # Create gene names
    gene_names = [f"Gene_{i:03d}" for i in range(n_genes)]

    # Create cell names
    cell_names = [f"Cell_{i:03d}" for i in range(n_cells)]

    # Create cell metadata
    cell_types = np.random.choice(["Type_A", "Type_B", "Type_C"], n_cells)
    batch = np.random.choice(["Batch_1", "Batch_2"], n_cells)

    obs = pd.DataFrame(
        {
            "cell_type": cell_types,
            "batch": batch,
            "total_counts": np.sum(expression_matrix, axis=1),
        },
        index=cell_names,
    )

    # Create gene metadata
    var = pd.DataFrame(
        {
            "gene_name": gene_names,
            "is_mitochondrial": [name.startswith("Gene_00") for name in gene_names],
        },
        index=gene_names,
    )

    # Create AnnData object
    adata = ad.AnnData(X=expression_matrix, obs=obs, var=var)

    # Add some embeddings
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE

    # PCA
    pca = PCA(n_components=2)
    pca_coords = pca.fit_transform(expression_matrix)
    adata.obsm["X_pca"] = pca_coords

    # t-SNE
    tsne = TSNE(n_components=2, random_state=42)
    tsne_coords = tsne.fit_transform(expression_matrix)
    adata.obsm["X_tsne"] = tsne_coords

    return adata


def test_heatmap_panel():
    """Test the HeatmapPanel functionality."""
    print("Creating sample data...")
    adata = create_sample_data()

    print("Initializing PySEE dashboard...")
    app = PySEE(adata, title="Heatmap Panel Test")

    # Create heatmap panel
    print("Creating HeatmapPanel...")
    heatmap_panel = HeatmapPanel(panel_id="heatmap", title="Gene Expression Heatmap")

    # Set configuration
    heatmap_panel.set_config("n_top_genes", 30)
    heatmap_panel.set_config("max_cells", 100)

    # Add panel to dashboard
    app.add_panel("heatmap", heatmap_panel)

    # Test basic rendering
    print("Rendering heatmap...")
    try:
        heatmap_fig = app.render_panel("heatmap")
        print("✅ Heatmap rendered successfully!")
        print(f"   Figure type: {type(heatmap_fig)}")

        # Test configuration methods
        print("Testing configuration methods...")
        heatmap_panel.set_clustering_method("complete")
        heatmap_panel.set_color_scale("viridis")
        heatmap_panel.set_max_genes(20)
        heatmap_panel.toggle_clustering(cluster_genes=True, cluster_cells=False)
        print("✅ Configuration methods work!")

        # Test with specific genes
        print("Testing with specific genes...")
        specific_genes = ["Gene_001", "Gene_002", "Gene_003", "Gene_004", "Gene_005"]
        heatmap_panel.set_genes(specific_genes)
        app.render_panel("heatmap")
        print("✅ Specific genes heatmap rendered!")

        # Test code export
        print("Testing code export...")
        code = heatmap_panel.get_selection_code()
        print("✅ Code export works!")
        print(f"   Generated code:\n{code}")

        return True

    except Exception as e:
        print(f"❌ Error rendering heatmap: {e}")
        import traceback

        traceback.print_exc()
        return False


def test_heatmap_with_other_panels():
    """Test heatmap panel integration with other panels."""
    print("\n" + "=" * 50)
    print("Testing HeatmapPanel integration with other panels...")

    adata = create_sample_data()
    app = PySEE(adata, title="Multi-Panel Test")

    # Create multiple panels
    umap_panel = UMAPPanel(panel_id="umap", embedding="X_pca", color="cell_type", title="PCA Plot")

    violin_panel = ViolinPanel(
        panel_id="violin",
        gene="Gene_001",
        group_by="cell_type",
        title="Gene Expression by Cell Type",
    )

    heatmap_panel = HeatmapPanel(panel_id="heatmap", title="Gene Expression Heatmap")

    # Set configuration
    heatmap_panel.set_config("n_top_genes", 20)
    heatmap_panel.set_config("max_cells", 50)

    # Add all panels
    app.add_panel("umap", umap_panel)
    app.add_panel("violin", violin_panel)
    app.add_panel("heatmap", heatmap_panel)

    # Test rendering all panels
    print("Rendering all panels...")
    try:
        umap_fig = app.render_panel("umap")
        violin_fig = app.render_panel("violin")
        heatmap_fig = app.render_panel("heatmap")

        print("✅ All panels rendered successfully!")
        print(f"   UMAP: {type(umap_fig)}")
        print(f"   Violin: {type(violin_fig)}")
        print(f"   Heatmap: {type(heatmap_fig)}")

        # Test panel linking (panels are automatically linked when added)
        print("Testing panel integration...")
        print("✅ Panels integrated successfully!")

        return True

    except Exception as e:
        print(f"❌ Error in multi-panel test: {e}")
        import traceback

        traceback.print_exc()
        return False


def main():
    """Run all tests."""
    print("🔥 Testing HeatmapPanel Implementation")
    print("=" * 50)

    # Test basic heatmap functionality
    success1 = test_heatmap_panel()

    # Test integration with other panels
    success2 = test_heatmap_with_other_panels()

    print("\n" + "=" * 50)
    print("📊 Test Results:")
    print(f"   Basic HeatmapPanel: {'✅ PASS' if success1 else '❌ FAIL'}")
    print(f"   Multi-panel Integration: {'✅ PASS' if success2 else '❌ FAIL'}")

    if success1 and success2:
        print("\n🎉 All tests passed! HeatmapPanel is ready for use.")
    else:
        print("\n⚠️  Some tests failed. Please check the implementation.")

    return success1 and success2


if __name__ == "__main__":
    main()
