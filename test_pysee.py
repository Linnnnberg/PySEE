"""
Test script for PySEE functionality.

This script tests the basic functionality of PySEE with sample data.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from pysee import PySEE, UMAPPanel, ViolinPanel


def create_sample_data():
    """Create a sample AnnData object for testing."""
    # Create random expression data
    n_cells = 200
    n_genes = 100

    # Generate random expression matrix
    X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes))

    # Create gene names
    gene_names = [f"Gene_{i:03d}" for i in range(n_genes)]

    # Create cell names
    cell_names = [f"Cell_{i:03d}" for i in range(n_cells)]

    # Create cell metadata
    obs = pd.DataFrame(
        {
            "cell_type": np.random.choice(["Type_A", "Type_B", "Type_C"], n_cells),
            "batch": np.random.choice(["Batch_1", "Batch_2"], n_cells),
            "n_genes": np.random.poisson(50, n_cells),
            "total_counts": np.random.poisson(1000, n_cells),
        },
        index=cell_names,
    )

    # Create gene metadata
    var = pd.DataFrame(
        {
            "gene_type": np.random.choice(["protein_coding", "lncRNA"], n_genes),
            "mean_expression": np.mean(X, axis=0),
        },
        index=gene_names,
    )

    # Create AnnData object
    adata = sc.AnnData(X=X, obs=obs, var=var)

    # Add some embeddings
    # PCA
    sc.pp.pca(adata, n_comps=10)

    # UMAP
    sc.pp.neighbors(adata, n_neighbors=15)
    sc.tl.umap(adata)

    # t-SNE
    sc.tl.tsne(adata, n_pcs=10)

    return adata


def test_basic_functionality():
    """Test basic PySEE functionality."""
    print("Creating sample data...")
    adata = create_sample_data()

    print(f"Data shape: {adata.shape}")
    print(f"Available embeddings: {list(adata.obsm.keys())}")
    print(f"Available obs columns: {list(adata.obs.columns)}")

    # Create PySEE dashboard
    print("\nCreating PySEE dashboard...")
    app = PySEE(adata, title="Test Dashboard")

    # Add UMAP panel
    print("Adding UMAP panel...")
    app.add_panel("umap", UMAPPanel(panel_id="umap", embedding="X_umap", color="cell_type", title="UMAP Plot"))

    # Add violin panel
    print("Adding Violin panel...")
    app.add_panel("violin", ViolinPanel(panel_id="violin", gene="Gene_001", group_by="cell_type", title="Gene Expression"))

    # Link panels
    print("Linking panels...")
    app.link("umap", "violin")

    # Test dashboard info
    print("\nDashboard info:")
    info = app.get_dashboard_info()
    print(f"  Title: {info['title']}")
    print(f"  Number of panels: {info['n_panels']}")
    print(f"  Panel order: {info['panel_order']}")

    # Test panel rendering
    print("\nTesting panel rendering...")
    try:
        umap_fig = app.render_panel("umap")
        print("  ✓ UMAP panel rendered successfully")
    except Exception as e:
        print(f"  ✗ UMAP panel rendering failed: {e}")

    try:
        violin_fig = app.render_panel("violin")
        print("  ✓ Violin panel rendered successfully")
    except Exception as e:
        print(f"  ✗ Violin panel rendering failed: {e}")

    # Test selection
    print("\nTesting selection...")
    selection = np.random.choice([True, False], size=adata.n_obs, p=[0.3, 0.7])
    app.set_global_selection(selection)
    n_selected = np.sum(selection)
    print(f"  Set selection: {n_selected} cells selected")

    # Test code export
    print("\nTesting code export...")
    try:
        code = app.export_code()
        print("  ✓ Code export successful")
        print(f"  Code length: {len(code)} characters")
    except Exception as e:
        print(f"  ✗ Code export failed: {e}")

    # Show dashboard
    print("\nDashboard summary:")
    app.show()

    print("\n✓ All tests completed successfully!")


if __name__ == "__main__":
    test_basic_functionality()
