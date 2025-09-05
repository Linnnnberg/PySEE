"""
PySEE Example - Interactive Bioinformatics Visualization

This example demonstrates how to use PySEE to create interactive
visualizations for single-cell data analysis.
"""

import numpy as np
import pandas as pd
import scanpy as sc

from pysee import PySEE, UMAPPanel, ViolinPanel

# Load sample data (you can replace this with your own data)
print("Loading sample data...")
adata = sc.datasets.pbmc3k()

# Basic preprocessing
print("Preprocessing data...")
adata.var_names_make_unique()
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.raw = adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

print(f"Data shape: {adata.shape}")
print(f"Available embeddings: {list(adata.obsm.keys())}")
print(f"Available obs columns: {list(adata.obs.columns)}")

# Create PySEE dashboard
print("\nCreating PySEE dashboard...")
app = PySEE(adata, title="PBMC3K Analysis")

# Add UMAP panel
print("Adding UMAP panel...")
app.add_panel(
    "umap", UMAPPanel(panel_id="umap", embedding="X_umap", color="leiden", title="UMAP Plot")
)

# Add violin panel for a specific gene
print("Adding Violin panel...")
app.add_panel(
    "violin",
    ViolinPanel(
        panel_id="violin", gene="CD3D", group_by="leiden", title="Gene Expression"
    ),  # T-cell marker
)

# Link panels so that UMAP selection affects violin plot
print("Linking panels...")
app.link("umap", "violin")

# Display dashboard information
print("\n" + "=" * 50)
print("PySEE Dashboard Information")
print("=" * 50)
app.show()

# Test panel rendering
print("\n" + "=" * 50)
print("Testing Panel Rendering")
print("=" * 50)

try:
    umap_fig = app.render_panel("umap")
    print("✓ UMAP panel rendered successfully")
    print(f"  Figure type: {type(umap_fig)}")
except Exception as e:
    print(f"✗ UMAP panel rendering failed: {e}")

try:
    violin_fig = app.render_panel("violin")
    print("✓ Violin panel rendered successfully")
    print(f"  Figure type: {type(violin_fig)}")
except Exception as e:
    print(f"✗ Violin panel rendering failed: {e}")

# Test selection functionality
print("\n" + "=" * 50)
print("Testing Selection Functionality")
print("=" * 50)

# Create a random selection
selection = np.random.choice([True, False], size=adata.n_obs, p=[0.2, 0.8])
app.set_global_selection(selection)
n_selected = np.sum(selection)
print(f"✓ Set global selection: {n_selected} cells selected")

# Test code export
print("\n" + "=" * 50)
print("Testing Code Export")
print("=" * 50)

try:
    code = app.export_code()
    print("✓ Code export successful")
    print(f"  Exported code length: {len(code)} characters")

    # Save code to file
    with open("exported_code.py", "w") as f:
        f.write(code)
    print("  Code saved to 'exported_code.py'")

except Exception as e:
    print(f"✗ Code export failed: {e}")

# Display dashboard summary
print("\n" + "=" * 50)
print("Dashboard Summary")
print("=" * 50)
info = app.get_dashboard_info()
print(f"Title: {info['title']}")
print(f"Number of panels: {info['n_panels']}")
print(f"Panel order: {info['panel_order']}")
print(f"Has global selection: {info['has_global_selection']}")
if info["has_global_selection"]:
    print(f"Selected cells: {info['n_selected_cells']}")

print("\n" + "=" * 50)
print("Example completed successfully!")
print("=" * 50)
print("\nTo use this in a Jupyter notebook:")
print("1. Run this script to generate the dashboard")
print("2. Use app.render_panel('umap') to get the UMAP plot")
print("3. Use app.render_panel('violin') to get the violin plot")
print("4. Use app.export_code() to get reproducible Python code")
