#!/usr/bin/env python3
"""
Test script for PySEE publication-quality export functionality.

This script demonstrates how to use the new export features for creating
publication-ready visualizations.
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

def create_test_data():
    """Create test single-cell data for demonstration."""
    print("Creating test data...")
    
    # Load the PBMC dataset
    adata = sc.datasets.pbmc3k()
    
    # Basic preprocessing
    adata.var_names_make_unique()
    adata.obs['n_counts'] = adata.X.sum(axis=1)
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
    
    # Calculate QC metrics
    mt_genes = adata.var_names.str.startswith('MT-')
    adata.obs['pct_counts_mt'] = np.array(adata[:, mt_genes].X.sum(axis=1)).flatten() / adata.obs['n_counts'] * 100
    
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
    sc.tl.pca(adata, svd_solver='arpack')
    
    # UMAP
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    
    # Clustering
    sc.tl.leiden(adata, resolution=0.5)
    
    # Find marker genes
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    
    print(f"Created dataset with {adata.n_obs} cells and {adata.n_vars} genes")
    return adata

def test_export_functionality():
    """Test the publication export functionality."""
    print("\n" + "="*60)
    print("TESTING PUBLICATION EXPORT FUNCTIONALITY")
    print("="*60)
    
    # Create test data
    adata = create_test_data()
    
    # Create PySEE dashboard
    print("\nCreating PySEE dashboard...")
    app = PySEE(adata, title="Publication Export Test")
    
    # Add panels
    print("Adding visualization panels...")
    app.add_panel(UMAPPanel("umap", color="leiden", title="UMAP Plot"))
    app.add_panel(ViolinPanel("violin", genes=["CD3D", "CD3E", "CD8A"], title="T Cell Markers"))
    app.add_panel(QCPanel("qc", title="Quality Control"))
    
    # Test individual panel export
    print("\n" + "-"*40)
    print("TESTING INDIVIDUAL PANEL EXPORT")
    print("-"*40)
    
    # Create output directory
    output_dir = Path("publication_exports")
    output_dir.mkdir(exist_ok=True)
    
    # Export UMAP panel in different formats
    print("Exporting UMAP panel...")
    try:
        # PNG export
        png_path = app.export_panel(
            "umap", 
            output_dir / "umap_nature.png", 
            format="png", 
            template="nature",
            title="UMAP Visualization",
            caption="Single-cell UMAP embedding colored by Leiden clusters"
        )
        print(f"✅ PNG exported: {png_path}")
        
        # SVG export
        svg_path = app.export_panel(
            "umap", 
            output_dir / "umap_science.svg", 
            format="svg", 
            template="science"
        )
        print(f"✅ SVG exported: {svg_path}")
        
        # PDF export
        pdf_path = app.export_panel(
            "umap", 
            output_dir / "umap_cell.pdf", 
            format="pdf", 
            template="cell"
        )
        print(f"✅ PDF exported: {pdf_path}")
        
    except Exception as e:
        print(f"❌ Error exporting UMAP panel: {e}")
    
    # Test batch export
    print("\n" + "-"*40)
    print("TESTING BATCH EXPORT")
    print("-"*40)
    
    try:
        exported_files = app.export_all_panels(
            output_dir / "batch_export",
            format="png",
            template="nature",
            prefix="figure"
        )
        print(f"✅ Batch exported {len(exported_files)} panels:")
        for file_path in exported_files:
            print(f"   - {file_path}")
    except Exception as e:
        print(f"❌ Error in batch export: {e}")
    
    # Test dashboard summary
    print("\n" + "-"*40)
    print("TESTING DASHBOARD SUMMARY")
    print("-"*40)
    
    try:
        summary_path = app.export_dashboard_summary(
            output_dir / "dashboard_summary.html",
            template="custom"
        )
        print(f"✅ Dashboard summary exported: {summary_path}")
    except Exception as e:
        print(f"❌ Error exporting dashboard summary: {e}")
    
    # Test export info
    print("\n" + "-"*40)
    print("EXPORT CAPABILITIES INFO")
    print("-"*40)
    
    try:
        export_info = app.get_export_info()
        print("Available templates:", export_info["available_templates"])
        print("Supported formats:", export_info["supported_formats"])
        print("Number of panels:", export_info["n_panels"])
        print("Exportable panels:", export_info["exportable_panels"])
        print("Default DPI:", export_info["default_dpi"])
    except Exception as e:
        print(f"❌ Error getting export info: {e}")
    
    # Test panel-level export methods
    print("\n" + "-"*40)
    print("TESTING PANEL-LEVEL EXPORT METHODS")
    print("-"*40)
    
    try:
        umap_panel = app.get_panel("umap")
        
        # Test export as bytes
        png_bytes = umap_panel.export_as_bytes(format="png", template="nature")
        print(f"✅ Panel exported as bytes: {len(png_bytes)} bytes")
        
        # Test export as base64
        base64_str = umap_panel.export_as_base64(format="png", template="nature")
        print(f"✅ Panel exported as base64: {len(base64_str)} characters")
        
        # Test publication info
        pub_info = umap_panel.get_publication_info()
        print(f"✅ Panel publication info: {pub_info['panel_type']} - {pub_info['title']}")
        
    except Exception as e:
        print(f"❌ Error testing panel-level export: {e}")
    
    print("\n" + "="*60)
    print("EXPORT TESTING COMPLETE")
    print("="*60)
    print(f"Check the '{output_dir}' directory for exported files")
    
    return output_dir

def main():
    """Main function to run the export tests."""
    print("PySEE Publication Export Test")
    print("=============================")
    
    try:
        output_dir = test_export_functionality()
        
        print(f"\n🎉 All tests completed successfully!")
        print(f"📁 Exported files are in: {output_dir.absolute()}")
        
        # List exported files
        if output_dir.exists():
            print("\n📋 Exported files:")
            for file_path in sorted(output_dir.rglob("*")):
                if file_path.is_file():
                    size = file_path.stat().st_size
                    print(f"   - {file_path.name} ({size:,} bytes)")
        
    except Exception as e:
        print(f"\n❌ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
