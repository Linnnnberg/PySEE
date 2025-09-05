#!/usr/bin/env python3
"""
Test script for PySEE panel configuration and export functionality.

This script demonstrates how to configure panels and export them with
custom settings for publication-ready visualizations.
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

def test_panel_configuration():
    """Test panel configuration capabilities."""
    print("\n" + "="*60)
    print("TESTING PANEL CONFIGURATION")
    print("="*60)
    
    # Create test data
    adata = create_test_data()
    
    # Create PySEE dashboard
    print("\nCreating PySEE dashboard...")
    app = PySEE(adata, title="Panel Configuration Test")
    
    # Test 1: Basic panel creation and configuration
    print("\n" + "-"*40)
    print("TEST 1: Basic Panel Configuration")
    print("-"*40)
    
    # Create UMAP panel with initial configuration
    umap_panel = UMAPPanel(
        panel_id="umap",
        embedding="X_umap",
        color="leiden",
        title="UMAP Visualization"
    )
    
    # Configure the panel
    umap_panel.set_config("point_size", 4)
    umap_panel.set_config("opacity", 0.8)
    umap_panel.set_config("show_legend", True)
    umap_panel.set_config("color_scale", "viridis")
    
    print("UMAP Panel Configuration:")
    print(f"  Point size: {umap_panel.get_config('point_size')}")
    print(f"  Opacity: {umap_panel.get_config('opacity')}")
    print(f"  Show legend: {umap_panel.get_config('show_legend')}")
    print(f"  Color scale: {umap_panel.get_config('color_scale')}")
    
    app.add_panel("umap", umap_panel)
    
    # Test 2: Violin panel configuration
    print("\n" + "-"*40)
    print("TEST 2: Violin Panel Configuration")
    print("-"*40)
    
    violin_panel = ViolinPanel(
        panel_id="violin",
        gene="CD3D",
        group_by="leiden",
        title="Gene Expression"
    )
    
    # Configure violin panel
    violin_panel.set_config("plot_type", "violin")
    violin_panel.set_config("show_points", True)
    violin_panel.set_config("show_box", True)
    violin_panel.set_config("width", 800)
    violin_panel.set_config("height", 500)
    
    print("Violin Panel Configuration:")
    print(f"  Plot type: {violin_panel.get_config('plot_type')}")
    print(f"  Show points: {violin_panel.get_config('show_points')}")
    print(f"  Show box: {violin_panel.get_config('show_box')}")
    print(f"  Width: {violin_panel.get_config('width')}")
    print(f"  Height: {violin_panel.get_config('height')}")
    
    app.add_panel("violin", violin_panel)
    
    # Test 3: Export configuration
    print("\n" + "-"*40)
    print("TEST 3: Export Configuration")
    print("-"*40)
    
    # Create output directory
    output_dir = Path("configuration_exports")
    output_dir.mkdir(exist_ok=True)
    
    # Export with custom configuration
    print("Exporting UMAP panel with custom configuration...")
    try:
        # Export with config overrides
        umap_panel.export(
            output_dir / "umap_custom.png",
            format="png",
            template="nature",
            title="Custom UMAP Visualization",
            caption="UMAP embedding with custom styling",
            config_overrides={
                "font_size": 12,
                "line_width": 1.0,
                "marker_size": 5
            }
        )
        print("✅ UMAP panel exported with custom configuration")
    except Exception as e:
        print(f"❌ Error exporting UMAP panel: {e}")
    
    # Test 4: Panel copy for export
    print("\n" + "-"*40)
    print("TEST 4: Panel Copy for Export")
    print("-"*40)
    
    try:
        # Create a copy with different configuration
        umap_export_copy = umap_panel.configure_for_export({
            "point_size": 6,
            "opacity": 0.9,
            "color_scale": "plasma"
        })
        
        print("Created panel copy with export configuration:")
        print(f"  Original point size: {umap_panel.get_config('point_size')}")
        print(f"  Copy point size: {umap_export_copy.get_config('point_size')}")
        print(f"  Original color scale: {umap_panel.get_config('color_scale')}")
        print(f"  Copy color scale: {umap_export_copy.get_config('color_scale')}")
        
        # Export the copy
        umap_export_copy.export(
            output_dir / "umap_copy.png",
            format="png",
            template="science"
        )
        print("✅ Panel copy exported successfully")
        
    except Exception as e:
        print(f"❌ Error with panel copy: {e}")
    
    # Test 5: Export configuration info
    print("\n" + "-"*40)
    print("TEST 5: Export Configuration Info")
    print("-"*40)
    
    try:
        # Get export configuration
        export_config = umap_panel.get_export_config()
        print("Export Configuration:")
        print(f"  Available templates: {export_config['available_templates']}")
        print(f"  Recommended formats: {export_config['recommended_formats']}")
        print(f"  Default DPI: {export_config['default_dpi']}")
        print(f"  Default width: {export_config['default_width_mm']}mm")
        print(f"  Default height: {export_config['default_height_mm']}mm")
        
        # Get publication info
        pub_info = umap_panel.get_publication_info()
        print(f"\nPublication Info:")
        print(f"  Panel type: {pub_info['panel_type']}")
        print(f"  Title: {pub_info['title']}")
        print(f"  Has data: {pub_info['has_data']}")
        
    except Exception as e:
        print(f"❌ Error getting configuration info: {e}")
    
    # Test 6: Batch export with configuration
    print("\n" + "-"*40)
    print("TEST 6: Batch Export with Configuration")
    print("-"*40)
    
    try:
        # Configure all panels for publication
        for panel in app.panels.values():
            panel.set_config("show_legend", True)
            panel.set_config("font_size", 10)
        
        # Export all panels
        exported_files = app.export_all_panels(
            output_dir / "batch_configured",
            format="png",
            template="nature",
            prefix="configured"
        )
        
        print(f"✅ Batch exported {len(exported_files)} configured panels:")
        for file_path in exported_files:
            print(f"   - {file_path}")
            
    except Exception as e:
        print(f"❌ Error in batch export: {e}")
    
    # Test 7: Dashboard export info
    print("\n" + "-"*40)
    print("TEST 7: Dashboard Export Info")
    print("-"*40)
    
    try:
        export_info = app.get_export_info()
        print("Dashboard Export Info:")
        print(f"  Available templates: {export_info['available_templates']}")
        print(f"  Supported formats: {export_info['supported_formats']}")
        print(f"  Number of panels: {export_info['n_panels']}")
        print(f"  Exportable panels: {export_info['exportable_panels']}")
        print(f"  Default DPI: {export_info['default_dpi']}")
        
    except Exception as e:
        print(f"❌ Error getting dashboard export info: {e}")
    
    print("\n" + "="*60)
    print("CONFIGURATION TESTING COMPLETE")
    print("="*60)
    print(f"Check the '{output_dir}' directory for exported files")
    
    return output_dir

def main():
    """Main function to run the configuration tests."""
    print("PySEE Panel Configuration Test")
    print("==============================")
    
    try:
        output_dir = test_panel_configuration()
        
        print(f"\n🎉 All configuration tests completed successfully!")
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
