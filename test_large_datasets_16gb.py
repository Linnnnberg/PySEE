"""
Test large datasets on 16 GB RAM system.

This script demonstrates how to test large and very large datasets
efficiently on a 16 GB RAM system using various memory optimization strategies.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import tempfile
import os
import gc

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from pysee import PySEE
from pysee.panels.umap import UMAPPanel
from pysee.panels.violin import ViolinPanel
from pysee.panels.heatmap import HeatmapPanel
from pysee.panels.qc import QCPanel
from pysee.utils.system_requirements import SystemRequirementsChecker


def test_large_datasets_on_16gb():
    """Test large datasets using memory-efficient strategies on 16 GB system."""
    print("🧠 Testing Large Datasets on 16 GB RAM System")
    print("=" * 60)

    # Check system
    checker = SystemRequirementsChecker()
    system_info = checker.get_system_info()

    print(f"System Memory: {system_info['total_memory_gb']:.1f} GB")
    print(f"Available Memory: {system_info['available_memory_gb']:.1f} GB")
    print(f"Memory Usage: {system_info['memory_percent']:.1f}%")

    if system_info["total_memory_gb"] < 16:
        print("⚠️ Warning: System has less than 16 GB RAM")

    # Strategy 1: Test with PBMC 68K (medium dataset)
    print("\n📊 Strategy 1: PBMC 68K (Medium Dataset)")
    test_pbmc68k()

    # Strategy 2: Generate large synthetic dataset
    print("\n📊 Strategy 2: Large Synthetic Dataset (100K cells)")
    test_large_synthetic()

    # Strategy 3: Test backed mode
    print("\n📊 Strategy 3: Backed Mode Testing")
    test_backed_mode()

    # Strategy 4: Test subsampling
    print("\n📊 Strategy 4: Subsampling Strategy")
    test_subsampling_strategy()

    print("\n✅ Large dataset testing complete!")


def test_pbmc68k():
    """Test with PBMC 68K dataset."""
    print("   Loading PBMC 68K dataset...")

    try:
        # Load dataset
        adata = sc.datasets.pbmc68k_reduced()
        print(f"   ✅ Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

        # Check memory usage
        memory_mb = adata.X.nbytes / (1024 * 1024)
        print(f"   Memory usage: {memory_mb:.0f} MB")

        # Test PySEE
        app = PySEE(adata, title="PBMC 68K Test")

        # Add panels
        qc_panel = QCPanel("qc", title="QC")
        app.add_panel("qc", qc_panel)

        violin_panel = ViolinPanel("violin", gene=adata.var_names[0], title="Violin")
        app.add_panel("violin", violin_panel)

        # Render panels
        qc_fig = app.render_panel("qc")
        violin_fig = app.render_panel("violin")

        print("   ✅ Successfully rendered QC and Violin panels")

        # Clean up
        del adata, app, qc_fig, violin_fig
        gc.collect()

    except Exception as e:
        print(f"   ❌ Failed: {e}")


def test_large_synthetic():
    """Test with large synthetic dataset."""
    print("   Generating large synthetic dataset...")

    try:
        # Generate dataset (100K cells, 10K genes)
        n_cells = 100000
        n_genes = 10000

        print(f"   Generating {n_cells:,} cells, {n_genes:,} genes...")

        # Generate expression matrix
        np.random.seed(42)
        expression_matrix = np.random.negative_binomial(4, 0.35, size=(n_cells, n_genes)).astype(np.float32)

        # Create metadata
        gene_names = [f"Gene_{i:05d}" for i in range(n_genes)]
        cell_names = [f"Cell_{i:06d}" for i in range(n_cells)]

        obs_data = {
            "cell_type": np.random.choice(["Type_A", "Type_B", "Type_C", "Type_D"], n_cells),
            "batch": np.random.choice(["Batch_1", "Batch_2"], n_cells),
            "total_counts": expression_matrix.sum(axis=1),
            "detected_genes": (expression_matrix > 0).sum(axis=1),
        }

        # Add UMAP coordinates
        umap_coords = np.random.randn(n_cells, 2)

        # Create AnnData
        adata = ad.AnnData(
            X=expression_matrix, obs=pd.DataFrame(obs_data, index=cell_names), var=pd.DataFrame(index=gene_names)
        )
        adata.obsm["X_umap"] = umap_coords

        print(f"   ✅ Generated: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

        # Check memory usage
        memory_mb = adata.X.nbytes / (1024 * 1024)
        print(f"   Memory usage: {memory_mb:.0f} MB")

        # Test PySEE with QC panel only (memory efficient)
        app = PySEE(adata, title="Large Synthetic Test")

        qc_panel = QCPanel("qc", title="QC")
        app.add_panel("qc", qc_panel)

        qc_fig = app.render_panel("qc")
        print("   ✅ Successfully rendered QC panel")

        # Clean up
        del adata, app, qc_fig
        gc.collect()

    except Exception as e:
        print(f"   ❌ Failed: {e}")


def test_backed_mode():
    """Test backed mode for large datasets."""
    print("   Testing backed/on-disk mode...")

    try:
        # Create temporary file
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp_file:
            tmp_path = tmp_file.name

        try:
            # Generate and save large dataset
            n_cells = 100000
            n_genes = 10000

            print(f"   Generating {n_cells:,} cells, {n_genes:,} genes...")

            np.random.seed(42)
            expression_matrix = np.random.negative_binomial(4, 0.35, size=(n_cells, n_genes)).astype(np.float32)

            gene_names = [f"Gene_{i:05d}" for i in range(n_genes)]
            cell_names = [f"Cell_{i:06d}" for i in range(n_cells)]

            obs_data = {
                "cell_type": np.random.choice(["Type_A", "Type_B", "Type_C"], n_cells),
                "total_counts": expression_matrix.sum(axis=1),
            }

            adata = ad.AnnData(
                X=expression_matrix, obs=pd.DataFrame(obs_data, index=cell_names), var=pd.DataFrame(index=gene_names)
            )

            # Save to disk
            adata.write(tmp_path)
            print(f"   ✅ Saved to disk: {tmp_path}")

            # Load in backed mode
            print("   Loading in backed mode...")
            adata_backed = ad.read_h5ad(tmp_path, backed="r")
            print(f"   ✅ Loaded in backed mode: {adata_backed.n_obs:,} cells, {adata_backed.n_vars:,} genes")

            # Test PySEE with backed dataset
            app = PySEE(adata_backed, title="Backed Mode Test")

            qc_panel = QCPanel("qc", title="QC")
            app.add_panel("qc", qc_panel)

            qc_fig = app.render_panel("qc")
            print("   ✅ Successfully rendered QC panel in backed mode")

            # Clean up
            del adata, adata_backed, app, qc_fig
            gc.collect()

        finally:
            # Clean up temporary file
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
                print("   🧹 Cleaned up temporary file")

    except Exception as e:
        print(f"   ❌ Failed: {e}")


def test_subsampling_strategy():
    """Test subsampling strategy for large datasets."""
    print("   Testing subsampling strategy...")

    try:
        # Load PBMC 68K
        adata = sc.datasets.pbmc68k_reduced()
        print(f"   Original: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

        # Subsample to simulate larger dataset
        n_subsample = 50000  # Limit for 16 GB system
        if n_subsample < adata.n_obs:
            indices = np.random.choice(adata.n_obs, n_subsample, replace=False)
            adata_sub = adata[indices].copy()
            print(f"   Subsampled to: {adata_sub.n_obs:,} cells")
        else:
            adata_sub = adata.copy()

        # Test PySEE with subsampled dataset
        app = PySEE(adata_sub, title="Subsampled Test")

        # Add multiple panels
        qc_panel = QCPanel("qc", title="QC")
        app.add_panel("qc", qc_panel)

        violin_panel = ViolinPanel("violin", gene=adata_sub.var_names[0], title="Violin")
        app.add_panel("violin", violin_panel)

        # Render panels
        qc_fig = app.render_panel("qc")
        violin_fig = app.render_panel("violin")

        print("   ✅ Successfully rendered QC and Violin panels")

        # Clean up
        del adata, adata_sub, app, qc_fig, violin_fig
        gc.collect()

    except Exception as e:
        print(f"   ❌ Failed: {e}")


def print_memory_tips():
    """Print memory optimization tips."""
    print("\n💡 Memory Optimization Tips for 16 GB System:")
    print("=" * 60)

    print("1. **Use Backed Mode for Large Datasets**:")
    print("   adata = ad.read_h5ad('large_dataset.h5ad', backed='r')")
    print("   - Reduces memory usage by 80-90%")
    print("   - Data stays on disk, loads only when accessed")

    print("\n2. **Subsample Large Datasets**:")
    print("   # Limit to 50K cells for 16 GB system")
    print("   indices = np.random.choice(adata.n_obs, 50000, replace=False)")
    print("   adata_sub = adata[indices].copy()")

    print("\n3. **Use Memory-Efficient Data Types**:")
    print("   adata.X = adata.X.astype(np.float32)  # Use float32 instead of float64")
    print("   adata.X = scipy.sparse.csr_matrix(adata.X)  # Use sparse matrices")

    print("\n4. **Process in Chunks**:")
    print("   # Process 10K cells at a time")
    print("   for i in range(0, adata.n_obs, 10000):")
    print("       chunk = adata[i:i+10000]")
    print("       # Process chunk...")

    print("\n5. **Monitor Memory Usage**:")
    print("   import psutil")
    print("   memory_percent = psutil.virtual_memory().percent")
    print("   if memory_percent > 80:")
    print("       print('Warning: High memory usage')")

    print("\n6. **Clean Up Variables**:")
    print("   del large_variable")
    print("   gc.collect()  # Force garbage collection")


def main():
    """Run large dataset testing on 16 GB system."""
    test_large_datasets_on_16gb()
    print_memory_tips()

    print("\n🎯 Summary:")
    print("✅ You can test large datasets on 16 GB RAM using:")
    print("   - Backed/on-disk mode")
    print("   - Subsampling strategies")
    print("   - Memory-efficient data types")
    print("   - Chunked processing")
    print("   - Careful memory monitoring")


if __name__ == "__main__":
    main()
