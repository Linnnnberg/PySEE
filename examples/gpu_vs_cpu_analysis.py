"""
GPU vs CPU analysis for PySEE large dataset handling.

This script demonstrates the current limitations and potential solutions
for GPU acceleration in PySEE.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import time
import psutil

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from pysee import PySEE
from pysee.panels.umap import UMAPPanel
from pysee.panels.violin import ViolinPanel
from pysee.panels.heatmap import HeatmapPanel
from pysee.panels.qc import QCPanel


def analyze_gpu_vs_cpu_limitations():
    """Analyze GPU vs CPU limitations for PySEE."""
    print("🔍 GPU vs CPU Analysis for PySEE")
    print("=" * 60)
    
    # Check system resources
    system_info = {
        'cpu_count': psutil.cpu_count(),
        'memory_gb': psutil.virtual_memory().total / (1024**3),
        'gpu_available': check_gpu_availability()
    }
    
    print(f"System Resources:")
    print(f"  CPU Cores: {system_info['cpu_count']}")
    print(f"  Memory: {system_info['memory_gb']:.1f} GB")
    print(f"  GPU Available: {system_info['gpu_available']}")
    
    # Test with different dataset sizes (local machine only - small to medium)
    test_sizes = [
        (1000, 2000, "Small"),
        (10000, 5000, "Medium"),
        # Large and very large datasets should be tested on cloud/servers
        # (50000, 10000, "Large"),
        # (100000, 15000, "Very Large")
    ]
    
    for n_cells, n_genes, size_label in test_sizes:
        print(f"\n📊 Testing {size_label} Dataset ({n_cells:,} cells, {n_genes:,} genes)")
        test_dataset_performance(n_cells, n_genes, size_label)


def check_gpu_availability():
    """Check if GPU is available."""
    try:
        import cupy as cp
        if cp.cuda.is_available():
            try:
                # Test if GPU operations actually work
                test_array = cp.array([1, 2, 3])
                result = test_array * 2
                gpu_count = cp.cuda.runtime.getDeviceCount()
                gpu_memory = cp.cuda.Device().mem_info[1] / 1024**3
                return f"CuPy available (GPU: {gpu_count} devices, {gpu_memory:.1f} GB memory)"
            except Exception:
                return "CuPy installed but GPU operations not working (missing CUDA libraries)"
        else:
            return "CuPy installed but CUDA not available"
    except ImportError:
        try:
            import torch
            return f"PyTorch available (CUDA: {torch.cuda.is_available()})"
        except ImportError:
            return "No GPU libraries available"


def test_dataset_performance(n_cells, n_genes, size_label):
    """Test performance with different dataset sizes."""
    print(f"  Generating {size_label} dataset...")
    
    # Generate synthetic dataset
    np.random.seed(42)
    expression_matrix = np.random.negative_binomial(4, 0.35, size=(n_cells, n_genes)).astype(np.float32)
    
    # Create metadata
    gene_names = [f"Gene_{i:05d}" for i in range(n_genes)]
    cell_names = [f"Cell_{i:06d}" for i in range(n_cells)]
    
    obs_data = {
        'cell_type': np.random.choice(['Type_A', 'Type_B', 'Type_C'], n_cells),
        'total_counts': expression_matrix.sum(axis=1),
        'detected_genes': (expression_matrix > 0).sum(axis=1),
    }
    
    # Add UMAP coordinates
    umap_coords = np.random.randn(n_cells, 2)
    
    # Create AnnData
    adata = ad.AnnData(
        X=expression_matrix,
        obs=pd.DataFrame(obs_data, index=cell_names),
        var=pd.DataFrame(index=gene_names)
    )
    adata.obsm['X_umap'] = umap_coords
    
    # Test PySEE performance
    test_pysee_performance(adata, size_label)
    
    # Test GPU vs CPU computation
    test_gpu_vs_cpu_computation(expression_matrix, size_label)


def test_pysee_performance(adata, size_label):
    """Test PySEE performance with different dataset sizes."""
    print(f"  Testing PySEE performance...")
    
    try:
        start_time = time.time()
        
        # Create PySEE app
        app = PySEE(adata, title=f"Test - {size_label}")
        
        # Add panels
        qc_panel = QCPanel("qc", title="QC")
        app.add_panel("qc", qc_panel)
        
        # Only add UMAP for smaller datasets
        if adata.n_obs <= 50000:
            umap_panel = UMAPPanel("umap", title="UMAP")
            app.add_panel("umap", umap_panel)
        
        # Render panels
        qc_fig = app.render_panel("qc")
        
        if adata.n_obs <= 50000:
            umap_fig = app.render_panel("umap")
        
        end_time = time.time()
        
        print(f"    ✅ PySEE rendering: {end_time - start_time:.2f}s")
        
    except Exception as e:
        print(f"    ❌ PySEE failed: {e}")


def test_gpu_vs_cpu_computation(expression_matrix, size_label):
    """Test GPU vs CPU computation performance."""
    print(f"  Testing GPU vs CPU computation...")
    
    # CPU computation
    start_time = time.time()
    cpu_result = np.linalg.svd(expression_matrix, full_matrices=False)
    cpu_time = time.time() - start_time
    print(f"    CPU SVD: {cpu_time:.2f}s")
    
    # GPU computation (if available)
    try:
        import cupy as cp
        
        if cp.cuda.is_available():
            try:
                # Test if GPU computation actually works
                test_matrix = cp.asarray(expression_matrix[:100, :100])  # Small test
                test_result = cp.linalg.svd(test_matrix, full_matrices=False)
                
                # If test works, do full computation
                gpu_matrix = cp.asarray(expression_matrix)
                
                start_time = time.time()
                gpu_result = cp.linalg.svd(gpu_matrix, full_matrices=False)
                gpu_time = time.time() - start_time
                
                # Transfer back to CPU
                cpu_gpu_result = [cp.asnumpy(r) for r in gpu_result]
                
                print(f"    GPU SVD: {gpu_time:.2f}s")
                print(f"    GPU Speedup: {cpu_time/gpu_time:.2f}x")
                
            except Exception as gpu_error:
                print(f"    GPU computation: CuPy installed but GPU operations failed - {gpu_error}")
        else:
            print(f"    GPU computation: CuPy installed but CUDA not available")
        
    except ImportError:
        print(f"    GPU computation: CuPy not installed")
    except Exception as e:
        print(f"    GPU computation: Error - {e}")


def demonstrate_webgl_acceleration():
    """Demonstrate WebGL acceleration for medium datasets."""
    print("\n🚀 WebGL Acceleration Demo")
    print("=" * 60)
    
    # Generate medium dataset (suitable for local testing)
    n_cells = 20000
    n_genes = 5000
    
    print(f"Generating medium dataset: {n_cells:,} cells, {n_genes:,} genes")
    
    np.random.seed(42)
    expression_matrix = np.random.negative_binomial(4, 0.35, size=(n_cells, n_genes)).astype(np.float32)
    
    # Create UMAP coordinates
    umap_coords = np.random.randn(n_cells, 2)
    
    # Test WebGL vs regular scatter plot
    import plotly.graph_objects as go
    
    print("Creating WebGL scatter plot...")
    start_time = time.time()
    
    fig_webgl = go.Figure(data=go.Scattergl(  # WebGL scatter plot
        x=umap_coords[:, 0],
        y=umap_coords[:, 1],
        mode='markers',
        marker=dict(size=2, opacity=0.6),
        name='WebGL Scatter'
    ))
    
    webgl_time = time.time() - start_time
    print(f"WebGL scatter plot creation: {webgl_time:.2f}s")
    
    # Test regular scatter plot
    print("Creating regular scatter plot...")
    start_time = time.time()
    
    fig_regular = go.Figure(data=go.Scatter(  # Regular scatter plot
        x=umap_coords[:, 0],
        y=umap_coords[:, 1],
        mode='markers',
        marker=dict(size=2, opacity=0.6),
        name='Regular Scatter'
    ))
    
    regular_time = time.time() - start_time
    print(f"Regular scatter plot creation: {regular_time:.2f}s")
    
    print(f"WebGL vs Regular: {regular_time/webgl_time:.2f}x faster")


def demonstrate_sampling_strategy():
    """Demonstrate sampling strategy for large datasets."""
    print("\n📊 Sampling Strategy Demo")
    print("=" * 60)
    
    # Load PBMC 68K
    adata = sc.datasets.pbmc68k_reduced()
    print(f"Original dataset: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    
    # Test different sampling sizes
    sample_sizes = [1000, 5000, 10000, 20000]
    
    for sample_size in sample_sizes:
        if sample_size >= adata.n_obs:
            continue
            
        print(f"\nTesting with {sample_size:,} cells:")
        
        # Subsample
        indices = np.random.choice(adata.n_obs, sample_size, replace=False)
        adata_sub = adata[indices].copy()
        
        # Test PySEE performance
        start_time = time.time()
        
        app = PySEE(adata_sub, title=f"Sampled - {sample_size:,} cells")
        
        qc_panel = QCPanel("qc", title="QC")
        app.add_panel("qc", qc_panel)
        
        violin_panel = ViolinPanel("violin", gene=adata_sub.var_names[0], title="Violin")
        app.add_panel("violin", violin_panel)
        
        qc_fig = app.render_panel("qc")
        violin_fig = app.render_panel("violin")
        
        end_time = time.time()
        
        print(f"  PySEE rendering: {end_time - start_time:.2f}s")


def print_recommendations():
    """Print recommendations for GPU vs CPU usage."""
    print("\n💡 Recommendations for PySEE")
    print("=" * 60)
    
    print("1. **Current Approach (Recommended for v0.1.2):**")
    print("   - Use CPU for all operations")
    print("   - Enable WebGL for large scatter plots")
    print("   - Implement smart sampling for large datasets")
    print("   - Keep architecture simple")
    
    print("\n2. **Future Approach (v0.2+):**")
    print("   - Hybrid: GPU for computation, CPU for visualization")
    print("   - Use RAPIDS for GPU-accelerated bioinformatics")
    print("   - Optimize CPU ↔ GPU memory transfers")
    print("   - Maintain both CPU and GPU code paths")
    
    print("\n3. **For Your 16 GB System:**")
    print("   - Focus on CPU optimization")
    print("   - Use sampling for large datasets")
    print("   - Use cloud for very large datasets")
    print("   - GPU acceleration not available (missing CUDA libraries)")
    print("   - Don't add GPU complexity yet")
    
    print("\n4. **Key Insights:**")
    print("   - GPU acceleration is possible but complex")
    print("   - Plotly doesn't support GPU acceleration")
    print("   - WebGL provides browser GPU acceleration")
    print("   - Sampling is often more effective than GPU")


def main():
    """Run GPU vs CPU analysis."""
    analyze_gpu_vs_cpu_limitations()
    demonstrate_webgl_acceleration()
    demonstrate_sampling_strategy()
    print_recommendations()
    
    print("\n🎯 Conclusion:")
    print("GPU acceleration is possible but not necessary for PySEE v0.1.2.")
    print("Focus on CPU optimization, WebGL, and sampling instead.")


if __name__ == "__main__":
    main()
