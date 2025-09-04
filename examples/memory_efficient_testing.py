"""
Memory-efficient testing strategies for large datasets on limited RAM.

This script demonstrates how to test large and very large datasets
even with limited system memory (e.g., 16 GB RAM).
"""

import sys
from pathlib import Path
import numpy as np
import anndata as ad
import scanpy as sc
from typing import Dict, Any, Optional

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from pysee import PySEE
from pysee.panels.umap import UMAPPanel
from pysee.panels.violin import ViolinPanel
from pysee.panels.heatmap import HeatmapPanel
from pysee.panels.qc import QCPanel
from pysee.utils.system_requirements import SystemRequirementsChecker


class MemoryEfficientTester:
    """Test large datasets efficiently on limited memory systems."""
    
    def __init__(self):
        self.checker = SystemRequirementsChecker()
        self.system_info = self.checker.get_system_info()
        
    def test_large_dataset_strategies(self):
        """Demonstrate different strategies for testing large datasets."""
        print("🧠 Memory-Efficient Testing Strategies")
        print("=" * 60)
        print(f"System Memory: {self.system_info['total_memory_gb']:.1f} GB")
        print(f"Available Memory: {self.system_info['available_memory_gb']:.1f} GB")
        
        # Strategy 1: Generate large synthetic dataset
        print("\n📊 Strategy 1: Large Synthetic Dataset")
        self._test_large_synthetic()
        
        # Strategy 2: Simulate large dataset with subsampling
        print("\n📊 Strategy 2: Subsampling Approach")
        self._test_subsampling_approach()
        
        # Strategy 3: Backed/on-disk mode simulation
        print("\n📊 Strategy 3: Backed Mode Simulation")
        self._test_backed_mode()
        
        # Strategy 4: Chunked processing
        print("\n📊 Strategy 4: Chunked Processing")
        self._test_chunked_processing()
    
    def _test_large_synthetic(self):
        """Test with a large synthetic dataset."""
        print("   Generating large synthetic dataset (100K cells)...")
        
        # Generate large synthetic dataset
        n_cells = 100000
        n_genes = 15000
        
        # Check memory requirements
        estimated_memory_mb = (n_cells * n_genes * 4) / (1024 * 1024)  # 4 bytes per float32
        print(f"   Estimated memory: {estimated_memory_mb:.0f} MB")
        
        if estimated_memory_mb > self.system_info['available_memory_gb'] * 1024 * 0.8:
            print("   ⚠️ Dataset too large for available memory, using smaller version")
            n_cells = 50000
            n_genes = 10000
            estimated_memory_mb = (n_cells * n_genes * 4) / (1024 * 1024)
            print(f"   Adjusted to: {n_cells:,} cells, {n_genes:,} genes ({estimated_memory_mb:.0f} MB)")
        
        # Generate dataset
        np.random.seed(42)
        expression_matrix = np.random.negative_binomial(4, 0.35, size=(n_cells, n_genes)).astype(np.float32)
        
        # Create gene and cell names
        gene_names = [f"Gene_{i:05d}" for i in range(n_genes)]
        cell_names = [f"Cell_{i:06d}" for i in range(n_cells)]
        
        # Create metadata
        obs_data = {
            'cell_type': np.random.choice(['Type_A', 'Type_B', 'Type_C', 'Type_D'], n_cells),
            'batch': np.random.choice(['Batch_1', 'Batch_2'], n_cells),
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
        
        print(f"   ✅ Generated dataset: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
        
        # Test PySEE with this dataset
        self._test_pysee_with_dataset(adata, "Large Synthetic")
    
    def _test_subsampling_approach(self):
        """Test with subsampling approach."""
        print("   Loading PBMC 68K and subsampling to simulate larger dataset...")
        
        # Load PBMC 68K
        adata = sc.datasets.pbmc68k_reduced()
        print(f"   Original: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
        
        # Subsample to create "larger" dataset simulation
        n_subsample = min(50000, adata.n_obs)  # Limit based on available memory
        if n_subsample < adata.n_obs:
            # Subsample cells
            indices = np.random.choice(adata.n_obs, n_subsample, replace=False)
            adata_sub = adata[indices].copy()
            print(f"   Subsampled to: {adata_sub.n_obs:,} cells")
        else:
            adata_sub = adata.copy()
        
        # Test PySEE with subsampled dataset
        self._test_pysee_with_dataset(adata_sub, "Subsampled PBMC 68K")
    
    def _test_backed_mode(self):
        """Test backed/on-disk mode simulation."""
        print("   Simulating backed/on-disk mode with memory mapping...")
        
        # Create a temporary file for backed mode
        import tempfile
        import os
        
        # Generate dataset and save to disk
        n_cells = 100000
        n_genes = 10000
        
        # Create temporary file
        with tempfile.NamedTemporaryFile(suffix='.h5ad', delete=False) as tmp_file:
            tmp_path = tmp_file.name
        
        try:
            # Generate and save dataset
            print(f"   Generating {n_cells:,} cells, {n_genes:,} genes...")
            np.random.seed(42)
            expression_matrix = np.random.negative_binomial(4, 0.35, size=(n_cells, n_genes)).astype(np.float32)
            
            gene_names = [f"Gene_{i:05d}" for i in range(n_genes)]
            cell_names = [f"Cell_{i:06d}" for i in range(n_cells)]
            
            obs_data = {
                'cell_type': np.random.choice(['Type_A', 'Type_B', 'Type_C'], n_cells),
                'total_counts': expression_matrix.sum(axis=1),
            }
            
            adata = ad.AnnData(
                X=expression_matrix,
                obs=pd.DataFrame(obs_data, index=cell_names),
                var=pd.DataFrame(index=gene_names)
            )
            
            # Save to disk
            adata.write(tmp_path)
            print(f"   ✅ Saved to disk: {tmp_path}")
            
            # Load in backed mode
            print("   Loading in backed mode...")
            adata_backed = ad.read_h5ad(tmp_path, backed='r')
            print(f"   ✅ Loaded in backed mode: {adata_backed.n_obs:,} cells, {adata_backed.n_vars:,} genes")
            
            # Test PySEE with backed dataset
            self._test_pysee_with_dataset(adata_backed, "Backed Mode Dataset")
            
        finally:
            # Clean up
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
                print("   🧹 Cleaned up temporary file")
    
    def _test_chunked_processing(self):
        """Test chunked processing approach."""
        print("   Testing chunked processing for large datasets...")
        
        # Load PBMC 68K
        adata = sc.datasets.pbmc68k_reduced()
        
        # Process in chunks
        chunk_size = 10000
        n_chunks = (adata.n_obs + chunk_size - 1) // chunk_size
        
        print(f"   Processing {adata.n_obs:,} cells in {n_chunks} chunks of {chunk_size:,} cells")
        
        results = []
        for i in range(n_chunks):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, adata.n_obs)
            
            chunk = adata[start_idx:end_idx].copy()
            print(f"   Processing chunk {i+1}/{n_chunks}: {chunk.n_obs:,} cells")
            
            # Test PySEE with chunk
            result = self._test_pysee_with_dataset(chunk, f"Chunk {i+1}")
            results.append(result)
        
        print(f"   ✅ Processed all {n_chunks} chunks successfully")
        return results
    
    def _test_pysee_with_dataset(self, adata: ad.AnnData, dataset_name: str) -> Dict[str, Any]:
        """Test PySEE with a given dataset."""
        print(f"   Testing PySEE with {dataset_name}...")
        
        try:
            # Create PySEE app
            app = PySEE(adata, title=f"Test - {dataset_name}")
            
            # Add panels (only if data supports them)
            panels_added = 0
            
            # QC Panel (always works)
            qc_panel = QCPanel("qc", title="QC")
            app.add_panel("qc", qc_panel)
            panels_added += 1
            
            # Violin Panel (if we have genes)
            if adata.n_vars > 0:
                violin_panel = ViolinPanel("violin", gene=adata.var_names[0], title="Violin")
                app.add_panel("violin", violin_panel)
                panels_added += 1
            
            # UMAP Panel (if we have UMAP coordinates)
            if 'X_umap' in adata.obsm:
                umap_panel = UMAPPanel("umap", title="UMAP")
                app.add_panel("umap", umap_panel)
                panels_added += 1
            
            # Heatmap Panel (for smaller datasets)
            if adata.n_obs <= 50000:  # Limit heatmap to reasonable size
                heatmap_panel = HeatmapPanel("heatmap", title="Heatmap")
                app.add_panel("heatmap", heatmap_panel)
                panels_added += 1
            
            # Render panels
            rendered_panels = 0
            for panel_id in app.panels.keys():
                try:
                    fig = app.render_panel(panel_id)
                    rendered_panels += 1
                except Exception as e:
                    print(f"     ⚠️ Failed to render {panel_id}: {e}")
            
            print(f"   ✅ Success: {panels_added} panels added, {rendered_panels} rendered")
            
            return {
                'dataset_name': dataset_name,
                'n_cells': adata.n_obs,
                'n_genes': adata.n_vars,
                'panels_added': panels_added,
                'panels_rendered': rendered_panels,
                'success': True
            }
            
        except Exception as e:
            print(f"   ❌ Failed: {e}")
            return {
                'dataset_name': dataset_name,
                'n_cells': adata.n_obs,
                'n_genes': adata.n_vars,
                'success': False,
                'error': str(e)
            }
    
    def print_memory_tips(self):
        """Print memory optimization tips."""
        print("\n💡 Memory Optimization Tips for Large Datasets:")
        print("=" * 60)
        
        print("1. **Backed/On-Disk Mode**:")
        print("   adata = ad.read_h5ad('large_dataset.h5ad', backed='r')")
        print("   - Data stays on disk, only loads when accessed")
        print("   - Reduces memory usage by 80-90%")
        
        print("\n2. **Subsampling**:")
        print("   # Subsample to 50K cells for testing")
        print("   indices = np.random.choice(adata.n_obs, 50000, replace=False)")
        print("   adata_sub = adata[indices].copy()")
        
        print("\n3. **Chunked Processing**:")
        print("   # Process in chunks of 10K cells")
        print("   for i in range(0, adata.n_obs, 10000):")
        print("       chunk = adata[i:i+10000]")
        print("       # Process chunk...")
        
        print("\n4. **Memory Monitoring**:")
        print("   # Monitor memory usage")
        print("   import psutil")
        print("   memory_percent = psutil.virtual_memory().percent")
        print("   if memory_percent > 80:")
        print("       print('Warning: High memory usage')")
        
        print("\n5. **Data Type Optimization**:")
        print("   # Use float32 instead of float64")
        print("   adata.X = adata.X.astype(np.float32)")
        print("   # Use sparse matrices")
        print("   adata.X = scipy.sparse.csr_matrix(adata.X)")


def main():
    """Run memory-efficient testing demonstration."""
    print("🧠 Memory-Efficient Testing for Large Datasets")
    print("=" * 60)
    
    tester = MemoryEfficientTester()
    
    # Run all strategies
    tester.test_large_dataset_strategies()
    
    # Print tips
    tester.print_memory_tips()
    
    print("\n✅ Memory-efficient testing demonstration complete!")
    print("💡 Use these strategies to test large datasets on your 16 GB system")


if __name__ == "__main__":
    main()
