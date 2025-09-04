"""
Dataset downloader for PySEE performance testing.

This script downloads real bioinformatics datasets for comprehensive
performance testing of PySEE components.
"""

import os
import sys
import requests
import gzip
import tempfile
from pathlib import Path
from typing import Optional, Dict, Any
import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc


class DatasetDownloader:
    """Download real datasets for performance testing."""
    
    def __init__(self, data_dir: str = "data"):
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(exist_ok=True)
    
    def download_pbmc3k(self) -> ad.AnnData:
        """Download PBMC 3K dataset."""
        print("📥 Downloading PBMC 3K dataset...")
        adata = sc.datasets.pbmc3k()
        
        # Save to file
        output_file = self.data_dir / "pbmc3k.h5ad"
        adata.write(output_file)
        print(f"✅ Saved to: {output_file}")
        
        return adata
    
    def download_pbmc68k(self) -> ad.AnnData:
        """Download PBMC 68K dataset."""
        print("📥 Downloading PBMC 68K dataset...")
        adata = sc.datasets.pbmc68k_reduced()
        
        # Save to file
        output_file = self.data_dir / "pbmc68k.h5ad"
        adata.write(output_file)
        print(f"✅ Saved to: {output_file}")
        
        return adata
    
    def download_tabula_sapiens_sample(self) -> Optional[ad.AnnData]:
        """Download a sample from Tabula Sapiens dataset."""
        print("📥 Downloading Tabula Sapiens sample...")
        
        # Tabula Sapiens download URL (example - actual URL may differ)
        url = "https://tabula-sapiens-portal.ds.czbiohub.org/static/downloads/tabula-sapiens.h5ad"
        
        try:
            output_file = self.data_dir / "tabula_sapiens_sample.h5ad"
            
            if output_file.exists():
                print(f"✅ File already exists: {output_file}")
                return ad.read_h5ad(output_file)
            
            print("⚠️ Tabula Sapiens download not implemented yet")
            print("   This would download a large dataset (~2GB)")
            print("   For now, using synthetic data instead")
            
            return None
            
        except Exception as e:
            print(f"❌ Failed to download Tabula Sapiens: {e}")
            return None
    
    def download_human_cell_atlas_sample(self) -> Optional[ad.AnnData]:
        """Download a sample from Human Cell Atlas dataset."""
        print("📥 Downloading Human Cell Atlas sample...")
        
        try:
            print("⚠️ Human Cell Atlas download not implemented yet")
            print("   This would download a very large dataset (~10GB+)")
            print("   For now, using synthetic data instead")
            
            return None
            
        except Exception as e:
            print(f"❌ Failed to download Human Cell Atlas: {e}")
            return None
    
    def generate_large_synthetic_dataset(self, n_cells: int = 100000, n_genes: int = 15000) -> ad.AnnData:
        """Generate a large synthetic dataset for stress testing."""
        print(f"🔬 Generating large synthetic dataset: {n_cells:,} cells, {n_genes:,} genes...")
        
        np.random.seed(42)
        
        # Create expression matrix with realistic structure
        expression_matrix = np.random.negative_binomial(4, 0.35, size=(n_cells, n_genes)).astype(np.float32)
        
        # Add some highly variable genes
        n_hvg = int(0.1 * n_genes)
        hvg_indices = np.random.choice(n_genes, n_hvg, replace=False)
        expression_matrix[:, hvg_indices] *= np.random.lognormal(0, 1, n_hvg)
        
        # Add mitochondrial genes
        n_mito_genes = int(0.1 * n_genes)
        mito_indices = np.random.choice(n_genes, n_mito_genes, replace=False)
        expression_matrix[:, mito_indices] *= 2.0
        
        # Create gene names
        gene_names = [f"Gene_{i:05d}" for i in range(n_genes)]
        for i, idx in enumerate(mito_indices):
            gene_names[idx] = f"MT-Gene_{i:04d}"
        
        # Create cell names
        cell_names = [f"Cell_{i:06d}" for i in range(n_cells)]
        
        # Create cell metadata
        obs_data = {
            'cell_type': np.random.choice([
                'T_cell', 'B_cell', 'NK_cell', 'Monocyte', 'Dendritic', 
                'Neutrophil', 'Eosinophil', 'Basophil', 'Mast_cell'
            ], n_cells),
            'batch': np.random.choice(['Batch_1', 'Batch_2', 'Batch_3', 'Batch_4'], n_cells),
            'total_counts': expression_matrix.sum(axis=1),
            'detected_genes': (expression_matrix > 0).sum(axis=1),
        }
        
        # Add mitochondrial percentage
        mito_genes = [name.startswith('MT-') for name in gene_names]
        mito_counts = expression_matrix[:, mito_genes].sum(axis=1)
        total_counts = expression_matrix.sum(axis=1)
        obs_data['mito_percent'] = mito_counts / total_counts * 100
        
        # Add some UMAP coordinates
        np.random.seed(42)
        umap_coords = np.random.randn(n_cells, 2)
        
        # Create AnnData
        adata = ad.AnnData(
            X=expression_matrix,
            obs=pd.DataFrame(obs_data, index=cell_names),
            var=pd.DataFrame(index=gene_names)
        )
        
        # Add UMAP coordinates
        adata.obsm['X_umap'] = umap_coords
        
        # Save to file
        output_file = self.data_dir / f"synthetic_large_{n_cells//1000}k_cells.h5ad"
        adata.write(output_file)
        print(f"✅ Saved to: {output_file}")
        
        return adata
    
    def get_dataset_info(self, adata: ad.AnnData) -> Dict[str, Any]:
        """Get information about a dataset."""
        return {
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars,
            'memory_usage_mb': adata.X.nbytes / (1024 * 1024),
            'sparsity': 1 - (adata.X > 0).sum() / (adata.n_obs * adata.n_vars),
            'has_umap': 'X_umap' in adata.obsm,
            'cell_types': list(adata.obs['cell_type'].unique()) if 'cell_type' in adata.obs else [],
            'batches': list(adata.obs['batch'].unique()) if 'batch' in adata.obs else [],
        }
    
    def download_all_datasets(self) -> Dict[str, ad.AnnData]:
        """Download all available datasets."""
        datasets = {}
        
        print("🚀 Downloading all test datasets...")
        print("=" * 50)
        
        # Download built-in datasets
        try:
            datasets['pbmc3k'] = self.download_pbmc3k()
        except Exception as e:
            print(f"❌ Failed to download PBMC 3K: {e}")
        
        try:
            datasets['pbmc68k'] = self.download_pbmc68k()
        except Exception as e:
            print(f"❌ Failed to download PBMC 68K: {e}")
        
        # Generate synthetic datasets
        try:
            datasets['synthetic_medium'] = self.generate_large_synthetic_dataset(10000, 5000)
        except Exception as e:
            print(f"❌ Failed to generate medium synthetic dataset: {e}")
        
        try:
            datasets['synthetic_large'] = self.generate_large_synthetic_dataset(50000, 10000)
        except Exception as e:
            print(f"❌ Failed to generate large synthetic dataset: {e}")
        
        try:
            datasets['synthetic_very_large'] = self.generate_large_synthetic_dataset(100000, 15000)
        except Exception as e:
            print(f"❌ Failed to generate very large synthetic dataset: {e}")
        
        # Print summary
        print("\n📊 Dataset Summary:")
        print("=" * 50)
        
        for name, adata in datasets.items():
            info = self.get_dataset_info(adata)
            print(f"{name:20} | {info['n_cells']:8,} cells | {info['n_genes']:6,} genes | {info['memory_usage_mb']:6.1f} MB")
        
        return datasets


def main():
    """Download all test datasets."""
    downloader = DatasetDownloader()
    
    try:
        datasets = downloader.download_all_datasets()
        
        print(f"\n✅ Successfully downloaded {len(datasets)} datasets")
        print("📁 Datasets saved in: data/")
        print("\n💡 Next steps:")
        print("   - Run performance tests: python tests/performance/run_performance_tests.py")
        print("   - Use datasets in your PySEE analysis")
        
        return 0
        
    except Exception as e:
        print(f"❌ Dataset download failed: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(main())
