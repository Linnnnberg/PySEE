"""
Dataset fixtures for performance testing.

This module provides various dataset sizes and types for comprehensive
performance testing of PySEE components.
"""

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from typing import Dict, Any, Optional
import os
import requests
import gzip
import tempfile


class DatasetFixtures:
    """Collection of dataset fixtures for performance testing."""

    @staticmethod
    def get_pbmc3k() -> ad.AnnData:
        """Get PBMC 3K dataset (small benchmark)."""
        return sc.datasets.pbmc3k()

    @staticmethod
    def get_pbmc68k() -> ad.AnnData:
        """Get PBMC 68K dataset (medium benchmark)."""
        return sc.datasets.pbmc68k_reduced()

    @staticmethod
    def generate_synthetic_small() -> ad.AnnData:
        """Generate small synthetic dataset (1K cells, 2K genes)."""
        return DatasetFixtures._generate_synthetic_dataset(n_cells=1000, n_genes=2000, complexity="low")

    @staticmethod
    def generate_synthetic_medium() -> ad.AnnData:
        """Generate medium synthetic dataset (10K cells, 5K genes)."""
        return DatasetFixtures._generate_synthetic_dataset(n_cells=10000, n_genes=5000, complexity="medium")

    @staticmethod
    def generate_synthetic_large() -> ad.AnnData:
        """Generate large synthetic dataset (50K cells, 10K genes)."""
        return DatasetFixtures._generate_synthetic_dataset(n_cells=50000, n_genes=10000, complexity="high")

    @staticmethod
    def generate_synthetic_very_large() -> ad.AnnData:
        """Generate very large synthetic dataset (100K cells, 15K genes)."""
        return DatasetFixtures._generate_synthetic_dataset(n_cells=100000, n_genes=15000, complexity="high")

    @staticmethod
    def _generate_synthetic_dataset(n_cells: int, n_genes: int, complexity: str = "medium") -> ad.AnnData:
        """Generate synthetic single-cell dataset for testing."""
        np.random.seed(42)

        # Base expression matrix
        if complexity == "low":
            expression_matrix = np.random.negative_binomial(3, 0.4, size=(n_cells, n_genes)).astype(float)
        elif complexity == "medium":
            expression_matrix = np.random.negative_binomial(4, 0.35, size=(n_cells, n_genes)).astype(float)
        else:  # high complexity
            expression_matrix = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(float)

        # Add some structure
        n_hvg = int(0.1 * n_genes)  # 10% highly variable genes
        hvg_indices = np.random.choice(n_genes, n_hvg, replace=False)
        expression_matrix[:, hvg_indices] *= np.random.lognormal(0, 1, n_hvg)

        # Add mitochondrial genes
        n_mito_genes = int(0.1 * n_genes)
        mito_indices = np.random.choice(n_genes, n_mito_genes, replace=False)
        expression_matrix[:, mito_indices] *= 2.0

        # Create gene names
        gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
        for i, idx in enumerate(mito_indices):
            gene_names[idx] = f"MT-Gene_{i:03d}"

        # Create cell names
        cell_names = [f"Cell_{i:06d}" for i in range(n_cells)]

        # Create cell metadata
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

        # Add some UMAP coordinates
        np.random.seed(42)
        umap_coords = np.random.randn(n_cells, 2)

        # Create AnnData
        adata = ad.AnnData(
            X=expression_matrix, obs=pd.DataFrame(obs_data, index=cell_names), var=pd.DataFrame(index=gene_names)
        )

        # Add UMAP coordinates
        adata.obsm["X_umap"] = umap_coords

        return adata

    @staticmethod
    def get_all_datasets() -> Dict[str, ad.AnnData]:
        """Get all available datasets for testing."""
        return {
            "pbmc3k": DatasetFixtures.get_pbmc3k(),
            "pbmc68k": DatasetFixtures.get_pbmc68k(),
            "synthetic_small": DatasetFixtures.generate_synthetic_small(),
            "synthetic_medium": DatasetFixtures.generate_synthetic_medium(),
            "synthetic_large": DatasetFixtures.generate_synthetic_large(),
            "synthetic_very_large": DatasetFixtures.generate_synthetic_very_large(),
        }

    @staticmethod
    def get_dataset_info(adata: ad.AnnData) -> Dict[str, Any]:
        """Get information about a dataset."""
        return {
            "n_cells": adata.n_obs,
            "n_genes": adata.n_vars,
            "memory_usage_mb": adata.X.nbytes / (1024 * 1024),
            "sparsity": 1 - (adata.X > 0).sum() / (adata.n_obs * adata.n_vars),
            "has_umap": "X_umap" in adata.obsm,
            "has_obs": len(adata.obs.columns),
            "has_var": len(adata.var.columns),
        }


class DatasetDownloader:
    """Download real datasets for performance testing."""

    @staticmethod
    def download_tabula_sapiens_sample() -> Optional[ad.AnnData]:
        """Download a sample from Tabula Sapiens dataset."""
        # This would implement actual download from Tabula Sapiens
        # For now, return None to indicate not implemented
        return None

    @staticmethod
    def download_human_cell_atlas_sample() -> Optional[ad.AnnData]:
        """Download a sample from Human Cell Atlas dataset."""
        # This would implement actual download from HCA
        # For now, return None to indicate not implemented
        return None

    @staticmethod
    def download_10x_genomics_sample() -> Optional[ad.AnnData]:
        """Download a sample from 10X Genomics public datasets."""
        # This would implement actual download from 10X
        # For now, return None to indicate not implemented
        return None
