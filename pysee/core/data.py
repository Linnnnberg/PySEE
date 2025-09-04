"""
AnnData integration and data handling for PySEE.

This module provides the AnnDataWrapper class that handles AnnData objects
with proper validation, preprocessing, and metadata extraction.
"""

import warnings
from typing import Optional, Union, List, Dict, Any, Tuple
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse


class AnnDataWrapper:
    """
    Wrapper class for AnnData objects with validation and preprocessing capabilities.
    
    This class provides a standardized interface for working with AnnData objects
    in PySEE, including data validation, preprocessing utilities, and metadata extraction.
    
    Parameters
    ----------
    adata : anndata.AnnData
        The AnnData object to wrap
    validate : bool, default True
        Whether to validate the AnnData object on initialization
    """
    
    def __init__(self, adata: ad.AnnData, validate: bool = True):
        self._adata = adata
        self._original_adata = adata.copy()
        
        if validate:
            self.validate()
    
    @property
    def adata(self) -> ad.AnnData:
        """Return the wrapped AnnData object."""
        return self._adata
    
    @property
    def original_adata(self) -> ad.AnnData:
        """Return the original AnnData object (before any modifications)."""
        return self._original_adata
    
    def validate(self) -> None:
        """
        Validate the AnnData object for PySEE compatibility.
        
        Raises
        ------
        ValueError
            If the AnnData object is not compatible with PySEE
        """
        if not isinstance(self._adata, ad.AnnData):
            raise ValueError("Input must be an AnnData object")
        
        if self._adata.n_obs == 0:
            raise ValueError("AnnData object has no observations (cells)")
        
        if self._adata.n_vars == 0:
            raise ValueError("AnnData object has no variables (genes)")
        
        # Check for required attributes
        if self._adata.X is None:
            raise ValueError("AnnData object must have expression data (X)")
        
        # Warn about potential issues
        if sparse.issparse(self._adata.X) and self._adata.X.nnz == 0:
            warnings.warn("Expression matrix appears to be empty", UserWarning)
        
        if self._adata.n_obs < 10:
            warnings.warn(f"Very few observations ({self._adata.n_obs}), some visualizations may not work well", UserWarning)
        
        if self._adata.n_vars < 10:
            warnings.warn(f"Very few variables ({self._adata.n_vars}), some visualizations may not work well", UserWarning)
    
    def get_obs_columns(self) -> List[str]:
        """Get list of observation (cell) metadata columns."""
        return list(self._adata.obs.columns)
    
    def get_var_columns(self) -> List[str]:
        """Get list of variable (gene) metadata columns."""
        return list(self._adata.var.columns)
    
    def get_obsm_keys(self) -> List[str]:
        """Get list of observation embedding keys (e.g., 'X_umap', 'X_pca')."""
        return list(self._adata.obsm.keys())
    
    def get_varm_keys(self) -> List[str]:
        """Get list of variable embedding keys."""
        return list(self._adata.varm.keys())
    
    def get_uns_keys(self) -> List[str]:
        """Get list of unstructured annotation keys."""
        return list(self._adata.uns.keys())
    
    def get_categorical_obs_columns(self) -> List[str]:
        """Get observation columns that are categorical."""
        categorical_cols = []
        for col in self._adata.obs.columns:
            if self._adata.obs[col].dtype == 'category' or self._adata.obs[col].dtype == 'object':
                categorical_cols.append(col)
        return categorical_cols
    
    def get_numerical_obs_columns(self) -> List[str]:
        """Get observation columns that are numerical."""
        numerical_cols = []
        for col in self._adata.obs.columns:
            if pd.api.types.is_numeric_dtype(self._adata.obs[col]):
                numerical_cols.append(col)
        return numerical_cols
    
    def get_embedding_data(self, key: str) -> np.ndarray:
        """
        Get embedding data from obsm.
        
        Parameters
        ----------
        key : str
            Key for the embedding (e.g., 'X_umap', 'X_pca')
            
        Returns
        -------
        np.ndarray
            Embedding coordinates
            
        Raises
        ------
        KeyError
            If the embedding key doesn't exist
        """
        if key not in self._adata.obsm:
            available_keys = list(self._adata.obsm.keys())
            raise KeyError(f"Embedding '{key}' not found. Available embeddings: {available_keys}")
        
        return self._adata.obsm[key]
    
    def get_expression_data(self, genes: Optional[Union[str, List[str]]] = None) -> np.ndarray:
        """
        Get expression data for specified genes.
        
        Parameters
        ----------
        genes : str or list of str, optional
            Gene names to extract. If None, returns all genes.
            
        Returns
        -------
        np.ndarray
            Expression data matrix
        """
        if genes is None:
            return self._adata.X.toarray() if sparse.issparse(self._adata.X) else self._adata.X
        
        if isinstance(genes, str):
            genes = [genes]
        
        # Find gene indices
        gene_indices = []
        for gene in genes:
            if gene in self._adata.var_names:
                gene_indices.append(self._adata.var_names.get_loc(gene))
            else:
                warnings.warn(f"Gene '{gene}' not found in dataset", UserWarning)
        
        if not gene_indices:
            raise ValueError("None of the specified genes found in dataset")
        
        # Extract expression data
        if sparse.issparse(self._adata.X):
            return self._adata.X[:, gene_indices].toarray()
        else:
            return self._adata.X[:, gene_indices]
    
    def get_cell_subset(self, obs_mask: Union[str, np.ndarray, List[bool]]) -> 'AnnDataWrapper':
        """
        Get a subset of cells based on observation mask.
        
        Parameters
        ----------
        obs_mask : str, np.ndarray, or list of bool
            Mask to select cells. If string, should be a column name in obs.
            
        Returns
        -------
        AnnDataWrapper
            New wrapper with subset of cells
        """
        if isinstance(obs_mask, str):
            if obs_mask not in self._adata.obs.columns:
                raise KeyError(f"Column '{obs_mask}' not found in obs")
            mask = self._adata.obs[obs_mask].values
        else:
            mask = np.array(obs_mask)
        
        if len(mask) != self._adata.n_obs:
            raise ValueError(f"Mask length ({len(mask)}) doesn't match number of observations ({self._adata.n_obs})")
        
        subset_adata = self._adata[mask].copy()
        return AnnDataWrapper(subset_adata, validate=False)
    
    def get_gene_subset(self, var_mask: Union[str, np.ndarray, List[bool]]) -> 'AnnDataWrapper':
        """
        Get a subset of genes based on variable mask.
        
        Parameters
        ----------
        var_mask : str, np.ndarray, or list of bool
            Mask to select genes. If string, should be a column name in var.
            
        Returns
        -------
        AnnDataWrapper
            New wrapper with subset of genes
        """
        if isinstance(var_mask, str):
            if var_mask not in self._adata.var.columns:
                raise KeyError(f"Column '{var_mask}' not found in var")
            mask = self._adata.var[var_mask].values
        else:
            mask = np.array(var_mask)
        
        if len(mask) != self._adata.n_vars:
            raise ValueError(f"Mask length ({len(mask)}) doesn't match number of variables ({self._adata.n_vars})")
        
        subset_adata = self._adata[:, mask].copy()
        return AnnDataWrapper(subset_adata, validate=False)
    
    def get_summary_stats(self) -> Dict[str, Any]:
        """
        Get summary statistics about the dataset.
        
        Returns
        -------
        dict
            Dictionary containing summary statistics
        """
        stats = {
            'n_obs': self._adata.n_obs,
            'n_vars': self._adata.n_vars,
            'obs_columns': self.get_obs_columns(),
            'var_columns': self.get_var_columns(),
            'obsm_keys': self.get_obsm_keys(),
            'varm_keys': self.get_varm_keys(),
            'uns_keys': self.get_uns_keys(),
            'categorical_obs_columns': self.get_categorical_obs_columns(),
            'numerical_obs_columns': self.get_numerical_obs_columns(),
        }
        
        # Add expression matrix info
        if sparse.issparse(self._adata.X):
            stats['sparse_matrix'] = True
            stats['nnz'] = self._adata.X.nnz
            stats['sparsity'] = 1 - (self._adata.X.nnz / (self._adata.n_obs * self._adata.n_vars))
        else:
            stats['sparse_matrix'] = False
        
        return stats
    
    def reset_to_original(self) -> None:
        """Reset the AnnData object to its original state."""
        self._adata = self._original_adata.copy()
    
    def copy(self) -> 'AnnDataWrapper':
        """Create a copy of the AnnDataWrapper."""
        return AnnDataWrapper(self._adata.copy(), validate=False)
    
    def __repr__(self) -> str:
        """String representation of the AnnDataWrapper."""
        return f"AnnDataWrapper(n_obs={self._adata.n_obs}, n_vars={self._adata.n_vars})"
    
    def __len__(self) -> int:
        """Return number of observations."""
        return self._adata.n_obs
