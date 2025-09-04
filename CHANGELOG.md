# Changelog

All notable changes to PySEE will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2025-01-04

### Added
- **Core Architecture**
  - `AnnDataWrapper` class for comprehensive AnnData handling and validation
  - `PySEE` dashboard engine with panel management and interaction system
  - `BasePanel` abstract class defining the panel interface

- **Visualization Panels**
  - `UMAPPanel` for interactive dimensionality reduction plots (UMAP, t-SNE, PCA)
  - `ViolinPanel` for gene expression distribution plots with grouping
  - Support for categorical and continuous color mapping
  - Interactive Plotly-based visualizations with hover tooltips

- **Panel Linking & Interaction**
  - Selection propagation system between linked panels
  - Global selection management across all panels
  - Code export functionality for reproducible analysis

- **Data Handling**
  - Comprehensive AnnData validation and preprocessing utilities
  - Support for different AnnData formats (h5ad, zarr)
  - Metadata extraction utilities (obs, var, obsm, varm, uns)
  - Data subsetting and filtering capabilities

- **CLI Interface**
  - Command-line interface for running PySEE dashboards
  - Support for data loading and panel configuration via CLI
  - Code export functionality from command line

- **Documentation & Testing**
  - Comprehensive README with installation and usage instructions
  - Working example scripts with real data (PBMC3K dataset)
  - Test suite with comprehensive functionality testing
  - Master TODO list with detailed project roadmap

- **Package Structure**
  - Complete Python package with proper setup.py
  - Requirements.txt with all necessary dependencies
  - GitHub repository with proper structure and documentation

### Features
- Interactive UMAP/t-SNE/PCA scatter plots with color mapping
- Gene expression violin/box/strip plots with grouping
- Linked panel selections (UMAP selection affects violin plots)
- Reproducible Python code export
- Jupyter notebook integration
- Data validation and error handling
- Configurable panel parameters and styling

### Technical Details
- Built on Plotly for interactive visualizations
- Integrates with Scanpy and AnnData ecosystem
- Supports sparse and dense expression matrices
- Handles categorical and continuous metadata
- Memory-efficient data handling
- Extensible panel architecture

### Dependencies
- Core: scanpy, anndata, pandas, numpy, scipy
- Visualization: plotly, matplotlib, seaborn
- Jupyter: jupyter, ipywidgets, ipykernel
- Bioinformatics: scvi-tools, mudata, zarr
- Clustering: python-igraph, leidenalg
- Development: pytest, black, flake8, mypy

---

## [Unreleased]

### Planned for v0.2
- Heatmap panel for gene expression matrices
- QC metrics panel for data quality assessment
- Dot plot panel for marker gene visualization
- Enhanced selection tools (lasso, polygon selection)
- Jupyter widget integration improvements

### Planned for v0.3
- Genome browser integration (IGV.js)
- Spatial transcriptomics viewer (Vitessce)
- Plugin system for custom panels
- Web deployment capabilities
- Cloud-scale data support

---

## Development Notes

### MVP Achievement
The MVP (v0.1) has been successfully completed with all core features implemented:
- ✅ AnnData integration and validation
- ✅ UMAP and Violin visualization panels
- ✅ Panel linking and selection propagation
- ✅ Code export for reproducibility
- ✅ Jupyter notebook integration
- ✅ CLI interface
- ✅ Comprehensive testing and documentation

### Next Steps
- Focus on additional panel types (heatmaps, QC plots)
- Implement plugin system for extensibility
- Add web deployment capabilities
- Optimize for large-scale datasets
- Build community and gather user feedback
