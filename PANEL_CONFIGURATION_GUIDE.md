# PySEE Panel Configuration Guide

This guide explains how to configure PySEE panels for optimal visualization and publication-ready export.

## Overview

PySEE panels can be configured in multiple ways:
1. **During Creation**: Set initial parameters when creating panels
2. **Runtime Configuration**: Modify settings after panel creation
3. **Export Configuration**: Override settings specifically for export
4. **Template Configuration**: Use journal-specific templates

## Panel Types and Configuration Options

### 1. UMAP Panel (`UMAPPanel`)

**Purpose**: Visualize dimensionality reduction embeddings (UMAP, t-SNE, PCA)

**Configuration Options**:
```python
# Create with configuration
umap_panel = UMAPPanel(
    panel_id="umap",
    embedding="X_umap",           # Key in adata.obsm
    color="leiden",               # Column in adata.obs for coloring
    title="UMAP Visualization"
)

# Runtime configuration
umap_panel.set_config("point_size", 5)        # Point size (default: 3)
umap_panel.set_config("opacity", 0.8)         # Point opacity (default: 0.7)
umap_panel.set_config("show_legend", True)    # Show legend (default: True)
umap_panel.set_config("color_scale", "viridis") # Color scale (default: "viridis")
```

**Available Settings**:
- `embedding`: Embedding key in `adata.obsm` (default: "X_umap")
- `color`: Column name in `adata.obs` for coloring (default: None)
- `point_size`: Size of points (default: 3)
- `opacity`: Point opacity 0-1 (default: 0.7)
- `show_legend`: Show color legend (default: True)
- `color_scale`: Color scale name (default: "viridis")

### 2. Violin Panel (`ViolinPanel`)

**Purpose**: Statistical plots for gene expression analysis

**Configuration Options**:
```python
# Create with configuration
violin_panel = ViolinPanel(
    panel_id="violin",
    genes=["CD3D", "CD3E", "CD8A"],  # List of genes to plot
    groupby="leiden",                # Grouping variable
    title="Gene Expression"
)

# Runtime configuration
violin_panel.set_config("plot_type", "violin")    # "violin", "box", "strip"
violin_panel.set_config("show_points", True)      # Show individual points
violin_panel.set_config("show_median", True)      # Show median line
violin_panel.set_config("jitter", 0.3)            # Jitter amount for points
```

**Available Settings**:
- `genes`: List of genes to plot (required)
- `groupby`: Column in `adata.obs` for grouping (default: None)
- `plot_type`: Type of plot ("violin", "box", "strip") (default: "violin")
- `show_points`: Show individual points (default: True)
- `show_median`: Show median line (default: True)
- `jitter`: Jitter amount for points (default: 0.3)

### 3. Heatmap Panel (`HeatmapPanel`)

**Purpose**: Interactive heatmaps for gene expression matrices

**Configuration Options**:
```python
# Create with configuration
heatmap_panel = HeatmapPanel(
    panel_id="heatmap",
    genes=["CD3D", "CD3E", "CD8A", "CD4"],
    groupby="leiden",
    title="Gene Expression Heatmap"
)

# Runtime configuration
heatmap_panel.set_config("cluster_genes", True)    # Cluster genes (default: True)
heatmap_panel.set_config("cluster_cells", True)    # Cluster cells (default: True)
heatmap_panel.set_config("color_scale", "RdBu_r")  # Color scale (default: "RdBu_r")
heatmap_panel.set_config("z_score", 1)             # Z-score normalization (default: 1)
```

**Available Settings**:
- `genes`: List of genes to plot (required)
- `groupby`: Column in `adata.obs` for grouping (default: None)
- `cluster_genes`: Cluster genes (default: True)
- `cluster_cells`: Cluster cells (default: True)
- `color_scale`: Color scale name (default: "RdBu_r")
- `z_score`: Z-score normalization (0=no, 1=genes, 2=cells) (default: 1)

### 4. QC Panel (`QCPanel`)

**Purpose**: Quality control visualizations for data assessment

**Configuration Options**:
```python
# Create with configuration
qc_panel = QCPanel(
    panel_id="qc",
    metrics=["n_genes", "n_counts", "pct_counts_mt"],
    title="Quality Control"
)

# Runtime configuration
qc_panel.set_config("show_thresholds", True)      # Show filtering thresholds
qc_panel.set_config("threshold_n_genes", 200)     # Threshold for n_genes
qc_panel.set_config("threshold_n_counts", 1000)   # Threshold for n_counts
qc_panel.set_config("threshold_pct_mt", 20)       # Threshold for pct_counts_mt
```

**Available Settings**:
- `metrics`: List of QC metrics to plot (default: ["n_genes", "n_counts", "pct_counts_mt"])
- `show_thresholds`: Show filtering thresholds (default: True)
- `threshold_n_genes`: Threshold for n_genes (default: 200)
- `threshold_n_counts`: Threshold for n_counts (default: 1000)
- `threshold_pct_mt`: Threshold for pct_counts_mt (default: 20)

### 5. Dot Plot Panel (`DotPlotPanel`)

**Purpose**: Dot plots for marker gene expression analysis

**Configuration Options**:
```python
# Create with configuration
dotplot_panel = DotPlotPanel(
    panel_id="dotplot",
    genes=["CD3D", "CD3E", "CD8A", "CD4"],
    groupby="leiden",
    title="Marker Gene Expression"
)

# Runtime configuration
dotplot_panel.set_config("dot_size", "pct_expressed")  # Size by "pct_expressed" or "mean_expr"
dotplot_panel.set_config("dot_color", "mean_expr")     # Color by "mean_expr" or "pct_expressed"
dotplot_panel.set_config("show_colorbar", True)        # Show color bar (default: True)
dotplot_panel.set_config("color_scale", "viridis")     # Color scale (default: "viridis")
```

**Available Settings**:
- `genes`: List of genes to plot (required)
- `groupby`: Column in `adata.obs` for grouping (default: None)
- `dot_size`: Size mapping ("pct_expressed" or "mean_expr") (default: "pct_expressed")
- `dot_color`: Color mapping ("mean_expr" or "pct_expressed") (default: "mean_expr")
- `show_colorbar`: Show color bar (default: True)
- `color_scale`: Color scale name (default: "viridis")

## Configuration Methods

### 1. Individual Configuration

```python
# Set single configuration value
panel.set_config("point_size", 5)

# Set multiple configuration values
panel.update_config({
    "point_size": 5,
    "opacity": 0.8,
    "show_legend": True
})

# Get configuration value
point_size = panel.get_config("point_size", default=3)
```

### 2. Export-Specific Configuration

```python
# Export with custom configuration
panel.export(
    "output.png",
    format="png",
    template="nature",
    width=1800,      # Override template width
    height=1200,     # Override template height
    dpi=300,         # Override template DPI
    title="Custom Title",
    caption="Custom caption"
)
```

### 3. Template Configuration

```python
# Use journal-specific templates
panel.export("output.png", template="nature")    # Nature journal style
panel.export("output.png", template="science")   # Science journal style
panel.export("output.png", template="cell")      # Cell journal style
panel.export("output.png", template="custom")    # Custom style (default)
```

## Journal Templates

### Nature Template
- **Width**: 180mm
- **Height**: 120mm
- **DPI**: 300
- **Font**: Arial, 8pt
- **Line Width**: 0.5pt
- **Marker Size**: 3pt

### Science Template
- **Width**: 180mm
- **Height**: 120mm
- **DPI**: 300
- **Font**: Arial, 7pt
- **Line Width**: 0.5pt
- **Marker Size**: 2.5pt

### Cell Template
- **Width**: 180mm
- **Height**: 120mm
- **DPI**: 300
- **Font**: Arial, 8pt
- **Line Width**: 0.5pt
- **Marker Size**: 3pt

## Best Practices

### 1. For Publication Figures
```python
# Use journal-specific templates
panel.export("figure.png", template="nature", dpi=300)

# Add descriptive titles and captions
panel.export(
    "figure.png",
    template="nature",
    title="UMAP Visualization of Single-Cell Data",
    caption="Cells are colored by Leiden clusters. Each point represents a single cell."
)
```

### 2. For Interactive Exploration
```python
# Configure for better visibility
panel.set_config("point_size", 4)
panel.set_config("opacity", 0.8)
panel.set_config("show_legend", True)
```

### 3. For Batch Export
```python
# Configure all panels consistently
for panel in app.panels.values():
    panel.set_config("show_legend", True)
    panel.set_config("font_size", 10)

# Export all panels
app.export_all_panels("output_dir", template="nature")
```

## Troubleshooting

### Common Issues

1. **Panel not rendering**: Check if data requirements are met
   ```python
   if not panel.validate_data():
       print("Panel data requirements not met")
   ```

2. **Export fails**: Check file permissions and directory existence
   ```python
   import os
   os.makedirs("output_dir", exist_ok=True)
   ```

3. **Low resolution**: Increase DPI and dimensions
   ```python
   panel.export("output.png", dpi=600, width=3600, height=2400)
   ```

### Getting Help

- Check panel information: `panel.get_panel_info()`
- Check publication info: `panel.get_publication_info()`
- Check export capabilities: `app.get_export_info()`

## Examples

### Complete Example
```python
import scanpy as sc
from pysee import PySEE
from pysee.panels import UMAPPanel, ViolinPanel, QCPanel

# Load data
adata = sc.datasets.pbmc3k()
# ... preprocessing ...

# Create dashboard
app = PySEE(adata, title="My Analysis")

# Add and configure panels
umap_panel = UMAPPanel("umap", color="leiden")
umap_panel.set_config("point_size", 4)
umap_panel.set_config("opacity", 0.8)
app.add_panel(umap_panel)

violin_panel = ViolinPanel("violin", genes=["CD3D", "CD3E"])
violin_panel.set_config("plot_type", "violin")
violin_panel.set_config("show_points", True)
app.add_panel(violin_panel)

# Export with custom settings
umap_panel.export(
    "umap_nature.png",
    template="nature",
    title="UMAP Visualization",
    caption="Single-cell UMAP embedding"
)

# Batch export
app.export_all_panels("publication_figures", template="nature")
```
