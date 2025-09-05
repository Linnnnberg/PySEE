# DotPlot Panel Chart Summary

## 📊 Chart Overview
The DotPlot panel displays **marker gene expression** across different cell clusters in the PBMC3K dataset.

## 🧬 Genes Analyzed
- **MS4A1** - B cell marker
- **NKG7** - NK cell marker
- **GNLY** - NK cell marker
- **CD19** - B cell marker
- **FCGR3A** - Monocyte marker

## 👥 Cell Clusters
The data shows **6 distinct cell clusters** (0, 1, 2, 3, 4, 5) identified by Leiden clustering.

## 📈 Chart Features

### Dot Size
- Represents **percentage of cells** expressing each gene
- Larger dots = higher percentage of cells expressing the gene

### Color Intensity
- Represents **mean expression level** of each gene
- Darker colors = higher expression levels

### Interactive Features
- **Hover** over dots to see detailed statistics
- **Zoom** and **pan** to explore the data
- **Legend** shows color scale for expression levels

## 🔍 Key Insights

### Cluster 3 (B Cells)
- **High MS4A1 expression** - Strong B cell marker
- **High CD19 expression** - Confirms B cell identity
- **Low NKG7/GNLY** - No NK cell markers

### Cluster 2 (NK Cells)
- **High NKG7 expression** - Strong NK cell marker
- **High GNLY expression** - Confirms NK cell identity
- **Low MS4A1/CD19** - No B cell markers

### Other Clusters
- **Clusters 0, 1, 4, 5** show mixed or low expression patterns
- Likely represent **T cells, monocytes, or other cell types**

## 🎯 Chart Purpose
This DotPlot helps identify **cell type-specific gene expression patterns** and validate the quality of cell clustering in single-cell RNA sequencing data.

## 📁 Generated Files
- `dotplot_final.html` - Latest working version
- `dotplot_debug.html` - Version with debug information
- `dotplot_fixed.html` - Version after fixing Index error
- `dotplot_panel_demo.html` - Initial demo version
