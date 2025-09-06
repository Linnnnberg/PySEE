# DotPlot Panel Guide

## 📊 What is a DotPlot?

A DotPlot is a visualization that shows **gene expression patterns** across different cell groups. It helps identify which genes are expressed in which cell types.

## 🎨 Visual Elements

### **Circle Size (Dot Size)**
- **Bigger circles** = **Higher percentage** of cells express this gene
- **Smaller circles** = **Lower percentage** of cells express this gene
- **Range**: 0% to 100% of cells in that group

### **Color Intensity**
- **Dark red** = **High gene expression** (gene is very active)
- **Yellow** = **Medium gene expression** (gene is moderately active)
- **Blue/White** = **Low gene expression** (gene is barely active)

## 🔍 How to Read the Chart

### **Rows (Y-axis)**
- Each row represents a **gene**
- Genes are listed vertically

### **Columns (X-axis)**
- Each column represents a **cell group/cluster**
- Groups are identified by numbers (0, 1, 2, 3, etc.)

### **Each Dot**
- **Position**: Intersection of gene (row) and group (column)
- **Size**: How many cells in that group express this gene
- **Color**: How much this gene is expressed in that group

## 📈 Interpreting Results

### **Cell Type Identification**
- **Big red dots** = **Strong markers** for that cell type
- **Small blue dots** = **Not markers** for that cell type

### **Gene Specificity**
- **Gene with big dots in one group only** = **Specific marker** for that cell type
- **Gene with big dots in multiple groups** = **Common gene** across cell types

### **Data Quality**
- **Clear patterns** = Good cell type separation
- **Mixed patterns** = Mixed cell populations or poor clustering

## 🧬 Example Interpretation

### **B Cell Markers (Cluster 3)**
- **MS4A1**: Big red dot = Most B cells express this gene strongly
- **CD19**: Big red dot = Most B cells express this gene strongly
- **NKG7**: Small blue dot = B cells don't express this NK marker

### **NK Cell Markers (Cluster 2)**
- **NKG7**: Big red dot = Most NK cells express this gene strongly
- **GNLY**: Big red dot = Most NK cells express this gene strongly
- **MS4A1**: Small blue dot = NK cells don't express this B cell marker

## 🎯 Key Questions to Ask

1. **Which genes are specific to which cell types?**
   - Look for genes with big dots in only one group

2. **How pure are my cell clusters?**
   - Good clusters have clear marker gene patterns

3. **Which genes are commonly expressed?**
   - Look for genes with big dots across multiple groups

4. **Are there any unexpected expression patterns?**
   - Look for genes with big dots in unexpected groups

## 🔧 Interactive Features

- **Hover** over dots to see exact values
- **Zoom** and **pan** to explore the data
- **Color bar** shows the expression level scale
- **Responsive design** works on different screen sizes

## 📊 Data Requirements

- **Genes**: List of genes to analyze
- **Groups**: Cell clusters or categories
- **Expression data**: Gene expression matrix
- **Metadata**: Group assignments for each cell

## 🎨 Customization Options

- **Color scale**: Choose different color schemes
- **Dot size range**: Adjust minimum and maximum dot sizes
- **Filtering**: Set minimum expression and percentage thresholds
- **Sorting**: Sort genes or groups by expression levels

## 💡 Tips for Analysis

1. **Start with known marker genes** to validate your clusters
2. **Look for clear patterns** rather than individual dots
3. **Compare relative sizes** within the same gene across groups
4. **Use hover information** for exact numerical values
5. **Consider both size and color** for complete interpretation

## 🚀 Getting Started

```python
from pysee.panels.dotplot import DotPlotPanel

# Create DotPlot panel
dotplot = DotPlotPanel(
    panel_id="my_dotplot",
    genes=["GENE1", "GENE2", "GENE3"],
    group_by="cell_type",
    title="My Gene Expression Analysis"
)

# Add to dashboard
app.add_panel("dotplot", dotplot)
```

---

**The DotPlot panel helps you understand which genes define which cell types! 🧬**
