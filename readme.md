# 🔬 PySEE — Interactive, Reproducible Bioinformatics Visualization for Python

**PySEE** is an open-source project bringing **iSEE-style linked dashboards** to the **Python bioinformatics ecosystem**.  

If you use **AnnData / Scanpy / MuData / Zarr**, you know the struggle of wiring up UMAP plots, violin plots, QC panels, and genome browsers by hand. R has [Shiny](https://shiny.posit.co/) and [iSEE](https://bioconductor.org/packages/release/bioc/html/iSEE.html).  

👉 **PySEE fills that gap in Python**: a lightweight, notebook-first toolkit for **interactive exploration + reproducible code export**.

---

## ✨ Features

### MVP (v0.1, <1 month)
- **AnnData support** out of the box  
- **Two linked panels**:
  - UMAP/t-SNE embedding (scatter)
  - Gene expression violin/box plot
- **Linked selection**: brushing on UMAP updates violin
- **Reproducible code export**: selections → Python snippet
- **Notebook-first UX** (Jupyter/VS Code, no server setup needed)

### Future Directions
- 🔗 More panels: heatmaps, dotplots, QC metrics  
- 🧬 Genome browser panels (IGV / JBrowse)  
- 🧩 Spatial viewer (Vitessce) and imaging viewer (napari)  
- ☁️ Cloud-scale rendering (Datashader, Zarr-backed data)  
- 🎛️ Plugin system for custom panels  
- 🌍 Deployment as shareable web apps (FastAPI/Dash backend)  

---

## 🚀 Why PySEE?

- **Python-native**: integrates directly with AnnData, Scanpy, scvi-tools, PyTorch  
- **Linked & interactive**: selections propagate across panels  
- **Reproducible**: every UI action can export a Python snippet  
- **Complementary**: works alongside projects like [OLAF](https://arxiv.org/abs/2504.03976) (LLM-based bioinformatics) and [OLSA](https://github.com/openlifescience-ai) (AI benchmarks) as the **visual exploration layer**

---

## 📊 Quickstart (MVP sketch)

```python
import scanpy as sc
from pysee import PySEE
from pysee.panels import UMAPPanel, ViolinPanel

# Load toy data
adata = sc.datasets.pbmc3k()
sc.pp.pca(adata); sc.pp.neighbors(adata); sc.tl.umap(adata)

# Build dashboard
app = PySEE(adata)
app.add_panel("umap", UMAPPanel(color="louvain"))
app.add_panel("violin", ViolinPanel(gene="CD3D"))

# Link: UMAP selection filters violin
app.link(source="umap", target="violin")

# Show in notebook
app.show()

# Export reproducible code for current selection
print(app.export_code())
