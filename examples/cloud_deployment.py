"""
Cloud deployment examples for PySEE.

This script demonstrates how to use PySEE in cloud environments
for large dataset analysis.
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import anndata as ad
import scanpy as sc

from pysee import PySEE
from pysee.panels.heatmap import HeatmapPanel
from pysee.panels.qc import QCPanel
from pysee.panels.umap import UMAPPanel
from pysee.panels.violin import ViolinPanel


def google_colab_example():
    """Example for Google Colab deployment."""
    print("🚀 Google Colab Deployment Example")
    print("=" * 50)

    print(
        """
# Google Colab Setup
!pip install pysee scanpy

import scanpy as sc
from pysee import PySEE
from pysee.panels.umap import UMAPPanel
from pysee.panels.violin import ViolinPanel
from pysee.panels.heatmap import HeatmapPanel
from pysee.panels.qc import QCPanel

# Load dataset (works well in Colab's 12-25 GB RAM)
adata = sc.datasets.pbmc68k_reduced()  # 68K cells, 8.5 GB
print(f"Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

# Create PySEE app
app = PySEE(adata, title="PBMC 68K Analysis")

# Add panels
qc_panel = QCPanel("qc", title="Quality Control")
app.add_panel("qc", qc_panel)

violin_panel = ViolinPanel("violin", gene="MS4A1", title="B Cell Marker")
app.add_panel("violin", violin_panel)

heatmap_panel = HeatmapPanel("heatmap", title="Gene Expression")
app.add_panel("heatmap", heatmap_panel)

# Render panels
qc_fig = app.render_panel("qc")
violin_fig = app.render_panel("violin")
heatmap_fig = app.render_panel("heatmap")

# Display results
qc_fig.show()
violin_fig.show()
heatmap_fig.show()
    """
    )


def aws_ec2_example():
    """Example for AWS EC2 deployment."""
    print("\n🚀 AWS EC2 Deployment Example")
    print("=" * 50)

    print(
        """
# AWS EC2 Setup (t3.xlarge or larger)
# 1. Launch EC2 instance with 16+ GB RAM
# 2. Install Python and dependencies

# Install PySEE
pip install pysee scanpy

# Load large dataset
import scanpy as sc
from pysee import PySEE

# For very large datasets, use backed mode
adata = ad.read_h5ad('large_dataset.h5ad', backed='r')
print(f"Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

# Create PySEE app
app = PySEE(adata, title="Large Dataset Analysis")

# Add panels (memory efficient)
qc_panel = QCPanel("qc", title="Quality Control")
app.add_panel("qc", qc_panel)

# Render and save results
qc_fig = app.render_panel("qc")
qc_fig.write_html("qc_results.html")

print("Analysis complete! Results saved to qc_results.html")
    """
    )


def local_server_example():
    """Example for local server deployment."""
    print("\n🚀 Local Server Deployment Example")
    print("=" * 50)

    print(
        """
# Local Server Setup (32+ GB RAM recommended)
# 1. Install PySEE on server
# 2. Use Jupyter Lab or similar

import scanpy as sc
from pysee import PySEE
import anndata as ad

# Load large dataset
adata = sc.datasets.pbmc68k_reduced()  # 68K cells
print(f"Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

# Create comprehensive analysis
app = PySEE(adata, title="Comprehensive Analysis")

# Add all panels
qc_panel = QCPanel("qc", title="Quality Control")
app.add_panel("qc", qc_panel)

violin_panel = ViolinPanel("violin", gene="MS4A1", title="B Cell Marker")
app.add_panel("violin", violin_panel)

heatmap_panel = HeatmapPanel("heatmap", title="Gene Expression")
app.add_panel("heatmap", heatmap_panel)

# Render all panels
results = {}
for panel_id in app.panels.keys():
    results[panel_id] = app.render_panel(panel_id)
    results[panel_id].write_html(f"{panel_id}_results.html")

print("All analyses complete! Results saved to HTML files.")
    """
    )


def memory_efficient_cloud_example():
    """Example for memory-efficient cloud usage."""
    print("\n🚀 Memory-Efficient Cloud Example")
    print("=" * 50)

    print(
        """
# Memory-Efficient Cloud Usage
# For datasets that are still too large for cloud RAM

import scanpy as sc
from pysee import PySEE
import anndata as ad
import numpy as np

# Load dataset in backed mode
adata = ad.read_h5ad('very_large_dataset.h5ad', backed='r')
print(f"Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

# Subsample for analysis
n_subsample = 50000  # Limit for cloud RAM
if adata.n_obs > n_subsample:
    indices = np.random.choice(adata.n_obs, n_subsample, replace=False)
    adata_sub = adata[indices].copy()
    print(f"Subsampled to: {adata_sub.n_obs:,} cells")
else:
    adata_sub = adata.copy()

# Create analysis
app = PySEE(adata_sub, title="Subsampled Analysis")

# Add panels
qc_panel = QCPanel("qc", title="Quality Control")
app.add_panel("qc", qc_panel)

# Render
qc_fig = app.render_panel("qc")
qc_fig.write_html("subsampled_qc_results.html")

print("Subsampled analysis complete!")
    """
    )


def main():
    """Run all cloud deployment examples."""
    print("☁️ PySEE Cloud Deployment Examples")
    print("=" * 60)

    google_colab_example()
    aws_ec2_example()
    local_server_example()
    memory_efficient_cloud_example()

    print("\n💡 Cloud Deployment Summary:")
    print("=" * 60)
    print("✅ Google Colab: Perfect for medium datasets (68K cells)")
    print("✅ AWS EC2: Great for large datasets (100K+ cells)")
    print("✅ Local Server: Optimal for very large datasets (1M+ cells)")
    print("✅ Memory-Efficient: Use subsampling and backed mode")

    print("\n🎯 Recommendation:")
    print("For large datasets, use cloud instead of complex local memory strategies!")


if __name__ == "__main__":
    main()
