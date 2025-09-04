"""
Quick preview of QC panel functionality and statistics.
"""

import numpy as np
import pandas as pd
import anndata as ad
from pysee import PySEE
from pysee.panels.qc import QCPanel


def create_preview_data():
    """Create a small dataset for quick preview."""
    n_cells = 500
    n_genes = 1000
    
    # Create expression matrix
    np.random.seed(42)
    expression_matrix = np.random.negative_binomial(3, 0.4, size=(n_cells, n_genes)).astype(float)
    
    # Add mitochondrial genes
    n_mito_genes = int(0.1 * n_genes)
    mito_gene_indices = np.random.choice(n_genes, n_mito_genes, replace=False)
    expression_matrix[:, mito_gene_indices] *= 2.0
    
    # Create gene names
    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
    for i, idx in enumerate(mito_gene_indices):
        gene_names[idx] = f"MT-Gene_{i:03d}"
    
    # Create cell names
    cell_names = [f"Cell_{i:04d}" for i in range(n_cells)]
    
    # Calculate QC metrics
    mito_genes = [name.startswith('MT-') for name in gene_names]
    mito_counts = expression_matrix[:, mito_genes].sum(axis=1)
    total_counts = expression_matrix.sum(axis=1)
    mito_percent = mito_counts / total_counts * 100
    
    # Create metadata
    obs_data = {
        'cell_type': np.random.choice(['T_cell', 'B_cell', 'NK_cell'], n_cells),
        'total_counts': total_counts,
        'detected_genes': (expression_matrix > 0).sum(axis=1),
        'mito_percent': mito_percent,
    }
    
    # Create AnnData
    adata = ad.AnnData(
        X=expression_matrix,
        obs=pd.DataFrame(obs_data, index=cell_names),
        var=pd.DataFrame(index=gene_names)
    )
    
    return adata


def show_qc_preview():
    """Show a preview of QC panel functionality."""
    print("🔍 QC Panel Preview")
    print("=" * 40)
    
    # Create data
    adata = create_preview_data()
    
    # Initialize PySEE and QC panel
    app = PySEE(adata)
    qc_panel = QCPanel("qc", title="QC Metrics Preview")
    app.add_panel("qc", qc_panel)
    
    # Calculate and show QC statistics
    metrics = qc_panel._calculate_qc_metrics()
    
    print("📊 QC Metrics Summary:")
    print(f"   Dataset: {adata.n_obs} cells, {adata.n_vars} genes")
    print()
    
    if 'mito_percent' in metrics:
        mito_stats = metrics['mito_percent']
        print("🧬 Mitochondrial Gene Percentage:")
        print(f"   Mean: {mito_stats.mean():.1f}%")
        print(f"   Median: {np.median(mito_stats):.1f}%")
        print(f"   Range: {mito_stats.min():.1f}% - {mito_stats.max():.1f}%")
        print(f"   Cells >20%: {(mito_stats > 20).sum()} ({(mito_stats > 20).mean()*100:.1f}%)")
        print()
    
    if 'total_counts' in metrics:
        count_stats = metrics['total_counts']
        print("📈 Total Gene Counts:")
        print(f"   Mean: {count_stats.mean():.0f}")
        print(f"   Median: {np.median(count_stats):.0f}")
        print(f"   Range: {count_stats.min():.0f} - {count_stats.max():.0f}")
        print()
    
    if 'detected_genes' in metrics:
        gene_stats = metrics['detected_genes']
        print("🔬 Detected Genes per Cell:")
        print(f"   Mean: {gene_stats.mean():.0f}")
        print(f"   Median: {np.median(gene_stats):.0f}")
        print(f"   Range: {gene_stats.min():.0f} - {gene_stats.max():.0f}")
        print()
    
    # Show filtering thresholds
    print("🎯 Current Filtering Thresholds:")
    print(f"   Mitochondrial %: < {qc_panel.get_config('mito_threshold')}%")
    print(f"   Total counts: {qc_panel.get_config('min_counts')} - {qc_panel.get_config('max_counts')}")
    print(f"   Detected genes: {qc_panel.get_config('min_genes')} - {qc_panel.get_config('max_genes')}")
    print()
    
    # Show code export preview
    print("💻 Generated Filtering Code Preview:")
    code = qc_panel.get_selection_code()
    code_lines = code.split('\n')[:10]  # Show first 10 lines
    for line in code_lines:
        print(f"   {line}")
    print("   ...")
    print()
    
    print("✅ QC Panel is ready for interactive visualization!")
    print("   Open the HTML files in your browser to see the full interactive charts.")


if __name__ == "__main__":
    show_qc_preview()
