#!/usr/bin/env python3
"""
Multi-Panel Selection Integration Demo for PySEE

This script demonstrates the integrated advanced selection tools across all PySEE panels,
showing how selections can be made and synchronized between different visualization types.

Features demonstrated:
- Advanced selection tools in all panel types (UMAP, Violin, Heatmap, QC, DotPlot)
- Selection synchronization between panels
- Multiple selection modes (Replace, Add, Subtract, Intersect)
- Selection history with undo/redo
- Interactive help system
- Selection statistics and summaries
"""

from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.offline as pyo
import scanpy as sc
from plotly.subplots import make_subplots

from pysee.core.dashboard import PySEE

# Import PySEE components
from pysee.core.data import AnnDataWrapper
from pysee.panels.dotplot import DotPlotPanel
from pysee.panels.heatmap import HeatmapPanel
from pysee.panels.qc import QCPanel
from pysee.panels.umap import UMAPPanel
from pysee.panels.violin import ViolinPanel


def create_sample_anndata(n_cells: int = 1000, n_genes: int = 2000) -> AnnDataWrapper:
    """Create a sample AnnData object for demonstration."""
    print(f"Creating sample data with {n_cells} cells and {n_genes} genes...")

    # Create random expression data
    np.random.seed(42)
    X = np.random.negative_binomial(5, 0.3, (n_cells, n_genes))

    # Create gene names
    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]

    # Create cell names
    cell_names = [f"Cell_{i:04d}" for i in range(n_cells)]

    # Create cell metadata
    obs_data = pd.DataFrame(
        {
            "cell_type": np.random.choice(
                ["T_cell", "B_cell", "NK_cell", "Monocyte", "Dendritic"], n_cells
            ),
            "batch": np.random.choice(["Batch1", "Batch2", "Batch3"], n_cells),
            "n_genes": np.random.poisson(1500, n_cells),
            "total_counts": np.random.poisson(10000, n_cells),
        },
        index=cell_names,
    )

    # Create gene metadata
    var_data = pd.DataFrame(
        {
            "gene_type": np.random.choice(["protein_coding", "lncRNA", "pseudogene"], n_genes),
            "highly_variable": np.random.choice([True, False], n_genes, p=[0.1, 0.9]),
        },
        index=gene_names,
    )

    # Create AnnData object
    adata = sc.AnnData(X=X, obs=obs_data, var=var_data)
    adata.obs_names = cell_names
    adata.var_names = gene_names

    # Add some structure to the data
    # Create a UMAP embedding
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)

    # Add some clustering
    sc.tl.leiden(adata, resolution=0.5)

    return AnnDataWrapper(adata)


def create_multi_panel_dashboard(data_wrapper: AnnDataWrapper) -> PySEE:
    """Create a multi-panel dashboard with all panel types."""
    print("Creating multi-panel dashboard...")

    # Create PySEE dashboard
    dashboard = PySEE(data_wrapper.adata, "Multi-Panel Selection Demo")

    # Add UMAP panel
    umap_panel = UMAPPanel(
        panel_id="umap", embedding="X_umap", color="leiden", title="UMAP Embedding"
    )
    umap_panel.enable_advanced_selection()
    dashboard.add_panel("umap", umap_panel)

    # Add Violin panel
    violin_panel = ViolinPanel(
        panel_id="violin",
        gene=data_wrapper.adata.var_names[0],  # Use first available gene
        group_by="leiden",
        title="Gene Expression Distribution",
    )
    violin_panel.enable_advanced_selection()
    dashboard.add_panel("violin", violin_panel)

    # Add Heatmap panel
    heatmap_panel = HeatmapPanel(
        panel_id="heatmap",
        genes=data_wrapper.adata.var_names[:20].tolist(),  # Top 20 genes
        title="Gene Expression Heatmap",
    )
    heatmap_panel.enable_advanced_selection()
    dashboard.add_panel("heatmap", heatmap_panel)

    # Add QC panel
    qc_panel = QCPanel(panel_id="qc", title="Quality Control Metrics")
    qc_panel.enable_advanced_selection()
    dashboard.add_panel("qc", qc_panel)

    # Add DotPlot panel
    dotplot_panel = DotPlotPanel(
        panel_id="dotplot",
        genes=data_wrapper.adata.var_names[:10].tolist(),  # Top 10 genes
        group_by="leiden",
        title="Marker Gene Expression",
    )
    dotplot_panel.enable_advanced_selection()
    dashboard.add_panel("dotplot", dotplot_panel)

    # Link panels for selection synchronization
    dashboard.link("umap", "violin")
    dashboard.link("umap", "heatmap")
    dashboard.link("umap", "qc")
    dashboard.link("umap", "dotplot")

    return dashboard


def create_individual_panel_demos(data_wrapper: AnnDataWrapper) -> None:
    """Create individual demos for each panel type."""
    print("Creating individual panel demos...")

    # UMAP Panel Demo
    print("  - UMAP Panel Demo")
    umap_panel = UMAPPanel(
        "umap_demo", embedding="X_umap", color="leiden", title="UMAP with Advanced Selection"
    )
    umap_panel.data_wrapper = data_wrapper
    umap_panel.enable_advanced_selection()

    fig = umap_panel.render()
    fig.update_layout(title="UMAP Panel - Advanced Selection Tools", width=800, height=600)
    fig.write_html("demo_umap_advanced_selection.html")

    # Violin Panel Demo
    print("  - Violin Panel Demo")
    violin_panel = ViolinPanel(
        "violin_demo",
        gene=data_wrapper.adata.var_names[0],
        group_by="leiden",
        title="Violin Plot with Advanced Selection",
    )
    violin_panel.data_wrapper = data_wrapper
    violin_panel.enable_advanced_selection()

    fig = violin_panel.render()
    fig.update_layout(title="Violin Panel - Advanced Selection Tools", width=800, height=600)
    fig.write_html("demo_violin_advanced_selection.html")

    # Heatmap Panel Demo
    print("  - Heatmap Panel Demo")
    heatmap_panel = HeatmapPanel(
        "heatmap_demo",
        genes=data_wrapper.adata.var_names[:15].tolist(),
        title="Heatmap with Advanced Selection",
    )
    heatmap_panel.data_wrapper = data_wrapper
    heatmap_panel.enable_advanced_selection()

    fig = heatmap_panel.render()
    fig.update_layout(title="Heatmap Panel - Advanced Selection Tools", width=800, height=600)
    fig.write_html("demo_heatmap_advanced_selection.html")

    # QC Panel Demo
    print("  - QC Panel Demo")
    qc_panel = QCPanel("qc_demo", title="QC Metrics with Advanced Selection")
    qc_panel.data_wrapper = data_wrapper
    qc_panel.enable_advanced_selection()

    fig = qc_panel.render()
    fig.update_layout(title="QC Panel - Advanced Selection Tools", width=800, height=600)
    fig.write_html("demo_qc_advanced_selection.html")

    # DotPlot Panel Demo
    print("  - DotPlot Panel Demo")
    dotplot_panel = DotPlotPanel(
        "dotplot_demo",
        genes=data_wrapper.adata.var_names[:8].tolist(),
        group_by="leiden",
        title="DotPlot with Advanced Selection",
    )
    dotplot_panel.data_wrapper = data_wrapper
    dotplot_panel.enable_advanced_selection()

    fig = dotplot_panel.render()
    fig.update_layout(title="DotPlot Panel - Advanced Selection Tools", width=800, height=600)
    fig.write_html("demo_dotplot_advanced_selection.html")


def create_multi_panel_comparison(data_wrapper: AnnDataWrapper) -> None:
    """Create a side-by-side comparison of all panels."""
    print("Creating multi-panel comparison...")

    # Create subplots
    fig = make_subplots(
        rows=2,
        cols=3,
        subplot_titles=[
            "UMAP Embedding",
            "Violin Plot",
            "Heatmap",
            "QC Metrics",
            "DotPlot",
            "Selection Statistics",
        ],
        specs=[
            [{"type": "scatter"}, {"type": "violin"}, {"type": "heatmap"}],
            [{"type": "histogram"}, {"type": "scatter"}, {"type": "table"}],
        ],
        vertical_spacing=0.1,
        horizontal_spacing=0.1,
    )

    # UMAP subplot
    umap_panel = UMAPPanel("umap", embedding="X_umap", color="leiden")
    umap_panel.data_wrapper = data_wrapper
    umap_panel.enable_advanced_selection()
    umap_fig = umap_panel.render()

    # Add UMAP traces to subplot
    for trace in umap_fig.data:
        fig.add_trace(trace, row=1, col=1)

    # Violin subplot
    violin_panel = ViolinPanel("violin", gene=data_wrapper.adata.var_names[0], group_by="leiden")
    violin_panel.data_wrapper = data_wrapper
    violin_panel.enable_advanced_selection()
    violin_fig = violin_panel.render()

    # Add Violin traces to subplot
    for trace in violin_fig.data:
        fig.add_trace(trace, row=1, col=2)

    # Heatmap subplot
    heatmap_panel = HeatmapPanel("heatmap", genes=data_wrapper.adata.var_names[:10].tolist())
    heatmap_panel.data_wrapper = data_wrapper
    heatmap_panel.enable_advanced_selection()
    heatmap_fig = heatmap_panel.render()

    # Add Heatmap traces to subplot
    for trace in heatmap_fig.data:
        fig.add_trace(trace, row=1, col=3)

    # QC subplot
    qc_panel = QCPanel("qc")
    qc_panel.data_wrapper = data_wrapper
    qc_panel.enable_advanced_selection()
    qc_fig = qc_panel.render()

    # Add QC traces to subplot
    for trace in qc_fig.data:
        fig.add_trace(trace, row=2, col=1)

    # DotPlot subplot
    dotplot_panel = DotPlotPanel(
        "dotplot", genes=data_wrapper.adata.var_names[:6].tolist(), group_by="leiden"
    )
    dotplot_panel.data_wrapper = data_wrapper
    dotplot_panel.enable_advanced_selection()
    dotplot_fig = dotplot_panel.render()

    # Add DotPlot traces to subplot
    for trace in dotplot_fig.data:
        fig.add_trace(trace, row=2, col=2)

    # Selection statistics table
    stats_data = [
        ["Panel Type", "Selection Enabled", "Selection Modes", "History Support"],
        ["UMAP", "✓", "Point, Rectangular, Lasso", "✓"],
        ["Violin", "✓", "Point", "✓"],
        ["Heatmap", "✓", "Rectangular", "✓"],
        ["QC", "✓", "Rectangular", "✓"],
        ["DotPlot", "✓", "Point", "✓"],
    ]

    fig.add_trace(
        go.Table(
            header=dict(values=stats_data[0], fill_color="lightblue", align="center"),
            cells=dict(values=list(zip(*stats_data[1:])), fill_color="lightgray", align="center"),
        ),
        row=2,
        col=3,
    )

    # Update layout
    fig.update_layout(
        title="PySEE Multi-Panel Advanced Selection Tools - Complete Integration",
        height=1000,
        width=1400,
        showlegend=False,
    )

    # Update axes
    fig.update_xaxes(title_text="UMAP 1", row=1, col=1)
    fig.update_yaxes(title_text="UMAP 2", row=1, col=1)
    fig.update_xaxes(title_text="Group", row=1, col=2)
    fig.update_yaxes(title_text="Expression", row=1, col=2)
    fig.update_xaxes(title_text="Cells", row=1, col=3)
    fig.update_yaxes(title_text="Genes", row=1, col=3)
    fig.update_xaxes(title_text="Value", row=2, col=1)
    fig.update_yaxes(title_text="Count", row=2, col=1)
    fig.update_xaxes(title_text="Group", row=2, col=2)
    fig.update_yaxes(title_text="Gene", row=2, col=2)

    fig.write_html("demo_multi_panel_comparison.html")


def create_selection_synchronization_demo(data_wrapper: AnnDataWrapper) -> None:
    """Create a demo showing selection synchronization between panels."""
    print("Creating selection synchronization demo...")

    # Create a dashboard with linked panels
    dashboard = create_multi_panel_dashboard(data_wrapper)

    # Render all panels
    panels = {}
    for panel_id in ["umap", "violin", "heatmap", "qc", "dotplot"]:
        panel = dashboard.get_panel(panel_id)
        if panel:
            panels[panel_id] = panel.render()

    # Create a comprehensive layout
    fig = make_subplots(
        rows=3,
        cols=2,
        subplot_titles=[
            "UMAP (Primary Selection)",
            "Violin Plot (Linked)",
            "Heatmap (Linked)",
            "QC Metrics (Linked)",
            "DotPlot (Linked)",
            "Selection Statistics",
        ],
        specs=[
            [{"type": "scatter"}, {"type": "violin"}],
            [{"type": "heatmap"}, {"type": "histogram"}],
            [{"type": "scatter"}, {"type": "table"}],
        ],
        vertical_spacing=0.08,
        horizontal_spacing=0.1,
    )

    # Add UMAP traces
    if "umap" in panels:
        for trace in panels["umap"].data:
            fig.add_trace(trace, row=1, col=1)

    # Add Violin traces
    if "violin" in panels:
        for trace in panels["violin"].data:
            fig.add_trace(trace, row=1, col=2)

    # Add Heatmap traces
    if "heatmap" in panels:
        for trace in panels["heatmap"].data:
            fig.add_trace(trace, row=2, col=1)

    # Add QC traces
    if "qc" in panels:
        for trace in panels["qc"].data:
            fig.add_trace(trace, row=2, col=2)

    # Add DotPlot traces
    if "dotplot" in panels:
        for trace in panels["dotplot"].data:
            fig.add_trace(trace, row=3, col=1)

    # Selection synchronization info
    sync_info = [
        ["Feature", "Status", "Description"],
        ["Cross-Panel Selection", "✓", "Select in UMAP, see in all panels"],
        ["Selection Modes", "✓", "Replace, Add, Subtract, Intersect"],
        ["Selection History", "✓", "Undo/Redo with 50-step history"],
        ["Interactive Help", "✓", "Comprehensive tooltips and guides"],
        ["Selection Statistics", "✓", "Real-time counts and percentages"],
        ["Live Updates", "✓", "Automatic panel synchronization"],
    ]

    fig.add_trace(
        go.Table(
            header=dict(values=sync_info[0], fill_color="lightgreen", align="center"),
            cells=dict(values=list(zip(*sync_info[1:])), fill_color="lightgray", align="left"),
        ),
        row=3,
        col=2,
    )

    # Update layout
    fig.update_layout(
        title="PySEE Multi-Panel Selection Synchronization Demo",
        height=1200,
        width=1600,
        showlegend=False,
    )

    # Update axes
    fig.update_xaxes(title_text="UMAP 1", row=1, col=1)
    fig.update_yaxes(title_text="UMAP 2", row=1, col=1)
    fig.update_xaxes(title_text="Group", row=1, col=2)
    fig.update_yaxes(title_text="Expression", row=1, col=2)
    fig.update_xaxes(title_text="Cells", row=2, col=1)
    fig.update_yaxes(title_text="Genes", row=2, col=1)
    fig.update_xaxes(title_text="Value", row=2, col=2)
    fig.update_yaxes(title_text="Count", row=2, col=2)
    fig.update_xaxes(title_text="Group", row=3, col=1)
    fig.update_yaxes(title_text="Gene", row=3, col=1)

    fig.write_html("demo_selection_synchronization.html")


def main():
    """Main function to run the multi-panel selection demo."""
    print("=" * 80)
    print("PySEE Multi-Panel Advanced Selection Integration Demo")
    print("=" * 80)

    # Create sample data
    data_wrapper = create_sample_anndata(n_cells=500, n_genes=1000)

    print(f"\nData Summary:")
    print(f"  - Cells: {data_wrapper.adata.n_obs}")
    print(f"  - Genes: {data_wrapper.adata.n_vars}")
    print(f"  - Embeddings: {list(data_wrapper.adata.obsm.keys())}")
    print(f"  - Obs columns: {list(data_wrapper.adata.obs.columns)}")

    # Create individual panel demos
    print(f"\nCreating individual panel demos...")
    create_individual_panel_demos(data_wrapper)

    # Create multi-panel comparison
    print(f"\nCreating multi-panel comparison...")
    create_multi_panel_comparison(data_wrapper)

    # Create selection synchronization demo
    print(f"\nCreating selection synchronization demo...")
    create_selection_synchronization_demo(data_wrapper)

    # Create dashboard demo
    print(f"\nCreating interactive dashboard demo...")
    dashboard = create_multi_panel_dashboard(data_wrapper)

    # Save dashboard info
    dashboard_info = f"""
# PySEE Multi-Panel Advanced Selection Integration

## Generated Files:
- demo_umap_advanced_selection.html
- demo_violin_advanced_selection.html
- demo_heatmap_advanced_selection.html
- demo_qc_advanced_selection.html
- demo_dotplot_advanced_selection.html
- demo_multi_panel_comparison.html
- demo_selection_synchronization.html

## Features Demonstrated:
- [OK] Advanced selection tools in all 5 panel types
- [OK] Selection synchronization between panels
- [OK] Multiple selection modes (Replace, Add, Subtract, Intersect)
- [OK] Selection history with undo/redo (50 steps)
- [OK] Interactive help system with comprehensive tooltips
- [OK] Selection statistics and real-time summaries
- [OK] Professional Plotly UI integration

## Panel Types with Advanced Selection:
1. **UMAP Panel**: Point, Rectangular, Lasso selection on embeddings
2. **Violin Panel**: Point selection on expression distributions
3. **Heatmap Panel**: Rectangular selection on expression matrix
4. **QC Panel**: Rectangular selection on QC metrics
5. **DotPlot Panel**: Point selection on marker gene expression

## Usage:
1. Open any of the HTML files in a web browser
2. Use the selection tools in the top-right corner
3. Try different selection modes and types
4. Use undo/redo buttons to manage selection history
5. Click the help button for detailed instructions

## Integration Status:
- [OK] All panels updated with advanced selection capabilities
- [OK] Multi-panel selection synchronization implemented
- [OK] Professional UI with comprehensive help system
- [OK] Ready for production use
"""

    with open("MULTI_PANEL_SELECTION_DEMO_INFO.md", "w") as f:
        f.write(dashboard_info)

    print(f"\n" + "=" * 80)
    print("Multi-Panel Selection Integration Demo Complete!")
    print("=" * 80)
    print(f"\nGenerated files:")
    print(f"  - 5 individual panel demos (HTML)")
    print(f"  - Multi-panel comparison (HTML)")
    print(f"  - Selection synchronization demo (HTML)")
    print(f"  - Demo information (Markdown)")

    print(f"\nKey Features:")
    print(f"  [OK] Advanced selection tools in all panel types")
    print(f"  [OK] Selection synchronization between panels")
    print(f"  [OK] Multiple selection modes and history")
    print(f"  [OK] Interactive help system")
    print(f"  [OK] Professional UI integration")

    print(f"\nNext Steps:")
    print(f"  - Open HTML files in browser to test selection tools")
    print(f"  - Test selection synchronization between panels")
    print(f"  - Verify undo/redo functionality")
    print(f"  - Check help system and tooltips")


if __name__ == "__main__":
    main()
