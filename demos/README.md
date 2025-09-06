# PySEE Demos

This directory contains demonstration scripts and outputs for PySEE features.

## Structure

```
demos/
├── html/                    # Interactive HTML demos
│   ├── demo_umap_advanced_selection.html
│   ├── demo_violin_advanced_selection.html
│   ├── demo_heatmap_advanced_selection.html
│   ├── demo_qc_advanced_selection.html
│   ├── demo_dotplot_advanced_selection.html
│   ├── demo_multi_panel_comparison.html
│   └── demo_selection_synchronization.html
├── selection/               # Advanced Selection Tools demos
│   ├── demo_advanced_selection.py
│   ├── demo_multi_panel_selection.py
│   └── create_png_exports.py
└── README.md               # This file
```

## Quick Start

### Interactive HTML Demos
1. Open any HTML file in `html/` directory in your web browser
2. Use the selection tools in the top-right corner
3. Try different selection modes and types
4. Use undo/redo buttons to manage selection history
5. Click the help button for detailed instructions

### Python Demos
1. Run individual panel demos:
   ```bash
   cd demos/selection
   python demo_advanced_selection.py
   ```

2. Run multi-panel integration demo:
   ```bash
   cd demos/selection
   python demo_multi_panel_selection.py
   ```

3. Generate static PNG exports:
   ```bash
   cd demos/selection
   python create_png_exports.py
   ```

## Features Demonstrated

### Advanced Selection Tools
- **5 Selection Types**: Point, Rectangular, Lasso, Polygon, Circle
- **4 Selection Modes**: Replace, Add, Subtract, Intersect
- **Selection History**: Undo/Redo with 50-step history
- **Interactive Help**: Comprehensive tooltips and guides
- **Selection Statistics**: Real-time counts and percentages
- **Professional UI**: Plotly integration with modern interface

### Multi-Panel Integration
- **Cross-Panel Selection**: Select in UMAP, see in all panels
- **Live Synchronization**: Automatic panel updates
- **Panel Types**: UMAP, Violin, Heatmap, QC, DotPlot
- **Selection Propagation**: Coordinated selection across panels

## Documentation

- [Advanced Selection Implementation Summary](../docs/features/ADVANCED_SELECTION_IMPLEMENTATION_SUMMARY.md)
- [Multi-Panel Selection Demo Info](../docs/features/MULTI_PANEL_SELECTION_DEMO_INFO.md)
- [iSEE Comparison Analysis](../docs/development/ISEE_COMPARISON_ANALYSIS.md)

## Requirements

- Python 3.9+
- PySEE package installed
- Web browser for HTML demos
- Optional: Kaleido for PNG export

## Troubleshooting

### HTML Demos Not Loading
- Ensure you're opening files in a modern web browser
- Check that all files are in the `html/` directory
- Try refreshing the page

### Python Demos Failing
- Ensure PySEE is properly installed
- Check that all dependencies are available
- Run from the correct directory

### PNG Export Issues
- Install Kaleido: `pip install kaleido`
- Or use the Matplotlib fallback in `create_png_exports.py`
