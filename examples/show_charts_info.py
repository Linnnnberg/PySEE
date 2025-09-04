"""
Script to display information about the generated PySEE charts.
"""

import os
import glob
from datetime import datetime


def show_charts_info():
    """Display information about generated chart files."""
    print("🎨 PySEE Generated Charts")
    print("=" * 50)
    
    # Find all HTML files
    html_files = glob.glob("*.html")
    
    if not html_files:
        print("❌ No HTML chart files found.")
        return
    
    print(f"📊 Found {len(html_files)} chart files:")
    print()
    
    for i, file in enumerate(sorted(html_files), 1):
        file_size = os.path.getsize(file)
        file_time = datetime.fromtimestamp(os.path.getmtime(file))
        
        print(f"{i:2d}. {file}")
        print(f"    Size: {file_size:,} bytes")
        print(f"    Created: {file_time.strftime('%Y-%m-%d %H:%M:%S')}")
        print()
    
    print("🌐 To view the charts:")
    print("   1. Open any of the HTML files in your web browser")
    print("   2. The charts are interactive - you can zoom, pan, and hover for details")
    print("   3. Each chart shows different aspects of the PySEE visualization system")
    print()
    
    print("📋 Chart Descriptions:")
    print("   • qc_panel_demo.html - Basic QC metrics visualization")
    print("   • qc_panel_custom_thresholds.html - QC panel with custom filtering thresholds")
    print("   • umap_panel_demo.html - UMAP dimensionality reduction visualization")
    print("   • violin_panel_demo.html - Gene expression violin plots by cell type")
    print("   • heatmap_panel_demo.html - Gene expression heatmap with clustering")
    print("   • qc_panel_multi_demo.html - QC panel in multi-panel context")
    print("   • qc_filtering_workflow.html - QC filtering workflow demonstration")
    print()
    
    print("🚀 PySEE Features Demonstrated:")
    print("   ✅ Interactive Plotly visualizations")
    print("   ✅ QC metrics with configurable thresholds")
    print("   ✅ Multi-panel dashboard integration")
    print("   ✅ Code export for reproducible analysis")
    print("   ✅ Realistic bioinformatics data simulation")
    print("   ✅ Professional-quality visualizations")


if __name__ == "__main__":
    show_charts_info()
