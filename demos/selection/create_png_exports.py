#!/usr/bin/env python3
"""
Create PNG exports of the Advanced Selection Tools demo using matplotlib.

This script creates static PNG versions of the selection tool demonstrations
as a fallback when kaleido PNG export fails.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
import os
from typing import Tuple, List

# Import our selection modules
from pysee.panels.selection_types import SelectionType, SelectionMode, SelectionStatistics
from pysee.panels.advanced_selection import AdvancedSelectionManager
from demo_advanced_selection import create_sample_data


def create_selection_demo_png(data_wrapper, output_dir: str = "selection_demo_output") -> None:
    """Create PNG versions of selection demonstrations using matplotlib."""
    print("🖼️  Creating PNG exports using matplotlib...")
    
    # Get data
    umap_coords = data_wrapper.get_embedding_data('X_umap')
    cell_types = data_wrapper.adata.obs['cell_type'].values
    
    # Create selection manager
    selection_manager = AdvancedSelectionManager(data_wrapper)
    
    # Set up the figure with subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('PySEE Advanced Selection Tools Demo', fontsize=16, fontweight='bold')
    
    # Color mapping for cell types
    type_colors = {'Type_A': 'red', 'Type_B': 'blue', 'Type_C': 'green'}
    colors = [type_colors[ct] for ct in cell_types]
    
    # 1. Rectangular Selection
    ax = axes[0, 0]
    ax.scatter(umap_coords[:, 0], umap_coords[:, 1], c=colors, alpha=0.6, s=20)
    
    # Add rectangle
    rect_bounds = (-1, -1, 1, 1)
    rect = patches.Rectangle((rect_bounds[0], rect_bounds[1]), 
                           rect_bounds[2] - rect_bounds[0], 
                           rect_bounds[3] - rect_bounds[1],
                           linewidth=2, edgecolor='black', facecolor='none', 
                           linestyle='--', alpha=0.8)
    ax.add_patch(rect)
    
    # Calculate selection statistics
    rect_selection = selection_manager.create_rectangular_selection(rect_bounds, umap_coords)
    n_selected = np.sum(rect_selection)
    percentage = (n_selected / len(rect_selection)) * 100
    
    ax.set_title(f'Rectangular Selection\n{n_selected}/{len(rect_selection)} cells ({percentage:.1f}%)')
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.grid(True, alpha=0.3)
    
    # 2. Lasso Selection
    ax = axes[0, 1]
    ax.scatter(umap_coords[:, 0], umap_coords[:, 1], c=colors, alpha=0.6, s=20)
    
    # Add lasso (circular)
    theta = np.linspace(0, 2*np.pi, 20)
    lasso_path = [(1.5 * np.cos(t), 1.5 * np.sin(t)) for t in theta]
    lasso_x, lasso_y = zip(*lasso_path)
    ax.plot(lasso_x, lasso_y, 'k--', linewidth=2, alpha=0.8)
    
    # Calculate selection
    lasso_selection = selection_manager.create_lasso_selection(lasso_path, umap_coords)
    n_selected = np.sum(lasso_selection)
    percentage = (n_selected / len(lasso_selection)) * 100
    
    ax.set_title(f'Lasso Selection\n{n_selected}/{len(lasso_selection)} cells ({percentage:.1f}%)')
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.grid(True, alpha=0.3)
    
    # 3. Polygon Selection
    ax = axes[0, 2]
    ax.scatter(umap_coords[:, 0], umap_coords[:, 1], c=colors, alpha=0.6, s=20)
    
    # Add polygon (triangle)
    polygon_vertices = [(-2, -1), (0, 2), (2, -1)]
    poly = patches.Polygon(polygon_vertices, closed=True, fill=False, 
                          edgecolor='black', linewidth=2, linestyle='--', alpha=0.8)
    ax.add_patch(poly)
    
    # Calculate selection
    poly_selection = selection_manager.create_polygon_selection(polygon_vertices, umap_coords)
    n_selected = np.sum(poly_selection)
    percentage = (n_selected / len(poly_selection)) * 100
    
    ax.set_title(f'Polygon Selection\n{n_selected}/{len(poly_selection)} cells ({percentage:.1f}%)')
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.grid(True, alpha=0.3)
    
    # 4. Selection Modes Comparison
    ax = axes[1, 0]
    
    # Base selection (rectangle)
    base_selection = selection_manager.create_rectangular_selection((-1, -1, 1, 1), umap_coords)
    
    # New selection (circle)
    new_selection = selection_manager.create_circle_selection((0.5, 0.5), 1.2, umap_coords)
    
    # Combined selection (ADD mode)
    combined_selection = base_selection | new_selection
    
    # Color points based on selection
    point_colors = ['red' if s else 'lightgray' for s in combined_selection]
    ax.scatter(umap_coords[:, 0], umap_coords[:, 1], c=point_colors, alpha=0.7, s=20)
    
    # Add shapes
    rect = patches.Rectangle((-1, -1), 2, 2, linewidth=1, edgecolor='blue', 
                           facecolor='none', linestyle=':')
    ax.add_patch(rect)
    
    circle = patches.Circle((0.5, 0.5), 1.2, linewidth=1, edgecolor='green', 
                          facecolor='none', linestyle=':')
    ax.add_patch(circle)
    
    n_selected = np.sum(combined_selection)
    percentage = (n_selected / len(combined_selection)) * 100
    ax.set_title(f'Selection Modes (ADD)\n{n_selected}/{len(combined_selection)} cells ({percentage:.1f}%)')
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.grid(True, alpha=0.3)
    
    # 5. Selection Statistics
    ax = axes[1, 1]
    
    # Create different selections for statistics
    selections = [
        ("Small Rect", selection_manager.create_rectangular_selection((-2, -2, -1, -1), umap_coords)),
        ("Large Rect", selection_manager.create_rectangular_selection((-3, -3, 3, 3), umap_coords)),
        ("Circle", selection_manager.create_circle_selection((0, 0), 1.5, umap_coords)),
        ("Polygon", selection_manager.create_polygon_selection([(-1, -1), (1, -1), (0, 1)], umap_coords))
    ]
    
    names = [name for name, _ in selections]
    counts = [np.sum(sel) for _, sel in selections]
    percentages = [(count / len(sel)) * 100 for count, (_, sel) in zip(counts, selections)]
    
    # Bar chart
    bars = ax.bar(names, counts, color=['skyblue', 'lightcoral', 'lightgreen', 'gold'], alpha=0.8)
    ax.set_title('Selection Statistics')
    ax.set_ylabel('Number of Selected Cells')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add percentage labels on bars
    for bar, pct in zip(bars, percentages):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 20,
                f'{pct:.1f}%', ha='center', va='bottom', fontweight='bold')
    
    # 6. Cell Type Distribution
    ax = axes[1, 2]
    
    # Count cell types
    type_counts = {}
    for ct in cell_types:
        type_counts[ct] = type_counts.get(ct, 0) + 1
    
    # Pie chart
    wedges, texts, autotexts = ax.pie(type_counts.values(), labels=type_counts.keys(), 
                                     autopct='%1.1f%%', startangle=90,
                                     colors=['red', 'blue', 'green'])
    ax.set_title('Cell Type Distribution')
    
    # Improve text visibility
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
    
    plt.tight_layout()
    
    # Save PNG
    png_path = os.path.join(output_dir, "advanced_selection_tools_demo.png")
    plt.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"✅ PNG saved: {png_path}")
    
    # Create individual selection type PNGs
    create_individual_pngs(data_wrapper, selection_manager, output_dir)


def create_individual_pngs(data_wrapper, selection_manager, output_dir: str) -> None:
    """Create individual PNG files for each selection type."""
    print("   Creating individual selection type PNGs...")
    
    umap_coords = data_wrapper.get_embedding_data('X_umap')
    cell_types = data_wrapper.adata.obs['cell_type'].values
    
    # Color mapping
    type_colors = {'Type_A': 'red', 'Type_B': 'blue', 'Type_C': 'green'}
    colors = [type_colors[ct] for ct in cell_types]
    
    # Individual selection demonstrations
    selections_to_demo = [
        {
            'name': 'rectangular_selection_static',
            'title': 'Rectangular Selection Tool',
            'bounds': (-1, -1, 1, 1),
            'type': 'rectangle'
        },
        {
            'name': 'lasso_selection_static', 
            'title': 'Lasso Selection Tool',
            'path': [(1.5 * np.cos(t), 1.5 * np.sin(t)) for t in np.linspace(0, 2*np.pi, 20)],
            'type': 'lasso'
        },
        {
            'name': 'polygon_selection_static',
            'title': 'Polygon Selection Tool', 
            'vertices': [(-2, -1), (0, 2), (2, -1)],
            'type': 'polygon'
        }
    ]
    
    for demo in selections_to_demo:
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        
        # Base scatter plot
        ax.scatter(umap_coords[:, 0], umap_coords[:, 1], c=colors, alpha=0.6, s=30)
        
        # Add selection area and calculate selection
        if demo['type'] == 'rectangle':
            bounds = demo['bounds']
            rect = patches.Rectangle((bounds[0], bounds[1]), 
                                   bounds[2] - bounds[0], bounds[3] - bounds[1],
                                   linewidth=3, edgecolor='black', facecolor='red',
                                   alpha=0.2, linestyle='--')
            ax.add_patch(rect)
            selection = selection_manager.create_rectangular_selection(bounds, umap_coords)
            
        elif demo['type'] == 'lasso':
            path = demo['path']
            lasso_x, lasso_y = zip(*path)
            ax.plot(lasso_x, lasso_y, 'k--', linewidth=3, alpha=0.8)
            ax.fill(lasso_x, lasso_y, color='red', alpha=0.2)
            selection = selection_manager.create_lasso_selection(path, umap_coords)
            
        elif demo['type'] == 'polygon':
            vertices = demo['vertices']
            poly = patches.Polygon(vertices, closed=True, fill=True, 
                                 facecolor='red', alpha=0.2,
                                 edgecolor='black', linewidth=3, linestyle='--')
            ax.add_patch(poly)
            selection = selection_manager.create_polygon_selection(vertices, umap_coords)
        
        # Highlight selected points
        selected_coords = umap_coords[selection]
        if len(selected_coords) > 0:
            ax.scatter(selected_coords[:, 0], selected_coords[:, 1], 
                      c='yellow', s=50, alpha=0.8, edgecolors='black', linewidth=1)
        
        # Calculate statistics
        n_selected = np.sum(selection)
        percentage = (n_selected / len(selection)) * 100
        
        ax.set_title(f'{demo["title"]}\nSelected: {n_selected:,} / {len(selection):,} cells ({percentage:.1f}%)',
                    fontsize=14, fontweight='bold')
        ax.set_xlabel('UMAP 1', fontsize=12)
        ax.set_ylabel('UMAP 2', fontsize=12)
        ax.grid(True, alpha=0.3)
        
        # Add legend
        legend_elements = [
            plt.scatter([], [], c='red', s=50, alpha=0.8, label='Type A'),
            plt.scatter([], [], c='blue', s=50, alpha=0.8, label='Type B'),
            plt.scatter([], [], c='green', s=50, alpha=0.8, label='Type C'),
            plt.scatter([], [], c='yellow', s=50, alpha=0.8, edgecolors='black', 
                       linewidth=1, label='Selected')
        ]
        ax.legend(handles=legend_elements, loc='upper right')
        
        plt.tight_layout()
        
        # Save
        png_path = os.path.join(output_dir, f"{demo['name']}.png")
        plt.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print(f"     ✅ {demo['name']}.png")


def main():
    """Main function to create PNG exports."""
    print("🖼️  Creating PNG exports for Advanced Selection Tools...")
    
    # Create sample data
    data_wrapper = create_sample_data(n_cells=2000, n_genes=100)
    
    # Create output directory
    output_dir = "selection_demo_output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Create PNG demonstrations
    create_selection_demo_png(data_wrapper, output_dir)
    
    print("\n📊 PNG Export Summary:")
    print("   ✅ Combined demonstration PNG created")
    print("   ✅ Individual selection tool PNGs created")
    print("   ✅ Static visualizations ready for review")
    print(f"   📁 Files saved in: {output_dir}/")
    
    return True


if __name__ == "__main__":
    main()
