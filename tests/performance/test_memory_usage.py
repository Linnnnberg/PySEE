"""
Memory usage tests for PySEE.

This module tests memory consumption patterns across different PySEE
components and dataset sizes.
"""

import pytest
import gc
from typing import Dict, Any
from pysee import PySEE
from pysee.panels.umap import UMAPPanel
from pysee.panels.violin import ViolinPanel
from pysee.panels.heatmap import HeatmapPanel
from pysee.panels.qc import QCPanel

from .fixtures.dataset_fixtures import DatasetFixtures
from .utils.performance_utils import PerformanceProfiler, PerformanceTargets


class TestMemoryUsage:
    """Test memory usage patterns for PySEE components."""
    
    @pytest.fixture(scope="class")
    def datasets(self):
        """Get all test datasets."""
        return DatasetFixtures.get_all_datasets()
    
    @pytest.fixture(scope="class")
    def dataset_info(self, datasets):
        """Get dataset information."""
        return {name: DatasetFixtures.get_dataset_info(adata) 
                for name, adata in datasets.items()}
    
    def test_data_loading_memory(self, datasets, dataset_info):
        """Test memory usage during data loading."""
        for name, adata in datasets.items():
            profiler = PerformanceProfiler()
            
            def load_data():
                # Simulate data loading
                app = PySEE(adata)
                return app
            
            metrics = profiler.measure_function(load_data)
            
            # Check memory targets
            dataset_size = self._get_dataset_size_category(dataset_info[name]['n_cells'])
            target_memory = PerformanceTargets.get_memory_target(dataset_size)
            
            assert metrics['memory_delta_mb'] <= target_memory, \
                f"Data loading memory too high for {name}: {metrics['memory_delta_mb']:.1f}MB > {target_memory}MB"
            
            print(f"Data loading memory for {name}: {metrics['memory_delta_mb']:.1f}MB")
    
    def test_panel_creation_memory(self, datasets, dataset_info):
        """Test memory usage during panel creation."""
        for name, adata in datasets.items():
            app = PySEE(adata)
            
            # Test each panel type
            panels = [
                ('umap', UMAPPanel("umap", title="UMAP")),
                ('violin', ViolinPanel("violin", gene="Gene_0001", title="Violin")),
                ('heatmap', HeatmapPanel("heatmap", title="Heatmap")),
                ('qc', QCPanel("qc", title="QC")),
            ]
            
            for panel_name, panel in panels:
                profiler = PerformanceProfiler()
                
                def create_panel():
                    app.add_panel(panel_name, panel)
                    return panel
                
                metrics = profiler.measure_function(create_panel)
                
                # Check memory targets (allow some overhead for panel creation)
                dataset_size = self._get_dataset_size_category(dataset_info[name]['n_cells'])
                target_memory = PerformanceTargets.get_memory_target(dataset_size) * 0.1  # 10% of total
                
                assert metrics['memory_delta_mb'] <= target_memory, \
                    f"Panel creation memory too high for {name}/{panel_name}: {metrics['memory_delta_mb']:.1f}MB"
                
                print(f"Panel creation memory for {name}/{panel_name}: {metrics['memory_delta_mb']:.1f}MB")
    
    def test_rendering_memory(self, datasets, dataset_info):
        """Test memory usage during panel rendering."""
        for name, adata in datasets.items():
            app = PySEE(adata)
            
            # Add all panels
            if dataset_info[name]['has_umap']:
                umap_panel = UMAPPanel("umap", title="UMAP")
                app.add_panel("umap", umap_panel)
            
            violin_panel = ViolinPanel("violin", gene="Gene_0001", title="Violin")
            app.add_panel("violin", violin_panel)
            
            heatmap_panel = HeatmapPanel("heatmap", title="Heatmap")
            app.add_panel("heatmap", heatmap_panel)
            
            qc_panel = QCPanel("qc", title="QC")
            app.add_panel("qc", qc_panel)
            
            # Test rendering memory for each panel
            panels_to_test = ['violin', 'heatmap', 'qc']
            if dataset_info[name]['has_umap']:
                panels_to_test.append('umap')
            
            for panel_name in panels_to_test:
                profiler = PerformanceProfiler()
                
                def render_panel():
                    return app.render_panel(panel_name)
                
                metrics = profiler.measure_function(render_panel)
                
                # Check memory targets (allow some overhead for rendering)
                dataset_size = self._get_dataset_size_category(dataset_info[name]['n_cells'])
                target_memory = PerformanceTargets.get_memory_target(dataset_size) * 0.2  # 20% of total
                
                assert metrics['memory_delta_mb'] <= target_memory, \
                    f"Rendering memory too high for {name}/{panel_name}: {metrics['memory_delta_mb']:.1f}MB"
                
                print(f"Rendering memory for {name}/{panel_name}: {metrics['memory_delta_mb']:.1f}MB")
    
    def test_memory_scaling(self, datasets, dataset_info):
        """Test memory scaling with dataset size."""
        memory_usage = {}
        
        for name, adata in datasets.items():
            profiler = PerformanceProfiler()
            
            def full_workflow():
                app = PySEE(adata)
                
                # Add all panels
                if dataset_info[name]['has_umap']:
                    umap_panel = UMAPPanel("umap", title="UMAP")
                    app.add_panel("umap", umap_panel)
                
                violin_panel = ViolinPanel("violin", gene="Gene_0001", title="Violin")
                app.add_panel("violin", violin_panel)
                
                heatmap_panel = HeatmapPanel("heatmap", title="Heatmap")
                app.add_panel("heatmap", heatmap_panel)
                
                qc_panel = QCPanel("qc", title="QC")
                app.add_panel("qc", qc_panel)
                
                # Render all panels
                results = {}
                if dataset_info[name]['has_umap']:
                    results['umap'] = app.render_panel("umap")
                results['violin'] = app.render_panel("violin")
                results['heatmap'] = app.render_panel("heatmap")
                results['qc'] = app.render_panel("qc")
                
                return results
            
            metrics = profiler.measure_function(full_workflow)
            memory_usage[name] = {
                'n_cells': dataset_info[name]['n_cells'],
                'memory_mb': metrics['memory_delta_mb'],
                'peak_memory_mb': metrics['peak_memory_mb'],
            }
        
        # Check memory scaling
        print("\nMemory Scaling Analysis:")
        print("Dataset\t\tCells\t\tMemory (MB)\tPeak (MB)")
        print("-" * 60)
        
        for name, usage in memory_usage.items():
            print(f"{name:15}\t{usage['n_cells']:8,}\t{usage['memory_mb']:8.1f}\t{usage['peak_memory_mb']:8.1f}")
        
        # Verify reasonable scaling
        small_memory = memory_usage.get('synthetic_small', {}).get('memory_mb', 0)
        large_memory = memory_usage.get('synthetic_large', {}).get('memory_mb', 0)
        
        if small_memory > 0 and large_memory > 0:
            scaling_factor = large_memory / small_memory
            cell_scaling_factor = memory_usage['synthetic_large']['n_cells'] / memory_usage['synthetic_small']['n_cells']
            
            # Memory should scale sub-linearly with cell count
            assert scaling_factor < cell_scaling_factor * 2, \
                f"Memory scaling too high: {scaling_factor:.2f}x for {cell_scaling_factor:.2f}x cells"
    
    def test_memory_cleanup(self, datasets, dataset_info):
        """Test memory cleanup after operations."""
        for name, adata in datasets.items():
            # Measure memory before
            profiler = PerformanceProfiler()
            profiler.start_profiling()
            initial_memory = profiler.stop_profiling()['end_memory_mb']
            
            # Create and destroy multiple times
            for i in range(3):
                app = PySEE(adata)
                
                # Add panels
                violin_panel = ViolinPanel("violin", gene="Gene_0001", title="Violin")
                app.add_panel("violin", violin_panel)
                
                # Render
                app.render_panel("violin")
                
                # Clean up
                del app
                del violin_panel
                gc.collect()
            
            # Measure memory after
            profiler.start_profiling()
            final_memory = profiler.stop_profiling()['end_memory_mb']
            
            memory_leak = final_memory - initial_memory
            
            # Check for memory leaks (allow small amount for Python overhead)
            assert memory_leak < 100, \
                f"Potential memory leak for {name}: {memory_leak:.1f}MB"
            
            print(f"Memory cleanup test for {name}: {memory_leak:.1f}MB leak")
    
    def _get_dataset_size_category(self, n_cells: int) -> str:
        """Get dataset size category based on number of cells."""
        if n_cells < 5000:
            return 'small'
        elif n_cells < 50000:
            return 'medium'
        elif n_cells < 200000:
            return 'large'
        else:
            return 'very_large'


if __name__ == "__main__":
    # Run memory tests
    pytest.main([__file__, "-v", "--tb=short"])
