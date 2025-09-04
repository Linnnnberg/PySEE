"""
Performance test runner for PySEE.

This script runs comprehensive performance tests and generates reports.
"""

import sys
import os
import json
from datetime import datetime
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from tests.performance.fixtures.dataset_fixtures import DatasetFixtures
from tests.performance.fixtures.dataset_registry import DatasetRegistry
from tests.performance.utils.performance_utils import (
    PerformanceBenchmark, PerformanceReporter, PerformanceTargets
)
from pysee import PySEE
from pysee.panels.umap import UMAPPanel
from pysee.panels.violin import ViolinPanel
from pysee.panels.heatmap import HeatmapPanel
from pysee.panels.qc import QCPanel


def run_quick_performance_test():
    """Run a quick performance test with small datasets."""
    print("🚀 PySEE Quick Performance Test")
    print("=" * 50)
    
    # Get small datasets for quick testing (using fixtures for now)
    datasets = {
        'pbmc3k': DatasetFixtures.get_pbmc3k(),
        'synthetic_small': DatasetFixtures.generate_synthetic_small(),
    }
    print(f"📊 Loaded {len(datasets)} small datasets from fixtures")
    
    results = []
    
    for name, adata in datasets.items():
        print(f"\n📊 Testing dataset: {name}")
        print(f"   Cells: {adata.n_obs:,}, Genes: {adata.n_vars:,}")
        
        # Test each panel type
        panels = [
            ('violin', ViolinPanel("violin", gene=adata.var_names[0], title="Violin")),
            ('heatmap', HeatmapPanel("heatmap", title="Heatmap")),
            ('qc', QCPanel("qc", title="QC")),
        ]
        
        # Add UMAP if available
        if 'X_umap' in adata.obsm:
            panels.append(('umap', UMAPPanel("umap", title="UMAP")))
        
        for panel_name, panel in panels:
            print(f"   Testing {panel_name} panel...")
            
            benchmark = PerformanceBenchmark(f"{name}_{panel_name}")
            
            def test_panel():
                app = PySEE(adata)
                app.add_panel(panel_name, panel)
                return app.render_panel(panel_name)
            
            result = benchmark.run_benchmark(test_panel, iterations=2)
            result['dataset'] = name
            result['panel'] = panel_name
            result['n_cells'] = adata.n_obs
            result['n_genes'] = adata.n_vars
            
            results.append(result)
            
            print(f"     Time: {result['mean_time']:.3f}s ± {result['std_time']:.3f}s")
            print(f"     Memory: {result['mean_memory_delta']:.1f}MB ± {result['std_memory_delta']:.1f}MB")
    
    # Generate report
    print("\n" + "=" * 50)
    print("📋 PERFORMANCE TEST SUMMARY")
    print("=" * 50)
    
    for result in results:
        dataset_size = get_dataset_size_category(result['n_cells'])
        time_target = PerformanceTargets.get_rendering_target(dataset_size)
        memory_target = PerformanceTargets.get_memory_target(dataset_size)
        
        time_ok = result['mean_time'] <= time_target
        memory_ok = result['mean_memory_delta'] <= memory_target
        
        status = "✅" if time_ok and memory_ok else "❌"
        
        print(f"{status} {result['dataset']}/{result['panel']}: "
              f"{result['mean_time']:.3f}s (target: {time_target}s), "
              f"{result['mean_memory_delta']:.1f}MB (target: {memory_target}MB)")
    
    # Save results
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_file = f"tests/performance/results/quick_test_{timestamp}.json"
    
    PerformanceReporter.save_results_to_json(results, results_file)
    print(f"\n💾 Results saved to: {results_file}")
    
    return results


def get_dataset_size_category(n_cells: int) -> str:
    """Get dataset size category based on number of cells."""
    if n_cells < 5000:
        return 'small'
    elif n_cells < 50000:
        return 'medium'
    elif n_cells < 200000:
        return 'large'
    else:
        return 'very_large'


def run_memory_stress_test():
    """Run memory stress test with larger datasets."""
    print("\n🧠 Memory Stress Test")
    print("=" * 50)
    
    # Test with progressively larger datasets
    datasets = {
        'synthetic_small': DatasetFixtures.generate_synthetic_small(),
        'synthetic_medium': DatasetFixtures.generate_synthetic_medium(),
    }
    
    memory_results = []
    
    for name, adata in datasets.items():
        print(f"\n📊 Memory test: {name} ({adata.n_obs:,} cells)")
        
        from tests.performance.utils.performance_utils import PerformanceProfiler
        
        profiler = PerformanceProfiler()
        
        def full_workflow():
            app = PySEE(adata)
            
            # Add all panels
            violin_panel = ViolinPanel("violin", gene="Gene_0001", title="Violin")
            app.add_panel("violin", violin_panel)
            
            heatmap_panel = HeatmapPanel("heatmap", title="Heatmap")
            app.add_panel("heatmap", heatmap_panel)
            
            qc_panel = QCPanel("qc", title="QC")
            app.add_panel("qc", qc_panel)
            
            # Render all panels
            violin_fig = app.render_panel("violin")
            heatmap_fig = app.render_panel("heatmap")
            qc_fig = app.render_panel("qc")
            
            return {
                'violin': violin_fig,
                'heatmap': heatmap_fig,
                'qc': qc_fig
            }
        
        metrics = profiler.measure_function(full_workflow)
        
        memory_results.append({
            'dataset': name,
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars,
            'memory_delta_mb': metrics['memory_delta_mb'],
            'peak_memory_mb': metrics['peak_memory_mb'],
            'execution_time': metrics['execution_time'],
        })
        
        print(f"   Memory delta: {metrics['memory_delta_mb']:.1f}MB")
        print(f"   Peak memory: {metrics['peak_memory_mb']:.1f}MB")
        print(f"   Execution time: {metrics['execution_time']:.3f}s")
    
    return memory_results


def main():
    """Run all performance tests."""
    print("🔬 PySEE Performance Test Suite")
    print("=" * 60)
    
    try:
        # Run quick performance test
        quick_results = run_quick_performance_test()
        
        # Run memory stress test
        memory_results = run_memory_stress_test()
        
        # Generate final summary
        print("\n" + "=" * 60)
        print("🎯 FINAL SUMMARY")
        print("=" * 60)
        
        print(f"✅ Quick performance test: {len(quick_results)} tests completed")
        print(f"✅ Memory stress test: {len(memory_results)} tests completed")
        
        # Check if all tests passed
        all_passed = True
        for result in quick_results:
            dataset_size = get_dataset_size_category(result['n_cells'])
            time_target = PerformanceTargets.get_rendering_target(dataset_size)
            memory_target = PerformanceTargets.get_memory_target(dataset_size)
            
            if result['mean_time'] > time_target or result['mean_memory_delta'] > memory_target:
                all_passed = False
                break
        
        if all_passed:
            print("🎉 All performance tests PASSED!")
        else:
            print("⚠️ Some performance tests FAILED - check results above")
        
        print("\n💡 Next steps:")
        print("   - Run full test suite: pytest tests/performance/ -v")
        print("   - Check results in tests/performance/results/")
        print("   - Optimize any failing components")
        
    except Exception as e:
        print(f"❌ Performance test failed: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
