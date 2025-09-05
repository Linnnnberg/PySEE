"""
Performance tests for PySEE panel rendering.

This module tests the rendering performance of different PySEE panels
across various dataset sizes.
"""

import pytest
import time
from typing import Dict, Any
from pysee import PySEE
from pysee.panels.umap import UMAPPanel
from pysee.panels.violin import ViolinPanel
from pysee.panels.heatmap import HeatmapPanel
from pysee.panels.qc import QCPanel

from .fixtures.dataset_fixtures import DatasetFixtures
from .utils.performance_utils import PerformanceBenchmark, PerformanceTargets


class TestRenderingPerformance:
    """Test rendering performance for all PySEE panels."""

    @pytest.fixture(scope="class")
    def datasets(self):
        """Get all test datasets."""
        return DatasetFixtures.get_all_datasets()

    @pytest.fixture(scope="class")
    def dataset_info(self, datasets):
        """Get dataset information."""
        return {name: DatasetFixtures.get_dataset_info(adata) for name, adata in datasets.items()}

    def test_umap_panel_rendering(self, datasets, dataset_info):
        """Test UMAP panel rendering performance."""
        benchmark = PerformanceBenchmark("UMAP Panel Rendering")

        for name, adata in datasets.items():
            if not dataset_info[name]["has_umap"]:
                continue

            app = PySEE(adata)
            umap_panel = UMAPPanel("umap", title=f"UMAP - {name}")
            app.add_panel("umap", umap_panel)

            def render_umap():
                return app.render_panel("umap")

            result = benchmark.run_benchmark(render_umap, iterations=3)
            result["dataset"] = name
            result["dataset_size"] = dataset_info[name]["n_cells"]
            benchmark.add_result(result)

            # Check performance targets
            dataset_size = self._get_dataset_size_category(dataset_info[name]["n_cells"])
            assert PerformanceTargets.check_rendering_performance(
                dataset_size, result["mean_time"]
            ), f"UMAP rendering too slow for {name}: {result['mean_time']:.3f}s"

    def test_violin_panel_rendering(self, datasets, dataset_info):
        """Test Violin panel rendering performance."""
        benchmark = PerformanceBenchmark("Violin Panel Rendering")

        for name, adata in datasets.items():
            app = PySEE(adata)
            violin_panel = ViolinPanel("violin", gene="Gene_0001", title=f"Violin - {name}")
            app.add_panel("violin", violin_panel)

            def render_violin():
                return app.render_panel("violin")

            result = benchmark.run_benchmark(render_violin, iterations=3)
            result["dataset"] = name
            result["dataset_size"] = dataset_info[name]["n_cells"]
            benchmark.add_result(result)

            # Check performance targets
            dataset_size = self._get_dataset_size_category(dataset_info[name]["n_cells"])
            assert PerformanceTargets.check_rendering_performance(
                dataset_size, result["mean_time"]
            ), f"Violin rendering too slow for {name}: {result['mean_time']:.3f}s"

    def test_heatmap_panel_rendering(self, datasets, dataset_info):
        """Test Heatmap panel rendering performance."""
        benchmark = PerformanceBenchmark("Heatmap Panel Rendering")

        for name, adata in datasets.items():
            app = PySEE(adata)
            heatmap_panel = HeatmapPanel("heatmap", title=f"Heatmap - {name}")
            app.add_panel("heatmap", heatmap_panel)

            def render_heatmap():
                return app.render_panel("heatmap")

            result = benchmark.run_benchmark(render_heatmap, iterations=3)
            result["dataset"] = name
            result["dataset_size"] = dataset_info[name]["n_cells"]
            benchmark.add_result(result)

            # Check performance targets
            dataset_size = self._get_dataset_size_category(dataset_info[name]["n_cells"])
            assert PerformanceTargets.check_rendering_performance(
                dataset_size, result["mean_time"]
            ), f"Heatmap rendering too slow for {name}: {result['mean_time']:.3f}s"

    def test_qc_panel_rendering(self, datasets, dataset_info):
        """Test QC panel rendering performance."""
        benchmark = PerformanceBenchmark("QC Panel Rendering")

        for name, adata in datasets.items():
            app = PySEE(adata)
            qc_panel = QCPanel("qc", title=f"QC - {name}")
            app.add_panel("qc", qc_panel)

            def render_qc():
                return app.render_panel("qc")

            result = benchmark.run_benchmark(render_qc, iterations=3)
            result["dataset"] = name
            result["dataset_size"] = dataset_info[name]["n_cells"]
            benchmark.add_result(result)

            # Check performance targets
            dataset_size = self._get_dataset_size_category(dataset_info[name]["n_cells"])
            assert PerformanceTargets.check_rendering_performance(
                dataset_size, result["mean_time"]
            ), f"QC rendering too slow for {name}: {result['mean_time']:.3f}s"

    def test_multi_panel_dashboard_rendering(self, datasets, dataset_info):
        """Test multi-panel dashboard rendering performance."""
        benchmark = PerformanceBenchmark("Multi-Panel Dashboard Rendering")

        for name, adata in datasets.items():
            app = PySEE(adata, title=f"Dashboard - {name}")

            # Add all panels
            if dataset_info[name]["has_umap"]:
                umap_panel = UMAPPanel("umap", title="UMAP")
                app.add_panel("umap", umap_panel)

            violin_panel = ViolinPanel("violin", gene="Gene_0001", title="Violin")
            app.add_panel("violin", violin_panel)

            heatmap_panel = HeatmapPanel("heatmap", title="Heatmap")
            app.add_panel("heatmap", heatmap_panel)

            qc_panel = QCPanel("qc", title="QC")
            app.add_panel("qc", qc_panel)

            def render_all_panels():
                results = {}
                if dataset_info[name]["has_umap"]:
                    results["umap"] = app.render_panel("umap")
                results["violin"] = app.render_panel("violin")
                results["heatmap"] = app.render_panel("heatmap")
                results["qc"] = app.render_panel("qc")
                return results

            result = benchmark.run_benchmark(render_all_panels, iterations=2)
            result["dataset"] = name
            result["dataset_size"] = dataset_info[name]["n_cells"]
            benchmark.add_result(result)

            # Check performance targets (allow more time for multi-panel)
            dataset_size = self._get_dataset_size_category(dataset_info[name]["n_cells"])
            target_time = PerformanceTargets.get_rendering_target(dataset_size) * 2  # 2x for multi-panel
            assert result["mean_time"] <= target_time, f"Multi-panel rendering too slow for {name}: {result['mean_time']:.3f}s"

    def _get_dataset_size_category(self, n_cells: int) -> str:
        """Get dataset size category based on number of cells."""
        if n_cells < 5000:
            return "small"
        elif n_cells < 50000:
            return "medium"
        elif n_cells < 200000:
            return "large"
        else:
            return "very_large"

    @pytest.fixture(scope="class", autouse=True)
    def generate_performance_report(self, datasets, dataset_info):
        """Generate performance report after all tests."""
        yield

        # This will run after all tests in the class
        print("\n" + "=" * 60)
        print("RENDERING PERFORMANCE TEST SUMMARY")
        print("=" * 60)

        for name, adata in datasets.items():
            info = dataset_info[name]
            print(f"\nDataset: {name}")
            print(f"  Cells: {info['n_cells']:,}")
            print(f"  Genes: {info['n_genes']:,}")
            print(f"  Memory: {info['memory_usage_mb']:.1f} MB")
            print(f"  Sparsity: {info['sparsity']:.3f}")


if __name__ == "__main__":
    # Run performance tests
    pytest.main([__file__, "-v", "--tb=short"])
