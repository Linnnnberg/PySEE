"""
Performance testing utilities for PySEE.

This module provides utilities for measuring and analyzing performance
metrics across different PySEE components.
"""

import time
import psutil
import os
import gc
from typing import Dict, Any, Callable, List, Optional
import numpy as np
import pandas as pd
from contextlib import contextmanager
import tracemalloc


class PerformanceProfiler:
    """Profiler for measuring PySEE performance metrics."""
    
    def __init__(self):
        self.process = psutil.Process(os.getpid())
        self.start_memory = None
        self.start_time = None
        self.peak_memory = 0
        self.tracemalloc_started = False
    
    def start_profiling(self):
        """Start profiling session."""
        gc.collect()  # Clean up before starting
        self.start_memory = self.process.memory_info().rss / (1024 * 1024)  # MB
        self.start_time = time.time()
        self.peak_memory = self.start_memory
        
        # Start memory tracing
        tracemalloc.start()
        self.tracemalloc_started = True
    
    def stop_profiling(self) -> Dict[str, Any]:
        """Stop profiling and return metrics."""
        if not self.tracemalloc_started:
            raise RuntimeError("Profiling not started")
        
        end_time = time.time()
        end_memory = self.process.memory_info().rss / (1024 * 1024)  # MB
        
        # Get memory trace
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        self.tracemalloc_started = False
        
        return {
            'execution_time': end_time - self.start_time,
            'start_memory_mb': self.start_memory,
            'end_memory_mb': end_memory,
            'memory_delta_mb': end_memory - self.start_memory,
            'peak_memory_mb': peak / (1024 * 1024),
            'traced_memory_mb': current / (1024 * 1024),
        }
    
    def measure_function(self, func: Callable, *args, **kwargs) -> Dict[str, Any]:
        """Measure performance of a function call."""
        self.start_profiling()
        try:
            result = func(*args, **kwargs)
            metrics = self.stop_profiling()
            metrics['result'] = result
            return metrics
        except Exception as e:
            self.stop_profiling()
            raise e


@contextmanager
def measure_performance(description: str = ""):
    """Context manager for measuring performance."""
    profiler = PerformanceProfiler()
    profiler.start_profiling()
    try:
        yield profiler
    finally:
        metrics = profiler.stop_profiling()
        if description:
            print(f"{description}: {metrics['execution_time']:.3f}s, "
                  f"Memory: {metrics['memory_delta_mb']:.1f}MB")


class PerformanceBenchmark:
    """Benchmark runner for PySEE performance tests."""
    
    def __init__(self, name: str):
        self.name = name
        self.results = []
    
    def run_benchmark(self, func: Callable, *args, iterations: int = 3, **kwargs) -> Dict[str, Any]:
        """Run a benchmark multiple times and return statistics."""
        times = []
        memory_deltas = []
        results = []
        
        for i in range(iterations):
            profiler = PerformanceProfiler()
            profiler.start_profiling()
            
            try:
                result = func(*args, **kwargs)
                metrics = profiler.stop_profiling()
                
                times.append(metrics['execution_time'])
                memory_deltas.append(metrics['memory_delta_mb'])
                results.append(result)
                
                # Clean up between iterations
                gc.collect()
                
            except Exception as e:
                profiler.stop_profiling()
                raise e
        
        return {
            'benchmark_name': self.name,
            'iterations': iterations,
            'mean_time': np.mean(times),
            'std_time': np.std(times),
            'min_time': np.min(times),
            'max_time': np.max(times),
            'mean_memory_delta': np.mean(memory_deltas),
            'std_memory_delta': np.std(memory_deltas),
            'min_memory_delta': np.min(memory_deltas),
            'max_memory_delta': np.max(memory_deltas),
            'results': results,
        }
    
    def add_result(self, result: Dict[str, Any]):
        """Add a benchmark result."""
        self.results.append(result)
    
    def get_summary(self) -> Dict[str, Any]:
        """Get summary of all benchmark results."""
        if not self.results:
            return {}
        
        return {
            'benchmark_name': self.name,
            'total_benchmarks': len(self.results),
            'results': self.results,
        }


class PerformanceTargets:
    """Performance targets and thresholds for PySEE."""
    
    # Memory usage targets (MB)
    MEMORY_TARGETS = {
        'small': 500,      # < 500 MB for small datasets
        'medium': 2000,    # < 2 GB for medium datasets
        'large': 8000,     # < 8 GB for large datasets
        'very_large': 32000,  # < 32 GB for very large datasets
    }
    
    # Rendering time targets (seconds)
    RENDERING_TARGETS = {
        'small': 2.0,      # < 2 seconds for small datasets
        'medium': 5.0,     # < 5 seconds for medium datasets
        'large': 10.0,     # < 10 seconds for large datasets
        'very_large': 30.0,  # < 30 seconds for very large datasets
    }
    
    # Interactive response time targets (seconds)
    INTERACTION_TARGETS = {
        'zoom_pan': 0.1,   # < 100ms for zoom/pan
        'selection': 0.5,  # < 500ms for selection
        'propagation': 0.2,  # < 200ms for selection propagation
        'code_export': 1.0,  # < 1 second for code export
    }
    
    @classmethod
    def get_memory_target(cls, dataset_size: str) -> float:
        """Get memory target for dataset size."""
        return cls.MEMORY_TARGETS.get(dataset_size, float('inf'))
    
    @classmethod
    def get_rendering_target(cls, dataset_size: str) -> float:
        """Get rendering time target for dataset size."""
        return cls.RENDERING_TARGETS.get(dataset_size, float('inf'))
    
    @classmethod
    def check_memory_performance(cls, dataset_size: str, actual_memory: float) -> bool:
        """Check if memory usage meets target."""
        target = cls.get_memory_target(dataset_size)
        return actual_memory <= target
    
    @classmethod
    def check_rendering_performance(cls, dataset_size: str, actual_time: float) -> bool:
        """Check if rendering time meets target."""
        target = cls.get_rendering_target(dataset_size)
        return actual_time <= target


class PerformanceReporter:
    """Generate performance reports and visualizations."""
    
    @staticmethod
    def generate_summary_report(results: List[Dict[str, Any]]) -> str:
        """Generate a text summary report."""
        report = []
        report.append("PySEE Performance Test Summary")
        report.append("=" * 50)
        
        for result in results:
            report.append(f"\nBenchmark: {result.get('benchmark_name', 'Unknown')}")
            report.append(f"  Mean Time: {result.get('mean_time', 0):.3f}s")
            report.append(f"  Memory Delta: {result.get('mean_memory_delta', 0):.1f}MB")
            report.append(f"  Iterations: {result.get('iterations', 0)}")
        
        return "\n".join(report)
    
    @staticmethod
    def save_results_to_json(results: List[Dict[str, Any]], filename: str):
        """Save results to JSON file."""
        import json
        
        # Convert numpy types to Python types for JSON serialization
        def convert_numpy(obj):
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            return obj
        
        # Clean results for JSON serialization
        clean_results = []
        for result in results:
            clean_result = {}
            for key, value in result.items():
                if key != 'results':  # Skip the actual results objects
                    clean_result[key] = convert_numpy(value)
            clean_results.append(clean_result)
        
        with open(filename, 'w') as f:
            json.dump(clean_results, f, indent=2)
    
    @staticmethod
    def create_performance_plot(results: List[Dict[str, Any]], output_file: str):
        """Create performance visualization plot."""
        import matplotlib.pyplot as plt
        
        # Extract data for plotting
        benchmark_names = [r.get('benchmark_name', 'Unknown') for r in results]
        mean_times = [r.get('mean_time', 0) for r in results]
        mean_memory = [r.get('mean_memory_delta', 0) for r in results]
        
        # Create subplots
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
        
        # Plot execution times
        ax1.bar(benchmark_names, mean_times)
        ax1.set_ylabel('Execution Time (seconds)')
        ax1.set_title('PySEE Performance Benchmarks - Execution Time')
        ax1.tick_params(axis='x', rotation=45)
        
        # Plot memory usage
        ax2.bar(benchmark_names, mean_memory)
        ax2.set_ylabel('Memory Delta (MB)')
        ax2.set_xlabel('Benchmark')
        ax2.set_title('PySEE Performance Benchmarks - Memory Usage')
        ax2.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
