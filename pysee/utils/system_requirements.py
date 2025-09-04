"""
System requirements checker for PySEE.

This module provides utilities to check system resources and warn users
about potential memory limitations when working with large datasets.
"""

import psutil
import os
from typing import Dict, Any, Optional, Tuple
import warnings


class SystemRequirementsChecker:
    """Check system requirements and provide memory warnings."""
    
    def __init__(self):
        self.total_memory_gb = psutil.virtual_memory().total / (1024**3)
        self.available_memory_gb = psutil.virtual_memory().available / (1024**3)
        self.cpu_count = psutil.cpu_count()
    
    def get_system_info(self) -> Dict[str, Any]:
        """Get comprehensive system information."""
        memory = psutil.virtual_memory()
        return {
            'total_memory_gb': memory.total / (1024**3),
            'available_memory_gb': memory.available / (1024**3),
            'used_memory_gb': memory.used / (1024**3),
            'memory_percent': memory.percent,
            'cpu_count': psutil.cpu_count(),
            'cpu_freq_mhz': psutil.cpu_freq().current if psutil.cpu_freq() else None,
            'platform': os.name,
        }
    
    def check_dataset_compatibility(self, dataset_size: str, dataset_memory_mb: float) -> Dict[str, Any]:
        """Check if a dataset is compatible with current system."""
        # System requirements by dataset size
        requirements = {
            'small': {'min_ram_gb': 8, 'recommended_ram_gb': 16},
            'medium': {'min_ram_gb': 16, 'recommended_ram_gb': 32},
            'large': {'min_ram_gb': 32, 'recommended_ram_gb': 64},
            'very_large': {'min_ram_gb': 64, 'recommended_ram_gb': 128},
        }
        
        if dataset_size not in requirements:
            raise ValueError(f"Unknown dataset size: {dataset_size}")
        
        req = requirements[dataset_size]
        dataset_memory_gb = dataset_memory_mb / 1024
        
        # Check compatibility
        compatible = self.total_memory_gb >= req['min_ram_gb']
        recommended = self.total_memory_gb >= req['recommended_ram_gb']
        
        # Calculate memory usage percentage
        memory_usage_percent = (dataset_memory_gb / self.total_memory_gb) * 100
        
        return {
            'compatible': compatible,
            'recommended': recommended,
            'dataset_memory_gb': dataset_memory_gb,
            'system_memory_gb': self.total_memory_gb,
            'memory_usage_percent': memory_usage_percent,
            'min_required_gb': req['min_ram_gb'],
            'recommended_gb': req['recommended_ram_gb'],
            'warnings': self._generate_warnings(compatible, recommended, memory_usage_percent, dataset_size)
        }
    
    def _generate_warnings(self, compatible: bool, recommended: bool, 
                          memory_usage_percent: float, dataset_size: str) -> list:
        """Generate appropriate warnings based on system capabilities."""
        warnings_list = []
        
        if not compatible:
            warnings_list.append(
                f"⚠️ CRITICAL: System has insufficient memory for {dataset_size} datasets. "
                f"Consider using backed/on-disk mode or upgrading to {self._get_min_ram(dataset_size)}+ GB RAM."
            )
        elif not recommended:
            warnings_list.append(
                f"⚠️ WARNING: System memory is below recommended for {dataset_size} datasets. "
                f"Performance may be degraded. Consider upgrading to {self._get_recommended_ram(dataset_size)}+ GB RAM."
            )
        
        if memory_usage_percent > 80:
            warnings_list.append(
                f"⚠️ WARNING: Dataset will use {memory_usage_percent:.1f}% of system memory. "
                "Consider closing other applications or using memory-efficient modes."
            )
        elif memory_usage_percent > 50:
            warnings_list.append(
                f"ℹ️ INFO: Dataset will use {memory_usage_percent:.1f}% of system memory. "
                "Monitor system performance during analysis."
            )
        
        return warnings_list
    
    def _get_min_ram(self, dataset_size: str) -> int:
        """Get minimum RAM requirement for dataset size."""
        requirements = {
            'small': 8, 'medium': 16, 'large': 32, 'very_large': 64
        }
        return requirements.get(dataset_size, 8)
    
    def _get_recommended_ram(self, dataset_size: str) -> int:
        """Get recommended RAM requirement for dataset size."""
        requirements = {
            'small': 16, 'medium': 32, 'large': 64, 'very_large': 128
        }
        return requirements.get(dataset_size, 16)
    
    def recommend_datasets(self) -> Dict[str, list]:
        """Recommend datasets based on current system capabilities."""
        recommendations = {
            'safe': [],
            'caution': [],
            'not_recommended': []
        }
        
        # Dataset memory requirements (in MB)
        dataset_requirements = {
            'pbmc3k': {'size': 'small', 'memory_mb': 350},
            'pbmc68k': {'size': 'medium', 'memory_mb': 8500},
            'mouse_brain_1_3m': {'size': 'large', 'memory_mb': 140000},
            'synthetic_small': {'size': 'small', 'memory_mb': 8},
            'synthetic_medium': {'size': 'medium', 'memory_mb': 200},
            'synthetic_large': {'size': 'large', 'memory_mb': 6000},
        }
        
        for dataset_id, req in dataset_requirements.items():
            compatibility = self.check_dataset_compatibility(req['size'], req['memory_mb'])
            
            if compatibility['compatible'] and compatibility['recommended']:
                recommendations['safe'].append(dataset_id)
            elif compatibility['compatible']:
                recommendations['caution'].append(dataset_id)
            else:
                recommendations['not_recommended'].append(dataset_id)
        
        return recommendations
    
    def print_system_report(self):
        """Print a comprehensive system report."""
        info = self.get_system_info()
        recommendations = self.recommend_datasets()
        
        print("🖥️ PySEE System Requirements Report")
        print("=" * 50)
        print(f"Total Memory: {info['total_memory_gb']:.1f} GB")
        print(f"Available Memory: {info['available_memory_gb']:.1f} GB")
        print(f"Used Memory: {info['used_memory_gb']:.1f} GB ({info['memory_percent']:.1f}%)")
        print(f"CPU Cores: {info['cpu_count']}")
        if info['cpu_freq_mhz']:
            print(f"CPU Frequency: {info['cpu_freq_mhz']:.0f} MHz")
        
        print(f"\n📊 Dataset Recommendations:")
        print(f"✅ Safe to use: {', '.join(recommendations['safe']) if recommendations['safe'] else 'None'}")
        print(f"⚠️ Use with caution: {', '.join(recommendations['caution']) if recommendations['caution'] else 'None'}")
        print(f"❌ Not recommended: {', '.join(recommendations['not_recommended']) if recommendations['not_recommended'] else 'None'}")
        
        # Memory recommendations
        if info['total_memory_gb'] < 16:
            print(f"\n💡 Memory Upgrade Recommendations:")
            print(f"   - Current: {info['total_memory_gb']:.1f} GB")
            print(f"   - Recommended: 16+ GB for medium datasets")
            print(f"   - Optimal: 32+ GB for large datasets")
    
    def warn_user(self, dataset_id: str, dataset_size: str, dataset_memory_mb: float):
        """Warn user about potential memory issues with a specific dataset."""
        compatibility = self.check_dataset_compatibility(dataset_size, dataset_memory_mb)
        
        if compatibility['warnings']:
            print(f"\n⚠️ Memory Warning for {dataset_id}:")
            for warning in compatibility['warnings']:
                print(f"   {warning}")
            
            # Suggest alternatives
            if not compatibility['compatible']:
                print(f"\n💡 Alternatives:")
                print(f"   - Use backed/on-disk mode: adata = ad.read_h5ad(file, backed='r')")
                print(f"   - Subsample the dataset to reduce memory usage")
                print(f"   - Use a smaller dataset for testing")
            
            return False
        
        return True


def check_system_requirements(dataset_id: str, dataset_size: str, dataset_memory_mb: float) -> bool:
    """Convenience function to check system requirements for a dataset."""
    checker = SystemRequirementsChecker()
    return checker.warn_user(dataset_id, dataset_size, dataset_memory_mb)


def get_system_info() -> Dict[str, Any]:
    """Get current system information."""
    checker = SystemRequirementsChecker()
    return checker.get_system_info()


def print_system_report():
    """Print a comprehensive system report."""
    checker = SystemRequirementsChecker()
    checker.print_system_report()


if __name__ == "__main__":
    # Test the system requirements checker
    print_system_report()
    
    # Test with a specific dataset
    print("\n" + "="*50)
    print("Testing dataset compatibility...")
    
    checker = SystemRequirementsChecker()
    checker.warn_user('pbmc68k', 'medium', 8500)
