"""
System requirements checker for PySEE.

This script checks your system's capabilities and provides recommendations
for which datasets you can safely use with PySEE.
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from pysee.utils.system_requirements import SystemRequirementsChecker, print_system_report


def main():
    """Check system requirements and provide recommendations."""
    print("🔍 PySEE System Requirements Checker")
    print("=" * 60)

    try:
        # Print comprehensive system report
        print_system_report()

        # Test specific datasets
        print("\n" + "=" * 60)
        print("🧪 Testing Dataset Compatibility")
        print("=" * 60)

        checker = SystemRequirementsChecker()

        # Test datasets
        test_datasets = [
            ("pbmc3k", "small", 350),
            ("pbmc68k", "medium", 8500),
            ("mouse_brain_1_3m", "large", 140000),
        ]

        for dataset_id, size, memory_mb in test_datasets:
            print(f"\n📊 Testing {dataset_id} ({size} dataset, {memory_mb} MB):")
            compatible = checker.warn_user(dataset_id, size, memory_mb)

            if compatible:
                print(f"   ✅ {dataset_id} is compatible with your system")
            else:
                print(f"   ⚠️ {dataset_id} may cause memory issues")

        # Provide upgrade recommendations
        print("\n" + "=" * 60)
        print("💡 Upgrade Recommendations")
        print("=" * 60)

        info = checker.get_system_info()

        if info["total_memory_gb"] < 16:
            print("🔧 Memory Upgrade Recommended:")
            print(f"   Current: {info['total_memory_gb']:.1f} GB")
            print(f"   Recommended: 16+ GB for medium datasets (PBMC 68K)")
            print(f"   Optimal: 32+ GB for large datasets (Mouse Brain 1.3M)")
            print(f"   Workstation: 64+ GB for very large datasets")
        elif info["total_memory_gb"] < 32:
            print("🔧 Memory Upgrade Optional:")
            print(f"   Current: {info['total_memory_gb']:.1f} GB")
            print(f"   Recommended: 32+ GB for large datasets (Mouse Brain 1.3M)")
            print(f"   Workstation: 64+ GB for very large datasets")
        else:
            print("✅ Your system has sufficient memory for most PySEE datasets!")

        print("\n💡 Tips for Memory-Efficient Usage:")
        print("   - Use backed/on-disk mode for large datasets")
        print("   - Subsample datasets to reduce memory usage")
        print("   - Close other applications when working with large datasets")
        print("   - Consider using cloud computing for very large datasets")

        return 0

    except Exception as e:
        print(f"❌ System check failed: {e}")
        import traceback

        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(main())
