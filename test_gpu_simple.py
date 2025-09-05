"""
Simple GPU test to verify CuPy installation and functionality.
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))


def test_gpu_installation():
    """Test GPU installation and functionality."""
    print("🔍 GPU Installation Test")
    print("=" * 40)

    # Test 1: Import CuPy
    try:
        import cupy as cp

        print("✅ CuPy import successful")
        print(f"   CuPy version: {cp.__version__}")
    except ImportError as e:
        print(f"❌ CuPy import failed: {e}")
        return False

    # Test 2: Check CUDA availability
    try:
        cuda_available = cp.cuda.is_available()
        print(f"✅ CUDA available: {cuda_available}")

        if cuda_available:
            gpu_count = cp.cuda.runtime.getDeviceCount()
            print(f"   GPU count: {gpu_count}")

            # Get GPU info
            device = cp.cuda.Device()
            mem_info = device.mem_info
            total_mem = mem_info[1] / 1024**3
            free_mem = mem_info[0] / 1024**3
            print(f"   GPU memory: {free_mem:.1f} GB free / {total_mem:.1f} GB total")
        else:
            print("❌ CUDA not available")
            return False

    except Exception as e:
        print(f"❌ CUDA check failed: {e}")
        return False

    # Test 3: Basic GPU operations
    try:
        print("\n🧪 Testing basic GPU operations...")

        # Test array creation
        cpu_array = cp.array([1, 2, 3, 4, 5])
        print("✅ GPU array creation successful")

        # Test simple computation
        result = cpu_array * 2
        print("✅ GPU computation successful")

        # Test transfer back to CPU
        cpu_result = cp.asnumpy(result)
        print("✅ CPU transfer successful")

    except Exception as e:
        print(f"❌ Basic GPU operations failed: {e}")
        return False

    # Test 4: Advanced operations (SVD)
    try:
        print("\n🧪 Testing advanced GPU operations (SVD)...")

        import numpy as np

        # Create small test matrix
        test_matrix = np.random.rand(50, 50).astype(np.float32)
        gpu_matrix = cp.asarray(test_matrix)

        # Test SVD
        u, s, v = cp.linalg.svd(gpu_matrix, full_matrices=False)
        print("✅ GPU SVD computation successful")

        # Test transfer back
        cpu_u = cp.asnumpy(u)
        cpu_s = cp.asnumpy(s)
        cpu_v = cp.asnumpy(v)
        print("✅ GPU SVD transfer to CPU successful")

    except Exception as e:
        print(f"❌ Advanced GPU operations failed: {e}")
        print("   This is likely due to missing CUDA libraries (cusolver, etc.)")
        return False

    print("\n🎉 All GPU tests passed!")
    return True


def test_gpu_vs_cpu_performance():
    """Test GPU vs CPU performance."""
    print("\n⚡ GPU vs CPU Performance Test")
    print("=" * 40)

    try:
        import cupy as cp
        import numpy as np
        import time

        if not cp.cuda.is_available():
            print("❌ CUDA not available, skipping performance test")
            return

        # Create test matrix
        size = 1000
        print(f"Testing with {size}x{size} matrix...")

        cpu_matrix = np.random.rand(size, size).astype(np.float32)

        # CPU computation
        start_time = time.time()
        cpu_result = np.linalg.svd(cpu_matrix, full_matrices=False)
        cpu_time = time.time() - start_time
        print(f"CPU SVD: {cpu_time:.2f}s")

        # GPU computation
        gpu_matrix = cp.asarray(cpu_matrix)

        start_time = time.time()
        gpu_result = cp.linalg.svd(gpu_matrix, full_matrices=False)
        gpu_time = time.time() - start_time
        print(f"GPU SVD: {gpu_time:.2f}s")

        speedup = cpu_time / gpu_time
        print(f"GPU Speedup: {speedup:.2f}x")

        if speedup > 1:
            print("✅ GPU is faster than CPU")
        else:
            print("⚠️ GPU is not faster than CPU (may be due to small matrix size)")

    except Exception as e:
        print(f"❌ Performance test failed: {e}")


def main():
    """Run GPU tests."""
    success = test_gpu_installation()

    if success:
        test_gpu_vs_cpu_performance()

        print("\n💡 GPU Status Summary:")
        print("✅ CuPy is properly installed and working")
        print("✅ GPU acceleration is available for PySEE")
        print("💡 Consider using GPU for heavy computations in PySEE")
    else:
        print("\n💡 GPU Status Summary:")
        print("⚠️ CuPy is installed but GPU operations are not working")
        print("💡 This is likely due to missing CUDA libraries")
        print("💡 For now, stick with CPU-based PySEE operations")
        print("💡 GPU acceleration can be added later when CUDA is properly configured")


if __name__ == "__main__":
    main()
