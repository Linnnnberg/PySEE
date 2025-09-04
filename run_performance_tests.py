"""
Quick performance test runner for PySEE.

This script runs a quick performance test to verify PySEE is working
efficiently with different dataset sizes.
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from tests.performance.run_performance_tests import main

if __name__ == "__main__":
    exit(main())
