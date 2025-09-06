#!/bin/bash
# Comprehensive check script for PySEE
# Run this before committing to catch issues early

set -e  # Exit on any error

echo "🔍 PySEE Comprehensive Check Script"
echo "=================================="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print status
print_status() {
    if [ $1 -eq 0 ]; then
        echo -e "${GREEN}✅ $2${NC}"
    else
        echo -e "${RED}❌ $2${NC}"
        exit 1
    fi
}

# Check if we're in the right directory
if [ ! -f "pyproject.toml" ]; then
    echo -e "${RED}❌ Please run this script from the PySEE root directory${NC}"
    exit 1
fi

echo "📦 Installing/updating dependencies..."
pip install -q -r requirements-ci.txt
pip install -q -e .

echo ""
echo "🎨 Running Black formatting check..."
black --check --diff --line-length=100 pysee/ test_pysee.py example.py
print_status $? "Black formatting"

echo ""
echo "📋 Running isort import check..."
isort --check-only --profile=black --line-length=100 pysee/ test_pysee.py example.py
print_status $? "Import sorting"

echo ""
echo "🔍 Running flake8 linting..."
flake8 pysee/ test_pysee.py example.py --count --max-line-length=100 --extend-ignore=E203,W503,E501,F401 --statistics
print_status $? "Flake8 linting"

echo ""
echo "🔬 Running mypy type checking..."
mypy pysee/ --ignore-missing-imports --no-strict-optional --warn-return-any --no-error-summary
print_status $? "MyPy type checking"

echo ""
echo "🧪 Running tests..."
pytest test_pysee.py -v --timeout=60
pytest tests/ -v --timeout=60 --maxfail=3
print_status $? "Tests"

echo ""
echo "📖 Testing example script..."
python example.py
print_status $? "Example script"

echo ""
echo "🏗️ Testing package build..."
python -m build
print_status $? "Package build"

echo ""
echo -e "${GREEN}🎉 All checks passed! Ready to commit.${NC}"
