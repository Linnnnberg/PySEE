@echo off
REM Comprehensive check script for PySEE (Windows)
REM Run this before committing to catch issues early

echo 🔍 PySEE Comprehensive Check Script
echo ==================================

REM Check if we're in the right directory
if not exist "pyproject.toml" (
    echo ❌ Please run this script from the PySEE root directory
    exit /b 1
)

echo 📦 Installing/updating dependencies...
pip install -q -r requirements-ci.txt
pip install -q -e .

echo.
echo 🎨 Running Black formatting check...
black --check --diff --line-length=100 pysee/ test_pysee.py example.py
if %errorlevel% neq 0 (
    echo ❌ Black formatting failed
    exit /b 1
)
echo ✅ Black formatting

echo.
echo 📋 Running isort import check...
isort --check-only --profile=black --line-length=100 pysee/ test_pysee.py example.py
if %errorlevel% neq 0 (
    echo ❌ Import sorting failed
    exit /b 1
)
echo ✅ Import sorting

echo.
echo 🔍 Running flake8 linting...
flake8 pysee/ test_pysee.py example.py --count --max-line-length=100 --extend-ignore=E203,W503,E501,F401 --statistics
if %errorlevel% neq 0 (
    echo ❌ Flake8 linting failed
    exit /b 1
)
echo ✅ Flake8 linting

echo.
echo 🔬 Running mypy type checking...
mypy pysee/ --ignore-missing-imports --no-strict-optional --warn-return-any --no-error-summary
if %errorlevel% neq 0 (
    echo ❌ MyPy type checking failed
    exit /b 1
)
echo ✅ MyPy type checking

echo.
echo 🧪 Running tests...
pytest test_pysee.py -v --timeout=60
if %errorlevel% neq 0 (
    echo ❌ Tests failed
    exit /b 1
)
pytest tests/ -v --timeout=60 --maxfail=3
if %errorlevel% neq 0 (
    echo ❌ Tests failed
    exit /b 1
)
echo ✅ Tests

echo.
echo 📖 Testing example script...
python example.py
if %errorlevel% neq 0 (
    echo ❌ Example script failed
    exit /b 1
)
echo ✅ Example script

echo.
echo 🏗️ Testing package build...
python -m build
if %errorlevel% neq 0 (
    echo ❌ Package build failed
    exit /b 1
)
echo ✅ Package build

echo.
echo 🎉 All checks passed! Ready to commit.
