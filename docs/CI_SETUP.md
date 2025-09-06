# CI/CD Setup for PySEE

This document explains the comprehensive CI/CD setup to avoid formatting and type checking issues.

## 🎯 Problem Solved

**Before**: CI failures due to:
- Black formatting inconsistencies (line-length mismatch)
- MyPy type checking errors
- Import sorting issues
- Missing type annotations

**After**: Automated checks catch issues early and provide clear feedback.

## 📁 Files Created

### 1. Pre-commit Configuration
- **`.pre-commit-config.yaml`**: Runs checks before each commit
- **Installation**: `pre-commit install`
- **Usage**: Automatically runs on `git commit`

### 2. CI Workflows
- **`.github/workflows/ci.yml`**: Main CI (existing, updated)
- **`.github/workflows/ci-simple.yml`**: Simplified CI for quick checks
- **`.github/workflows/ci-comprehensive.yml`**: Full CI with security scans
- **`.github/workflows/pre-commit.yml`**: Pre-commit CI validation

### 3. Local Development Scripts
- **`scripts/check-all.sh`**: Unix/Linux check script
- **`scripts/check-all.bat`**: Windows check script
- **`Makefile`**: Easy commands for development

## 🚀 Quick Start

### Install Pre-commit Hooks
```bash
# Install pre-commit
pip install pre-commit

# Install hooks
pre-commit install

# Run on all files
pre-commit run --all-files
```

### Use Make Commands
```bash
# Install development dependencies
make install-dev

# Run all checks
make check-all

# Format code
make format

# Run tests
make test

# Clean build artifacts
make clean
```

### Use Check Scripts
```bash
# Unix/Linux
bash scripts/check-all.sh

# Windows
scripts/check-all.bat
```

## 🔧 Configuration Details

### Black Formatting
- **Line length**: 100 characters (matches CI)
- **Target Python**: 3.8+
- **Files**: All `.py` files

### isort Import Sorting
- **Profile**: Black-compatible
- **Line length**: 100 characters
- **Files**: All `.py` files

### flake8 Linting
- **Max line length**: 100 characters
- **Ignores**: E203, W503, E501, F401 (Black compatibility)
- **Max complexity**: 15

### MyPy Type Checking
- **Mode**: Non-strict (ignore missing imports)
- **Warnings**: Return any, unused configs
- **Target**: `pysee/` directory only

## 📋 Pre-commit Hooks

| Hook | Purpose | Files |
|------|---------|-------|
| `black` | Code formatting | `*.py` |
| `isort` | Import sorting | `*.py` |
| `flake8` | Linting | `*.py` |
| `mypy` | Type checking | `pysee/` |
| `trailing-whitespace` | Remove trailing spaces | All |
| `end-of-file-fixer` | Ensure newline at EOF | All |
| `check-yaml` | Validate YAML files | `*.yml`, `*.yaml` |
| `check-json` | Validate JSON files | `*.json` |
| `check-toml` | Validate TOML files | `*.toml` |
| `check-merge-conflict` | Detect merge conflicts | All |
| `check-added-large-files` | Prevent large files | All |
| `debug-statements` | Remove debug statements | `*.py` |
| `check-docstring-first` | Ensure docstrings | `*.py` |

## 🎨 CI Workflow Strategy

### 1. Fast Feedback (ci-simple.yml)
- **Purpose**: Quick checks for every PR
- **Duration**: ~5-10 minutes
- **Includes**: Format, lint, type, test, build

### 2. Comprehensive Checks (ci-comprehensive.yml)
- **Purpose**: Full validation with security scans
- **Duration**: ~15-20 minutes
- **Includes**: All checks + security + performance

### 3. Pre-commit Validation (pre-commit.yml)
- **Purpose**: Validate pre-commit configuration
- **Duration**: ~5 minutes
- **Includes**: Pre-commit hooks only

## 🛠️ Development Workflow

### Before Committing
```bash
# Option 1: Use pre-commit (automatic)
git add .
git commit -m "your message"  # Pre-commit runs automatically

# Option 2: Use make command
make check-all

# Option 3: Use check script
bash scripts/check-all.sh
```

### Before Pushing
```bash
# Run local CI
make ci-local

# Or use the check script
bash scripts/check-all.sh
```

### Fix Common Issues
```bash
# Fix formatting
make fix-format

# Fix imports
make fix-imports

# Check specific directory
make check-pysee
```

## 🔍 Troubleshooting

### Black Formatting Issues
```bash
# Check current formatting
black --check --diff --line-length=100 pysee/ test_pysee.py example.py

# Fix formatting
black --line-length=100 pysee/ test_pysee.py example.py
```

### MyPy Type Issues
```bash
# Check types
mypy pysee/ --ignore-missing-imports --no-strict-optional --warn-return-any --no-error-summary

# Common fixes:
# 1. Add type annotations: **kwargs: Any
# 2. Use type: ignore comments for complex cases
# 3. Cast types: int(self.config["value"])  # type: ignore[arg-type]
```

### Import Sorting Issues
```bash
# Check imports
isort --check-only --profile=black --line-length=100 pysee/ test_pysee.py example.py

# Fix imports
isort --profile=black --line-length=100 pysee/ test_pysee.py example.py
```

## 📊 Benefits

1. **Early Detection**: Issues caught before CI
2. **Consistent Formatting**: All code follows same style
3. **Type Safety**: MyPy catches type errors
4. **Import Organization**: Clean, sorted imports
5. **Fast Feedback**: Quick local checks
6. **CI Reliability**: Fewer CI failures
7. **Developer Experience**: Clear error messages

## 🎯 Best Practices

1. **Always run checks locally** before pushing
2. **Use pre-commit hooks** for automatic checking
3. **Fix formatting issues** immediately
4. **Add type annotations** for new functions
5. **Keep imports organized** with isort
6. **Use make commands** for common tasks
7. **Check CI logs** for detailed error messages

This setup ensures consistent code quality and prevents the formatting/type checking issues that were causing CI failures.
