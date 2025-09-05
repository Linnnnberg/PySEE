# PySEE Makefile
# Easy commands for development and CI

.PHONY: help install test lint format type-check check-all clean build

# Default target
help:
	@echo "PySEE Development Commands"
	@echo "========================="
	@echo ""
	@echo "Setup:"
	@echo "  install     Install dependencies"
	@echo "  install-dev Install development dependencies"
	@echo ""
	@echo "Code Quality:"
	@echo "  format      Format code with Black and isort"
	@echo "  lint        Run flake8 linting"
	@echo "  type-check  Run mypy type checking"
	@echo "  check-all   Run all quality checks"
	@echo ""
	@echo "Testing:"
	@echo "  test        Run tests"
	@echo "  test-fast   Run fast tests only"
	@echo "  test-perf   Run performance tests"
	@echo ""
	@echo "Build:"
	@echo "  build       Build package"
	@echo "  clean       Clean build artifacts"
	@echo ""
	@echo "CI:"
	@echo "  ci-local    Run local CI checks"
	@echo "  pre-commit  Run pre-commit hooks"

# Installation
install:
	pip install -e .

install-dev:
	pip install -r requirements-ci.txt
	pip install -e .
	pre-commit install

# Code formatting
format:
	@echo "🎨 Formatting code with Black..."
	black --line-length=100 pysee/ test_pysee.py example.py
	@echo "📋 Sorting imports with isort..."
	isort --profile=black --line-length=100 pysee/ test_pysee.py example.py
	@echo "✅ Code formatted!"

# Linting
lint:
	@echo "🔍 Running flake8 linting..."
	flake8 pysee/ test_pysee.py example.py --count --max-line-length=100 --extend-ignore=E203,W503,E501,F401 --statistics

# Type checking
type-check:
	@echo "🔬 Running mypy type checking..."
	mypy pysee/ --ignore-missing-imports --no-strict-optional --warn-return-any --no-error-summary

# All quality checks
check-all: format lint type-check test
	@echo "🎉 All quality checks passed!"

# Testing
test:
	@echo "🧪 Running tests..."
	pytest test_pysee.py -v --timeout=60
	pytest tests/ -v --timeout=60 --maxfail=3

test-fast:
	@echo "🧪 Running fast tests..."
	pytest test_pysee.py -v --timeout=30

test-perf:
	@echo "🧪 Running performance tests..."
	pytest tests/performance/ -v --timeout=300

# Build
build:
	@echo "🏗️ Building package..."
	python -m build

clean:
	@echo "🧹 Cleaning build artifacts..."
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

# CI
ci-local:
	@echo "🔍 Running local CI checks..."
	@if [ -f "scripts/check-all.sh" ]; then \
		bash scripts/check-all.sh; \
	else \
		echo "❌ Check script not found. Running individual checks..."; \
		$(MAKE) check-all; \
	fi

pre-commit:
	@echo "🔍 Running pre-commit hooks..."
	pre-commit run --all-files

# Development helpers
dev-setup: install-dev
	@echo "🚀 Development environment ready!"
	@echo "Run 'make check-all' before committing"

# Quick fix for common issues
fix-format:
	@echo "🔧 Fixing formatting issues..."
	$(MAKE) format
	@echo "✅ Formatting fixed!"

fix-imports:
	@echo "🔧 Fixing import issues..."
	isort --profile=black --line-length=100 pysee/ test_pysee.py example.py
	@echo "✅ Imports fixed!"

# Check specific files
check-pysee:
	@echo "🔍 Checking pysee/ directory..."
	black --check --diff --line-length=100 pysee/
	isort --check-only --profile=black --line-length=100 pysee/
	flake8 pysee/ --max-line-length=100 --extend-ignore=E203,W503,E501,F401
	mypy pysee/ --ignore-missing-imports --no-strict-optional --warn-return-any --no-error-summary
