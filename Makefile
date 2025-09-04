.PHONY: help install test lint format clean build

help: ## Show this help message
	@echo "Available commands:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

install: ## Install package in development mode
	pip install -e .

install-dev: ## Install development dependencies
	pip install -r requirements.txt
	pip install pre-commit

test: ## Run tests
	python test_pysee.py
	python example.py

lint: ## Run linting
	flake8 . --max-line-length=127
	black --check .
	mypy pysee/ --ignore-missing-imports

format: ## Format code
	black .
	isort .

clean: ## Clean build artifacts
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	find . -type d -name __pycache__ -delete
	find . -type f -name "*.pyc" -delete

build: ## Build package
	python -m build

check-build: ## Check built package
	twine check dist/*

setup-pre-commit: ## Setup pre-commit hooks
	pre-commit install

run-pre-commit: ## Run pre-commit on all files
	pre-commit run --all-files

ci: lint test ## Run CI checks locally
