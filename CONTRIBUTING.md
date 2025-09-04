# Contributing to PySEE

Thank you for your interest in contributing to PySEE! This document provides guidelines and information for contributors.

## Getting Started

### Prerequisites
- Python 3.8 or higher
- Git
- Basic knowledge of single-cell analysis and Python

### Development Setup

1. **Fork and clone the repository**
   ```bash
   git clone https://github.com/YOUR_USERNAME/PySEE.git
   cd PySEE
   ```

2. **Create a virtual environment**
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   pip install -e .
   ```

4. **Run tests to ensure everything works**
   ```bash
   python test_pysee.py
   python example.py
   ```

## Development Workflow

### Branch Strategy
- `main`: Production-ready code
- `develop`: Integration branch for features
- `feature/*`: Feature branches
- `bugfix/*`: Bug fix branches
- `hotfix/*`: Critical bug fixes

### Making Changes

1. **Create a feature branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes**
   - Follow the coding style guidelines
   - Add tests for new functionality
   - Update documentation as needed

3. **Test your changes**
   ```bash
   # Run linting
   flake8 .
   black --check .
   mypy pysee/

   # Run tests
   python test_pysee.py
   ```

4. **Commit your changes**
   ```bash
   git add .
   git commit -m "Add: brief description of changes"
   ```

5. **Push and create a pull request**
   ```bash
   git push origin feature/your-feature-name
   ```

## Coding Standards

### Python Style
- Follow PEP 8
- Use Black for code formatting
- Use type hints where appropriate
- Keep functions small and focused
- Write descriptive docstrings

### Code Formatting
```bash
# Format code with Black
black .

# Check formatting
black --check .
```

### Linting
```bash
# Run flake8
flake8 . --max-line-length=127

# Run mypy for type checking
mypy pysee/ --ignore-missing-imports
```

### Documentation
- Use Google-style docstrings
- Include type hints in function signatures
- Update README.md for user-facing changes
- Update CHANGELOG.md for significant changes

## Testing

### Running Tests
```bash
# Run basic functionality test
python test_pysee.py

# Run example with real data
python example.py

# Run with pytest (when available)
pytest tests/ -v
```

### Writing Tests
- Test new functionality thoroughly
- Include edge cases and error conditions
- Use descriptive test names
- Keep tests simple and focused

## Pull Request Process

### Before Submitting
- [ ] Code follows style guidelines
- [ ] Tests pass locally
- [ ] Documentation is updated
- [ ] CHANGELOG.md is updated (if applicable)
- [ ] No merge conflicts

### PR Description
- Clearly describe what the PR does
- Reference any related issues
- Include screenshots for UI changes
- List any breaking changes

### Review Process
1. Automated CI checks must pass
2. At least one maintainer review required
3. Address all review comments
4. Maintainer will merge when ready

## Issue Guidelines

### Bug Reports
- Use the bug report template
- Include steps to reproduce
- Provide environment details
- Include error messages and screenshots

### Feature Requests
- Use the feature request template
- Describe the use case clearly
- Consider implementation complexity
- Check for existing similar requests

## Release Process

### Version Numbering
- Follow semantic versioning (MAJOR.MINOR.PATCH)
- Update version in `setup.py` and `pysee/__init__.py`
- Update CHANGELOG.md

### Creating a Release
1. Update version numbers
2. Update CHANGELOG.md
3. Create a git tag: `git tag v0.2.0`
4. Push tag: `git push origin v0.2.0`
5. GitHub Actions will automatically create a release

## Project Structure

```
pysee/
├── core/           # Core dashboard and data handling
├── panels/         # Visualization panels
├── cli/            # Command-line interface
├── utils/          # Utility functions
└── __init__.py     # Package initialization

.github/
├── workflows/      # GitHub Actions CI/CD
├── ISSUE_TEMPLATE/ # Issue templates
└── pull_request_template.md

docs/               # Documentation (future)
tests/              # Test suite (future)
```

## Communication

### Getting Help
- Check existing issues and discussions
- Create a new issue for bugs or feature requests
- Use GitHub Discussions for questions

### Code of Conduct
- Be respectful and inclusive
- Focus on constructive feedback
- Help others learn and grow
- Follow the [Contributor Covenant](https://www.contributor-covenant.org/)

## Recognition

Contributors will be recognized in:
- CONTRIBUTORS.md file
- Release notes
- Project documentation

Thank you for contributing to PySEE! 🎉
