# Git Workflow Strategy for PySEE

## Overview

PySEE follows a **feature branch workflow** with **protected main branch** and **automated CI/CD**. This ensures code quality, collaboration safety, and maintainable development practices.

## Branch Strategy

### Core Branches

- **`main`** - Production-ready code, always stable
- **`develop`** - Integration branch for features, staging area
- **`feature/*`** - Individual feature development branches
- **`hotfix/*`** - Critical bug fixes for production
- **`release/*`** - Release preparation branches

### Branch Naming Convention

```
feature/task-description
hotfix/critical-bug-description
release/v0.2.0
```

Examples:
- `feature/heatmap-panel`
- `feature/igv-integration`
- `hotfix/umap-selection-bug`
- `release/v0.2.0`

## Development Workflow

### 1. Starting a New Task

```bash
# Ensure you're on develop and it's up to date
git checkout develop
git pull origin develop

# Create a new feature branch
git checkout -b feature/task-description

# Example:
git checkout -b feature/heatmap-panel
```

### 2. Development Process

```bash
# Make your changes
# ... code changes ...

# Stage and commit with descriptive messages
git add .
git commit -m "Add heatmap panel with interactive selection

- Implement HeatmapPanel class
- Add gene expression heatmap visualization
- Support for row/column clustering
- Integrate with existing panel linking system"

# Push the feature branch
git push origin feature/heatmap-panel
```

### 3. Creating Pull Requests

1. **GitHub Web Interface**:
   - Go to your repository
   - Click "Compare & pull request"
   - Set base: `develop` ← compare: `feature/heatmap-panel`
   - Add detailed description
   - Request reviews from collaborators

2. **Pull Request Template**:
   ```markdown
   ## Description
   Brief description of changes

   ## Type of Change
   - [ ] Bug fix
   - [ ] New feature
   - [ ] Breaking change
   - [ ] Documentation update

   ## Testing
   - [ ] Local tests pass
   - [ ] CI tests pass
   - [ ] Manual testing completed

   ## Checklist
   - [ ] Code follows project style
   - [ ] Self-review completed
   - [ ] Documentation updated
   - [ ] No breaking changes (or documented)
   ```

### 4. Code Review Process

1. **Automated Checks** (must pass):
   - ✅ CI/CD pipeline (linting, testing, building)
   - ✅ All tests pass
   - ✅ Code coverage maintained

2. **Manual Review**:
   - Code quality and style
   - Architecture decisions
   - Performance implications
   - Documentation completeness

3. **Approval Requirements**:
   - At least 1 approval from maintainer
   - All CI checks green
   - No merge conflicts

### 5. Merging Strategy

```bash
# After PR approval, merge via GitHub (recommended)
# OR merge locally:

git checkout develop
git pull origin develop
git merge --no-ff feature/heatmap-panel
git push origin develop

# Clean up feature branch
git branch -d feature/heatmap-panel
git push origin --delete feature/heatmap-panel
```

## Release Workflow

### 1. Release Preparation

```bash
# Create release branch from develop
git checkout develop
git pull origin develop
git checkout -b release/v0.2.0

# Update version numbers, changelog
# ... make release preparations ...

git add .
git commit -m "Prepare release v0.2.0"
git push origin release/v0.2.0
```

### 2. Release Process

```bash
# Create PR: release/v0.2.0 → main
# After approval and merge:

# Tag the release
git checkout main
git pull origin main
git tag -a v0.2.0 -m "Release version 0.2.0"
git push origin v0.2.0

# Merge back to develop
git checkout develop
git merge main
git push origin develop
```

## Hotfix Workflow

For critical production bugs:

```bash
# Create hotfix from main
git checkout main
git pull origin main
git checkout -b hotfix/critical-bug

# Fix the bug
# ... make minimal changes ...

git add .
git commit -m "Fix critical selection bug in UMAP panel"
git push origin hotfix/critical-bug

# Create PR: hotfix/critical-bug → main
# After merge, also merge to develop
```

## Branch Protection Rules

### Main Branch Protection
- ✅ Require pull request reviews
- ✅ Require status checks (CI/CD)
- ✅ Require branches to be up to date
- ✅ Restrict pushes to main
- ✅ Allow force pushes: ❌
- ✅ Allow deletions: ❌

### Develop Branch Protection
- ✅ Require pull request reviews
- ✅ Require status checks (CI/CD)
- ✅ Allow force pushes: ❌

## CI/CD Integration

### Automated Checks
Every PR automatically runs:
- **Linting**: flake8, black, mypy
- **Testing**: pytest with coverage
- **Building**: Package installation test
- **Multi-Python**: 3.9, 3.10, 3.11, 3.12

### Required Status Checks
- `ci / lint-and-test (3.9)`
- `ci / lint-and-test (3.10)`
- `ci / lint-and-test (3.11)`
- `ci / lint-and-test (3.12)`

## Best Practices

### Commit Messages
Follow [Conventional Commits](https://www.conventionalcommits.org/):

```
feat: add heatmap panel with clustering
fix: resolve UMAP selection propagation bug
docs: update installation instructions
test: add unit tests for ViolinPanel
refactor: simplify panel linking logic
```

### Branch Hygiene
- ✅ Keep feature branches focused and small
- ✅ Rebase before creating PR (optional)
- ✅ Delete merged branches promptly
- ✅ Use descriptive branch names
- ❌ Don't commit directly to main/develop
- ❌ Don't force push to shared branches

### Code Quality
- ✅ Write tests for new features
- ✅ Update documentation
- ✅ Follow existing code style
- ✅ Keep commits atomic and focused
- ✅ Write clear commit messages

## Emergency Procedures

### Reverting a Bad Merge
```bash
# Find the merge commit
git log --oneline --merges

# Revert the merge
git revert -m 1 <merge-commit-hash>
git push origin main
```

### Recovering from Force Push
```bash
# Find the lost commit
git reflog

# Recover the commit
git checkout <commit-hash>
git checkout -b recovery-branch
```

## Tools and Automation

### Pre-commit Hooks
Automatically run on every commit:
- Code formatting (black)
- Linting (flake8)
- Type checking (mypy)
- Security checks

### GitHub Actions
- **CI Pipeline**: Automated testing and quality checks
- **Release Pipeline**: Automated package building and publishing
- **Dependency Updates**: Automated security updates

## Collaboration Guidelines

### For Contributors
1. Fork the repository
2. Create feature branch from `develop`
3. Make changes and test locally
4. Create PR to `develop`
5. Address review feedback
6. Wait for approval and merge

### For Maintainers
1. Review PRs promptly
2. Ensure CI passes before merging
3. Maintain code quality standards
4. Update documentation as needed
5. Coordinate releases

## Example Workflow

```bash
# 1. Start new feature
git checkout develop
git pull origin develop
git checkout -b feature/igv-integration

# 2. Develop feature
# ... make changes ...
git add .
git commit -m "feat: add IGV genome browser integration"
git push origin feature/igv-integration

# 3. Create PR on GitHub
# 4. Address review feedback
# 5. After approval, merge via GitHub
# 6. Clean up
git checkout develop
git pull origin develop
git branch -d feature/igv-integration
```

This workflow ensures:
- ✅ **Code Quality**: All changes reviewed and tested
- ✅ **Collaboration**: Clear process for multiple developers
- ✅ **Stability**: Main branch always production-ready
- ✅ **Traceability**: Clear history of all changes
- ✅ **Automation**: CI/CD handles quality checks
