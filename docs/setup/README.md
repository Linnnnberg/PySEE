# GitHub Branch Protection Setup

This directory contains scripts and configuration files for setting up GitHub branch protection for the PySEE project.

## Files

### `branch-protection.json`
JSON configuration file for GitHub branch protection rules. This file defines:
- Required status checks (test suites for Python 3.9-3.12, build)
- Pull request review requirements (1 approval, code owner reviews)
- Branch restrictions (no force pushes, no deletions, linear history)

### `setup-branch-protection.sh`
Bash script to automatically apply branch protection rules using GitHub CLI.

## Usage

### Manual Setup
1. Authenticate with GitHub CLI:
   ```bash
   gh auth login
   ```

2. Apply branch protection:
   ```bash
   gh api repos/Linnnnberg/PySEE/branches/main/protection \
     --method PUT \
     --input branch-protection.json
   ```

### Automated Setup
1. Make the script executable:
   ```bash
   chmod +x setup-branch-protection.sh
   ```

2. Run the setup script:
   ```bash
   ./setup-branch-protection.sh
   ```

## Current Protection Rules

✅ **Status Checks Required:**
- test (3.9)
- test (3.10)
- test (3.11)
- test (3.12)
- build

✅ **Pull Request Requirements:**
- 1 approval required
- Code owner reviews required
- Stale reviews dismissed on push

✅ **Branch Restrictions:**
- Force pushes: disabled
- Deletions: disabled
- Linear history: required
- Direct pushes: disabled (must use PRs)

## CODEOWNERS

The `.github/CODEOWNERS` file specifies:
- `@Linnnnberg` as owner for all files
- Special protection for critical files:
  - `pyproject.toml`
  - `.github/workflows/`
  - `README.md`
  - `pysee/` (core code)
  - `tests/` (test code)

## Verification

Check current protection status:
```bash
gh api repos/Linnnnberg/PySEE/branches/main/protection
```

## Notes

- Branch protection is now active on the `main` branch
- All changes must go through pull requests
- Status checks must pass before merging
- Code owner approval required for protected files
