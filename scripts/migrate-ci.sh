#!/bin/bash
# CI Migration Script
# Migrates from 5 overlapping workflows to 2 streamlined workflows

set -e

echo "🔄 PySEE CI Migration Script"
echo "============================="

# Check if we're in the right directory
if [ ! -f "pyproject.toml" ]; then
    echo "❌ Please run this script from the PySEE root directory"
    exit 1
fi

# Create backup directory
echo "📁 Creating backup directory..."
mkdir -p .github/workflows/backup
BACKUP_DIR=".github/workflows/backup/$(date +%Y%m%d_%H%M%S)"
mkdir -p "$BACKUP_DIR"

# Backup current workflows
echo "💾 Backing up current workflows..."
cp .github/workflows/ci.yml "$BACKUP_DIR/ci-backup.yml" 2>/dev/null || true
cp .github/workflows/ci-simple.yml "$BACKUP_DIR/ci-simple-backup.yml" 2>/dev/null || true
cp .github/workflows/ci-comprehensive.yml "$BACKUP_DIR/ci-comprehensive-backup.yml" 2>/dev/null || true
cp .github/workflows/pre-commit.yml "$BACKUP_DIR/pre-commit-backup.yml" 2>/dev/null || true

echo "✅ Backed up to: $BACKUP_DIR"

# Remove old workflows
echo "🗑️ Removing old workflows..."
rm -f .github/workflows/ci.yml
rm -f .github/workflows/ci-simple.yml
rm -f .github/workflows/ci-comprehensive.yml
rm -f .github/workflows/pre-commit.yml

# Deploy new workflows
echo "🚀 Deploying streamlined workflows..."
mv .github/workflows/ci-streamlined.yml .github/workflows/ci.yml
mv .github/workflows/pre-commit-streamlined.yml .github/workflows/pre-commit.yml

echo "✅ Streamlined workflows deployed"

# Verify workflows
echo "🔍 Verifying workflows..."
if [ -f ".github/workflows/ci.yml" ] && [ -f ".github/workflows/pre-commit.yml" ]; then
    echo "✅ Workflows verified"
else
    echo "❌ Workflow verification failed"
    exit 1
fi

# Show current workflow structure
echo ""
echo "📋 Current workflow structure:"
echo "=============================="
ls -la .github/workflows/

echo ""
echo "🎉 CI Migration Complete!"
echo "========================="
echo ""
echo "New structure:"
echo "  - ci.yml: Main CI workflow (quality, tests, build, security, release)"
echo "  - pre-commit.yml: Pre-commit hook validation"
echo "  - release.yml: Unchanged (tag-based releases)"
echo ""
echo "Backup location: $BACKUP_DIR"
echo ""
echo "Next steps:"
echo "  1. Commit and push changes"
echo "  2. Create a test PR to verify workflows"
echo "  3. Check CI execution times and results"
echo "  4. Remove backup files after confirmation"
echo ""
echo "To rollback if needed:"
echo "  cp $BACKUP_DIR/* .github/workflows/"
