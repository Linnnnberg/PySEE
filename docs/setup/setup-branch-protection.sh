#!/bin/bash

# GitHub CLI Branch Protection Setup for PySEE
echo "🔐 Setting up branch protection for PySEE..."

# Check if authenticated
if ! gh auth status >/dev/null 2>&1; then
    echo "❌ Please authenticate first: gh auth login"
    exit 1
fi

echo "✅ GitHub CLI authenticated"

# Apply branch protection
echo "🛡️ Applying branch protection rules..."

gh api repos/Linnnnberg/PySEE/branches/main/protection \
  --method PUT \
  --field required_status_checks='{"strict":true,"contexts":["test (3.9)","test (3.10)","test (3.11)","test (3.12)","build"]}' \
  --field enforce_admins=true \
  --field required_pull_request_reviews='{"required_approving_review_count":1,"dismiss_stale_reviews":true,"require_code_owner_reviews":true}' \
  --field restrictions=null \
  --field allow_force_pushes=false \
  --field allow_deletions=false \
  --field required_linear_history=true

if [ $? -eq 0 ]; then
    echo "✅ Branch protection applied successfully!"
else
    echo "❌ Failed to apply branch protection"
    exit 1
fi

# Commit and push CODEOWNERS
echo "📝 Adding CODEOWNERS file..."
git add .github/CODEOWNERS
git commit -m "feat: add CODEOWNERS for branch protection"
git push origin main

echo "🎉 Branch protection setup complete!"
echo "📋 Summary:"
echo "  - Status checks required: test (3.9-3.12), build"
echo "  - PR reviews required: 1 approval"
echo "  - CODEOWNERS: @Linnnnberg"
echo "  - Force pushes: disabled"
echo "  - Direct pushes: disabled"
