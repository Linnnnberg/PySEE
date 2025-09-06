@echo off
REM CI Migration Script (Windows)
REM Migrates from 5 overlapping workflows to 2 streamlined workflows

echo 🔄 PySEE CI Migration Script
echo =============================

REM Check if we're in the right directory
if not exist "pyproject.toml" (
    echo ❌ Please run this script from the PySEE root directory
    exit /b 1
)

REM Create backup directory
echo 📁 Creating backup directory...
mkdir .github\workflows\backup 2>nul
for /f "tokens=2 delims==" %%a in ('wmic OS Get localdatetime /value') do set "dt=%%a"
set "YY=%dt:~2,2%" & set "YYYY=%dt:~0,4%" & set "MM=%dt:~4,2%" & set "DD=%dt:~6,2%"
set "HH=%dt:~8,2%" & set "Min=%dt:~10,2%" & set "Sec=%dt:~12,2%"
set "BACKUP_DIR=.github\workflows\backup\%YYYY%%MM%%DD%_%HH%%Min%%Sec%"
mkdir "%BACKUP_DIR%" 2>nul

REM Backup current workflows
echo 💾 Backing up current workflows...
copy .github\workflows\ci.yml "%BACKUP_DIR%\ci-backup.yml" >nul 2>&1
copy .github\workflows\ci-simple.yml "%BACKUP_DIR%\ci-simple-backup.yml" >nul 2>&1
copy .github\workflows\ci-comprehensive.yml "%BACKUP_DIR%\ci-comprehensive-backup.yml" >nul 2>&1
copy .github\workflows\pre-commit.yml "%BACKUP_DIR%\pre-commit-backup.yml" >nul 2>&1

echo ✅ Backed up to: %BACKUP_DIR%

REM Remove old workflows
echo 🗑️ Removing old workflows...
del .github\workflows\ci.yml 2>nul
del .github\workflows\ci-simple.yml 2>nul
del .github\workflows\ci-comprehensive.yml 2>nul
del .github\workflows\pre-commit.yml 2>nul

REM Deploy new workflows
echo 🚀 Deploying streamlined workflows...
move .github\workflows\ci-streamlined.yml .github\workflows\ci.yml
move .github\workflows\pre-commit-streamlined.yml .github\workflows\pre-commit.yml

echo ✅ Streamlined workflows deployed

REM Verify workflows
echo 🔍 Verifying workflows...
if exist ".github\workflows\ci.yml" if exist ".github\workflows\pre-commit.yml" (
    echo ✅ Workflows verified
) else (
    echo ❌ Workflow verification failed
    exit /b 1
)

REM Show current workflow structure
echo.
echo 📋 Current workflow structure:
echo ==============================
dir .github\workflows\ /b

echo.
echo 🎉 CI Migration Complete!
echo =========================
echo.
echo New structure:
echo   - ci.yml: Main CI workflow (quality, tests, build, security, release)
echo   - pre-commit.yml: Pre-commit hook validation
echo   - release.yml: Unchanged (tag-based releases)
echo.
echo Backup location: %BACKUP_DIR%
echo.
echo Next steps:
echo   1. Commit and push changes
echo   2. Create a test PR to verify workflows
echo   3. Check CI execution times and results
echo   4. Remove backup files after confirmation
echo.
echo To rollback if needed:
echo   copy "%BACKUP_DIR%\*" .github\workflows\
