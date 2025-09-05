# PySEE Version Strategy

## Overview

PySEE follows [Semantic Versioning (SemVer)](https://semver.org/) with a clear release strategy for scientific software development.

## Version Format: `MAJOR.MINOR.PATCH`

- **MAJOR**: Breaking changes that require user code modifications
- **MINOR**: New features, backward compatible
- **PATCH**: Bug fixes, backward compatible

## Current Version Status

### **v0.1.0** - MVP Release ✅ COMPLETED
- **Status**: Production-ready MVP
- **Features**: Core architecture, UMAP panel, Violin panel, linking system
- **Stability**: Stable for basic bioinformatics workflows

### **v0.1.1** - Infrastructure Optimization ✅ COMPLETED
- **Status**: Infrastructure improvements
- **Features**: CI/CD optimization, Git workflow, documentation
- **Stability**: Enhanced development process

### **v0.1.2** - Ready for Release NOW! 🚀 READY
- **Status**: Ready for immediate release
- **Features**: MVP + Heatmap panel (completed)
- **Target**: Immediate release to PyPI

### **v0.2.0** - Next Major Feature Release 🔄 IN DEVELOPMENT
- **Status**: Planned for next development phase
- **Features**: QC metrics, Dot plot panel, advanced selection
- **Target**: Q1 2025

## Release Strategy

### Development Phases

#### **v0.x Series** - Alpha/Beta Development
- **Purpose**: Rapid iteration and feature development
- **Stability**: API may change between minor versions
- **Users**: Early adopters, developers, researchers
- **Release Frequency**: Every 2-4 weeks

#### **v1.x Series** - Stable Production
- **Purpose**: Production-ready, stable API
- **Stability**: Backward compatibility maintained
- **Users**: Production users, research labs, institutions
- **Release Frequency**: Every 2-3 months

#### **v2.x Series** - Future Evolution
- **Purpose**: Major architectural improvements
- **Stability**: May include breaking changes
- **Users**: Advanced users, new adopters
- **Release Frequency**: Every 6-12 months

### Version Numbering Rules

#### **Patch Versions (0.1.x)**
- Bug fixes
- Documentation updates
- Performance improvements
- Security patches
- **No new features**

#### **Minor Versions (0.x.0)**
- New visualization panels
- New features
- API enhancements
- **Backward compatible**

#### **Major Versions (x.0.0)**
- Breaking API changes
- Major architectural changes
- Removal of deprecated features
- **Requires user code updates**

## Release Process

### 1. Development Branch Strategy
```
main (stable releases)
├── develop (integration branch)
    ├── feature/heatmap-panel
    ├── feature/qc-metrics-panel
    └── feature/dotplot-panel
```

### 2. Release Workflow

#### **Patch Release (0.1.2)**
```bash
# 1. Create hotfix branch from main
git checkout main
git checkout -b hotfix/bug-fix-description

# 2. Fix the issue
# ... make changes ...

# 3. Update version in setup.py and pyproject.toml
# 4. Update CHANGELOG.md
# 5. Create PR: hotfix/bug-fix-description → main
# 6. After merge, tag release
git tag -a v0.1.2 -m "Release version 0.1.2"
git push origin v0.1.2
```

#### **Minor Release (0.2.0)**
```bash
# 1. Create release branch from develop
git checkout develop
git checkout -b release/v0.2.0

# 2. Update version numbers
# 3. Update CHANGELOG.md with all new features
# 4. Final testing and documentation updates
# 5. Create PR: release/v0.2.0 → main
# 6. After merge, tag release
git tag -a v0.2.0 -m "Release version 0.2.0"
git push origin v0.2.0

# 7. Merge back to develop
git checkout develop
git merge main
git push origin develop
```

### 3. Version Update Locations

When releasing, update version in:
- `setup.py` - `version="0.2.0"`
- `pyproject.toml` - `version = "0.2.0"`
- `CHANGELOG.md` - Add new version section
- `README.md` - Update version references if needed

## Version Compatibility

### Python Version Support

#### **v0.1.x - v0.2.x**
- **Python**: 3.9, 3.10, 3.11, 3.12
- **Rationale**: Modern scientific Python ecosystem

#### **v1.x+**
- **Python**: 3.10, 3.11, 3.12 (drop 3.9)
- **Rationale**: Focus on actively supported Python versions

### Dependency Compatibility

#### **Core Dependencies**
- **AnnData**: >=0.8.0 (stable API)
- **Scanpy**: >=1.9.0 (modern features)
- **Plotly**: >=5.0.0 (interactive features)
- **NumPy**: >=1.21.0 (scientific computing)

#### **Version Pinning Strategy**
- **Major versions**: Pin to avoid breaking changes
- **Minor versions**: Allow updates for bug fixes
- **Patch versions**: Always allow latest

## Release Schedule

### **2025 Roadmap**

#### **Q1 2025**
- **v0.2.0**: Heatmap, QC Metrics, Dot Plot panels
- **v0.2.1**: Advanced selection tools
- **v0.2.2**: Jupyter widget integration

#### **Q2 2025**
- **v0.3.0**: Genome browser integration
- **v0.3.1**: Spatial transcriptomics viewer
- **v0.3.2**: Plugin system foundation

#### **Q3 2025**
- **v1.0.0**: First stable release
- **v1.0.1**: Production bug fixes
- **v1.1.0**: Web deployment capabilities

#### **Q4 2025**
- **v1.2.0**: Cloud-scale data support
- **v1.3.0**: Advanced analytics integration

## Version Communication

### **Release Announcements**
- **GitHub Releases**: Detailed changelog and download links
- **PyPI**: Package updates with metadata
- **Documentation**: Version-specific documentation
- **Community**: GitHub Discussions, Twitter, scientific forums

### **Deprecation Policy**
- **Notice Period**: 2 minor versions before removal
- **Documentation**: Clear migration guides
- **Warnings**: Runtime warnings for deprecated features
- **Timeline**: 6-month minimum deprecation period

## Quality Gates

### **Release Criteria**

#### **Patch Release (0.1.x)**
- ✅ All tests pass
- ✅ No new linting errors
- ✅ Documentation updated
- ✅ Backward compatibility maintained

#### **Minor Release (0.x.0)**
- ✅ All tests pass
- ✅ New features documented
- ✅ Performance benchmarks met
- ✅ User feedback incorporated
- ✅ Migration guide provided (if needed)

#### **Major Release (x.0.0)**
- ✅ All tests pass
- ✅ Breaking changes documented
- ✅ Migration guide provided
- ✅ Community feedback incorporated
- ✅ Long-term support plan

## Version Strategy Benefits

### **For Users**
- **Predictable**: Clear version meaning
- **Stable**: Backward compatibility in minor versions
- **Transparent**: Clear deprecation and migration paths

### **For Developers**
- **Flexible**: Rapid iteration in 0.x series
- **Structured**: Clear release process
- **Collaborative**: Community-driven feature development

### **For Scientific Community**
- **Reliable**: Stable APIs for research workflows
- **Evolutive**: Continuous improvement and new features
- **Reproducible**: Version-specific documentation and examples

## Release Readiness Assessment

### **✅ READY FOR PYPI RELEASE NOW!**

#### **v0.1.2 Release Package Contents:**
- **Core Architecture**: AnnData integration, dashboard engine
- **Visualization Panels**: UMAP, Violin, Heatmap (3 panels)
- **Panel Linking**: Selection propagation and interaction system
- **Code Export**: Reproducible Python code generation
- **Jupyter Integration**: Seamless notebook experience
- **Professional Quality**: Comprehensive testing, CI/CD, documentation

#### **Release Benefits:**
- **Early User Feedback**: Get community input on core features
- **Package Presence**: Establish PySEE in the Python ecosystem
- **Research Impact**: Enable immediate use in bioinformatics workflows
- **Community Building**: Attract contributors and users

#### **Release Process:**
1. Complete QC Metrics Panel (current task)
2. Update version to v0.1.2 in setup.py and pyproject.toml
3. Update CHANGELOG.md with new features
4. Create GitHub release with automated PyPI upload
5. Announce to bioinformatics community

## Current Status

- **Active Version**: v0.1.1
- **Next Release**: v0.1.2 (Ready for immediate release!)
- **Development Branch**: `develop`
- **Release Process**: GitHub Actions automated
- **Documentation**: Version-specific guides

---

*This version strategy ensures PySEE evolves systematically while maintaining stability for scientific research workflows.*
