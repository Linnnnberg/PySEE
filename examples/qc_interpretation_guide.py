"""
QC Panel Interpretation Guide for Single-Cell RNA-seq Data.

This script explains typical QC ranges and how to interpret QC metrics.
"""

import numpy as np
import pandas as pd
import anndata as ad
from pysee import PySEE
from pysee.panels.qc import QCPanel


def create_qc_interpretation_demo():
    """Create a demo showing different QC scenarios."""
    print("📊 QC Panel Interpretation Guide")
    print("=" * 60)

    # Create different quality scenarios
    scenarios = {
        "High Quality": create_high_quality_data(),
        "Medium Quality": create_medium_quality_data(),
        "Low Quality": create_low_quality_data(),
        "Mixed Quality": create_mixed_quality_data(),
    }

    for scenario_name, adata in scenarios.items():
        print(f"\n🔍 {scenario_name} Dataset:")
        print("-" * 40)
        analyze_qc_metrics(adata, scenario_name)


def create_high_quality_data():
    """Create high-quality single-cell data."""
    n_cells = 1000
    n_genes = 2000

    np.random.seed(42)
    # High expression, low mitochondrial content
    expression_matrix = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(float)

    # Add mitochondrial genes (10% of total)
    n_mito_genes = int(0.1 * n_genes)
    mito_indices = np.random.choice(n_genes, n_mito_genes, replace=False)
    expression_matrix[:, mito_indices] *= 1.5  # Low mitochondrial expression

    # Create gene names
    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
    for i, idx in enumerate(mito_indices):
        gene_names[idx] = f"MT-Gene_{i:03d}"

    # Create cell names
    cell_names = [f"Cell_{i:04d}" for i in range(n_cells)]

    # Calculate QC metrics
    mito_genes = [name.startswith("MT-") for name in gene_names]
    mito_counts = expression_matrix[:, mito_genes].sum(axis=1)
    total_counts = expression_matrix.sum(axis=1)
    mito_percent = mito_counts / total_counts * 100

    # Create metadata
    obs_data = {
        "cell_type": np.random.choice(["T_cell", "B_cell", "NK_cell"], n_cells),
        "total_counts": total_counts,
        "detected_genes": (expression_matrix > 0).sum(axis=1),
        "mito_percent": mito_percent,
    }

    return ad.AnnData(X=expression_matrix, obs=pd.DataFrame(obs_data, index=cell_names), var=pd.DataFrame(index=gene_names))


def create_medium_quality_data():
    """Create medium-quality single-cell data."""
    n_cells = 1000
    n_genes = 2000

    np.random.seed(123)
    # Medium expression, moderate mitochondrial content
    expression_matrix = np.random.negative_binomial(3, 0.4, size=(n_cells, n_genes)).astype(float)

    # Add mitochondrial genes (15% of total)
    n_mito_genes = int(0.15 * n_genes)
    mito_indices = np.random.choice(n_genes, n_mito_genes, replace=False)
    expression_matrix[:, mito_indices] *= 2.5  # Moderate mitochondrial expression

    # Create gene names
    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
    for i, idx in enumerate(mito_indices):
        gene_names[idx] = f"MT-Gene_{i:03d}"

    # Create cell names
    cell_names = [f"Cell_{i:04d}" for i in range(n_cells)]

    # Calculate QC metrics
    mito_genes = [name.startswith("MT-") for name in gene_names]
    mito_counts = expression_matrix[:, mito_genes].sum(axis=1)
    total_counts = expression_matrix.sum(axis=1)
    mito_percent = mito_counts / total_counts * 100

    # Create metadata
    obs_data = {
        "cell_type": np.random.choice(["T_cell", "B_cell", "NK_cell"], n_cells),
        "total_counts": total_counts,
        "detected_genes": (expression_matrix > 0).sum(axis=1),
        "mito_percent": mito_percent,
    }

    return ad.AnnData(X=expression_matrix, obs=pd.DataFrame(obs_data, index=cell_names), var=pd.DataFrame(index=gene_names))


def create_low_quality_data():
    """Create low-quality single-cell data."""
    n_cells = 1000
    n_genes = 2000

    np.random.seed(456)
    # Low expression, high mitochondrial content
    expression_matrix = np.random.negative_binomial(2, 0.5, size=(n_cells, n_genes)).astype(float)

    # Add mitochondrial genes (20% of total)
    n_mito_genes = int(0.2 * n_genes)
    mito_indices = np.random.choice(n_genes, n_mito_genes, replace=False)
    expression_matrix[:, mito_indices] *= 4.0  # High mitochondrial expression

    # Add some cells with very high mitochondrial content (dying cells)
    high_mito_cells = np.random.choice(n_cells, int(0.2 * n_cells), replace=False)
    expression_matrix[high_mito_cells, mito_indices] *= 3.0

    # Create gene names
    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
    for i, idx in enumerate(mito_indices):
        gene_names[idx] = f"MT-Gene_{i:03d}"

    # Create cell names
    cell_names = [f"Cell_{i:04d}" for i in range(n_cells)]

    # Calculate QC metrics
    mito_genes = [name.startswith("MT-") for name in gene_names]
    mito_counts = expression_matrix[:, mito_genes].sum(axis=1)
    total_counts = expression_matrix.sum(axis=1)
    mito_percent = mito_counts / total_counts * 100

    # Create metadata
    obs_data = {
        "cell_type": np.random.choice(["T_cell", "B_cell", "NK_cell"], n_cells),
        "total_counts": total_counts,
        "detected_genes": (expression_matrix > 0).sum(axis=1),
        "mito_percent": mito_percent,
    }

    return ad.AnnData(X=expression_matrix, obs=pd.DataFrame(obs_data, index=cell_names), var=pd.DataFrame(index=gene_names))


def create_mixed_quality_data():
    """Create mixed quality data (realistic scenario)."""
    n_cells = 1000
    n_genes = 2000

    np.random.seed(789)
    # Mixed quality with some good and some bad cells
    expression_matrix = np.random.negative_binomial(4, 0.35, size=(n_cells, n_genes)).astype(float)

    # Add mitochondrial genes (12% of total)
    n_mito_genes = int(0.12 * n_genes)
    mito_indices = np.random.choice(n_genes, n_mito_genes, replace=False)
    expression_matrix[:, mito_indices] *= 2.0

    # Add some low-quality cells
    low_quality_cells = np.random.choice(n_cells, int(0.15 * n_cells), replace=False)
    expression_matrix[low_quality_cells] *= 0.3  # Low expression
    expression_matrix[low_quality_cells, mito_indices] *= 5.0  # High mitochondrial

    # Create gene names
    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
    for i, idx in enumerate(mito_indices):
        gene_names[idx] = f"MT-Gene_{i:03d}"

    # Create cell names
    cell_names = [f"Cell_{i:04d}" for i in range(n_cells)]

    # Calculate QC metrics
    mito_genes = [name.startswith("MT-") for name in gene_names]
    mito_counts = expression_matrix[:, mito_genes].sum(axis=1)
    total_counts = expression_matrix.sum(axis=1)
    mito_percent = mito_counts / total_counts * 100

    # Create metadata
    obs_data = {
        "cell_type": np.random.choice(["T_cell", "B_cell", "NK_cell"], n_cells),
        "total_counts": total_counts,
        "detected_genes": (expression_matrix > 0).sum(axis=1),
        "mito_percent": mito_percent,
    }

    return ad.AnnData(X=expression_matrix, obs=pd.DataFrame(obs_data, index=cell_names), var=pd.DataFrame(index=gene_names))


def analyze_qc_metrics(adata, scenario_name):
    """Analyze and interpret QC metrics."""
    # Calculate metrics
    mito_genes = adata.var_names.str.startswith(("MT-", "mt-", "Mt-"))
    mito_counts = adata[:, mito_genes].X.sum(axis=1)
    total_counts = adata.X.sum(axis=1)
    mito_percent = mito_counts / total_counts * 100
    detected_genes = (adata.X > 0).sum(axis=1)

    print(f"Dataset: {adata.n_obs} cells, {adata.n_vars} genes")
    print()

    # Mitochondrial percentage analysis
    print("🧬 Mitochondrial Gene Percentage:")
    print(f"   Mean: {mito_percent.mean():.1f}%")
    print(f"   Median: {np.median(mito_percent):.1f}%")
    print(f"   Range: {mito_percent.min():.1f}% - {mito_percent.max():.1f}%")

    # Interpret mitochondrial percentage
    mito_interpretation = interpret_mito_percent(mito_percent)
    print(f"   Interpretation: {mito_interpretation}")
    print()

    # Total counts analysis
    print("📈 Total Gene Counts:")
    print(f"   Mean: {total_counts.mean():.0f}")
    print(f"   Median: {np.median(total_counts):.0f}")
    print(f"   Range: {total_counts.min():.0f} - {total_counts.max():.0f}")

    # Interpret total counts
    counts_interpretation = interpret_total_counts(total_counts)
    print(f"   Interpretation: {counts_interpretation}")
    print()

    # Detected genes analysis
    print("🔬 Detected Genes per Cell:")
    print(f"   Mean: {detected_genes.mean():.0f}")
    print(f"   Median: {np.median(detected_genes):.0f}")
    print(f"   Range: {detected_genes.min():.0f} - {detected_genes.max():.0f}")

    # Interpret detected genes
    genes_interpretation = interpret_detected_genes(detected_genes)
    print(f"   Interpretation: {genes_interpretation}")
    print()

    # Overall quality assessment
    overall_quality = assess_overall_quality(mito_percent, total_counts, detected_genes)
    print(f"🎯 Overall Quality Assessment: {overall_quality}")
    print()


def interpret_mito_percent(mito_percent):
    """Interpret mitochondrial percentage."""
    mean_mito = mito_percent.mean()
    high_mito_cells = (mito_percent > 20).sum()
    high_mito_pct = (mito_percent > 20).mean() * 100

    if mean_mito < 10:
        return "✅ Excellent - Low mitochondrial content"
    elif mean_mito < 15:
        return "✅ Good - Acceptable mitochondrial content"
    elif mean_mito < 20:
        return "⚠️ Moderate - Some cells may need filtering"
    elif mean_mito < 30:
        return "❌ Poor - High mitochondrial content, consider filtering"
    else:
        return "❌ Very Poor - Very high mitochondrial content, extensive filtering needed"


def interpret_total_counts(total_counts):
    """Interpret total gene counts."""
    mean_counts = total_counts.mean()
    low_counts = (total_counts < 1000).sum()
    high_counts = (total_counts > 50000).sum()

    if mean_counts > 10000:
        return "✅ Excellent - High expression levels"
    elif mean_counts > 5000:
        return "✅ Good - Adequate expression levels"
    elif mean_counts > 2000:
        return "⚠️ Moderate - Low expression, may need filtering"
    else:
        return "❌ Poor - Very low expression, extensive filtering needed"


def interpret_detected_genes(detected_genes):
    """Interpret number of detected genes."""
    mean_genes = detected_genes.mean()
    low_genes = (detected_genes < 500).sum()

    if mean_genes > 2000:
        return "✅ Excellent - High gene detection"
    elif mean_genes > 1000:
        return "✅ Good - Adequate gene detection"
    elif mean_genes > 500:
        return "⚠️ Moderate - Low gene detection, may need filtering"
    else:
        return "❌ Poor - Very low gene detection, extensive filtering needed"


def assess_overall_quality(mito_percent, total_counts, detected_genes):
    """Assess overall data quality."""
    # Scoring system
    score = 0

    # Mitochondrial score (0-3)
    mean_mito = mito_percent.mean()
    if mean_mito < 10:
        score += 3
    elif mean_mito < 15:
        score += 2
    elif mean_mito < 20:
        score += 1

    # Counts score (0-3)
    mean_counts = total_counts.mean()
    if mean_counts > 10000:
        score += 3
    elif mean_counts > 5000:
        score += 2
    elif mean_counts > 2000:
        score += 1

    # Genes score (0-3)
    mean_genes = detected_genes.mean()
    if mean_genes > 2000:
        score += 3
    elif mean_genes > 1000:
        score += 2
    elif mean_genes > 500:
        score += 1

    # Overall assessment
    if score >= 8:
        return "🟢 High Quality - Ready for analysis"
    elif score >= 6:
        return "🟡 Medium Quality - Some filtering recommended"
    elif score >= 4:
        return "🟠 Low Quality - Extensive filtering needed"
    else:
        return "🔴 Very Low Quality - Consider re-sequencing"


def show_qc_guidelines():
    """Show general QC guidelines."""
    print("\n📋 QC Guidelines for Single-Cell RNA-seq Data")
    print("=" * 60)

    print("\n🧬 Mitochondrial Gene Percentage:")
    print("   ✅ Excellent: < 10%")
    print("   ✅ Good: 10-15%")
    print("   ⚠️ Acceptable: 15-20%")
    print("   ❌ Poor: 20-30%")
    print("   ❌ Very Poor: > 30%")
    print("   📝 Note: High mitochondrial % indicates dying cells")

    print("\n📈 Total Gene Counts (UMI counts):")
    print("   ✅ Excellent: > 10,000")
    print("   ✅ Good: 5,000-10,000")
    print("   ⚠️ Acceptable: 2,000-5,000")
    print("   ❌ Poor: 1,000-2,000")
    print("   ❌ Very Poor: < 1,000")
    print("   📝 Note: Low counts indicate poor cell capture")

    print("\n🔬 Detected Genes per Cell:")
    print("   ✅ Excellent: > 2,000 genes")
    print("   ✅ Good: 1,000-2,000 genes")
    print("   ⚠️ Acceptable: 500-1,000 genes")
    print("   ❌ Poor: 200-500 genes")
    print("   ❌ Very Poor: < 200 genes")
    print("   📝 Note: Low gene detection indicates poor RNA quality")

    print("\n🎯 Recommended Filtering Thresholds:")
    print("   • Mitochondrial %: < 20% (or < 15% for high-quality data)")
    print("   • Total counts: 1,000 - 50,000")
    print("   • Detected genes: 200 - 5,000")
    print("   • Adjust based on your specific dataset and cell type")

    print("\n💡 Tips for QC Interpretation:")
    print("   • Compare to published studies of similar cell types")
    print("   • Consider your experimental conditions")
    print("   • Look for batch effects in QC metrics")
    print("   • Use multiple metrics together, not just one")
    print("   • Start with lenient thresholds and tighten gradually")


def main():
    """Run the QC interpretation guide."""
    create_qc_interpretation_demo()
    show_qc_guidelines()

    print("\n🎉 QC Interpretation Guide Complete!")
    print("   Use these guidelines to interpret your QC panel results.")


if __name__ == "__main__":
    main()
