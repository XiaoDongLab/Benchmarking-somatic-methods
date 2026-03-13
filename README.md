# Sequencing Quality and Cross-Test Somatic Benchmark Pipelines

Pipelines for sequencing quality assessment and cross-test somatic
variant benchmarking.

This repository contains **two independent analysis workflows**:

1.  **Sequencing Metrics Pipeline** -- calculates sequencing
    quality metrics and generates coverage, Lorenz, depth, and VAF plots
    from BAM and VCF files.
2.  **Cross‑Test Somatic Mutation Benchmark Pipeline** -- evaluates
    somatic mutation calls using cross‑sample comparison and generates
    benchmarking statistics and plots.

------------------------------------------------------------------------

# Table of Contents

-   [Repository Structure](#repository-structure)
-   [Pipeline 1: Sequencing Metrics](#pipeline-1-sequencing-metrics)
    -   [Inputs](#inputs)
    -   [Run](#run)
    -   [Outputs](#outputs)
    -   [Software Requirements](#software-requirements)
    -   [Conda Environment](#conda-environment)
-   [Pipeline 2: Cross‑Test Somatic
    Benchmark](#pipeline-2-cross-test-somatic-benchmark)
    -   [Inputs](#inputs-1)
    -   [Run](#run-1)
    -   [Outputs](#outputs-1)
-   [License](#license)

------------------------------------------------------------------------

# Repository Structure

Example layout:

    repo/
    │
    ├── run_sequencing_metrics_pipeline.sh
    ├── plot_sequencing_metrics.R
    │
    ├── run_sm.sh
    ├── plot_sm.R
    │
    ├── sequencing_metrics.yml
    │
    └── README.md

The repository provides **two separate pipelines**. They can be used
independently.

------------------------------------------------------------------------

# Pipeline 1: Sequencing Metrics

This pipeline computes sequencing statistics and generates plots used to
evaluate sequencing data quality.

Metrics include:

-   reads vs covered bases
-   Lorenz curve and Gini coefficient
-   depth distribution
-   chromosome coverage windows
-   variant allele frequency (VAF)

The pipeline processes **multiple samples automatically**.

------------------------------------------------------------------------

# Inputs

Create a text file listing sample IDs:

    SampleA
    SampleB
    SampleC

Each sample must have a directory:

    ../lab_data/<sample_id>/

Each directory must contain:

    <sample_id>_WGNS_filtered.bam
    ssnv.results.muts.vcf.gz
    sindel.results.indel.vcf.gz

Reference genome:

    reference.fa
    reference.fa.fai

------------------------------------------------------------------------

# Run

Basic usage:

    bash run_sequencing_metrics_pipeline.sh sample_ids.txt

Full usage:

    bash run_sequencing_metrics_pipeline.sh \
        sample_ids.txt \
        8 \
        ../lab_data \
        ./sequencing_metrics_results \
        reference.fa

------------------------------------------------------------------------

# Outputs

## Per‑sample outputs

Inside each sample directory:

    cov_<sample>.txt
    <sample>_lorenz.txt
    <sample>_coverage.txt
    <sample>_distribution.txt
    chrom_cov_table.txt

    pass.ssnv.vcf.gz
    pass.sindel.vcf.gz

    <sample>_ssnv_AD.vcf
    <sample>_sindel_AD.vcf

## Aggregated tables

    sequencing_metrics_results/bridge_tsv/

Example files:

    depth_curve.tsv
    lorenz_curve.tsv
    depth_distribution.tsv
    chrom_coverage_1Mb.tsv
    vaf_values.tsv

## Final plots

    sequencing_metrics_results/plots/

Examples:

    reads_coverage_curve.pdf
    lorenz_curve_and_gini.pdf
    depth_distribution_0_100.pdf
    chrom_window_coverage_all_samples.pdf
    vaf_raw_including_zero.pdf
    vaf_raw_excluding_zero.pdf

------------------------------------------------------------------------

# Software Requirements

Command‑line tools:

    samtools
    bcftools
    bam-lorenz-coverage
    Rscript
    awk
    sed
    bc

Required R packages:

    ggplot2
    dplyr
    tidyr
    readr
    patchwork
    stringr
    ggbeeswarm
    scales

------------------------------------------------------------------------

# Conda Environment

Create environment:

    conda env create -f sequencing_metrics.yml
    conda activate sequencing_metrics

------------------------------------------------------------------------

# Pipeline 2: Cross‑Test Somatic Benchmark

This pipeline evaluates somatic mutation calls by comparing two mutation
call sets and estimating performance metrics.

Rationale:

This analysis uses cross-test comparisons to separate signal that is
reproducible across call sets from signal that is likely technical
noise. By benchmarking overlap with germline-derived expectations and a
pseudo-false-positive reference, it provides a practical way to estimate
precision, pseudo-true-positive rates, and burden correction factors for
somatic call sets.

The pipeline produces:

-   true positive vs false positive comparisons
-   mutation burden estimation
-   mixture inference plots

------------------------------------------------------------------------

# Inputs

Required VCF files:

    somatic.vcf.gz
    cross_test_somatic.vcf.gz
    germline_sample.vcf.gz
    germline_test.vcf.gz
    pseudo_fp_reference.vcf.gz

Optional:

    reference.fa

All VCF files should be **compressed and indexed**.

------------------------------------------------------------------------

# Run

Minimal example:

    bash run_sm.sh \
      --sample SAMPLE \
      --callable 2500000000 \
      --somatic som.vcf.gz \
      --cross cross.vcf.gz \
      --germS gS.vcf.gz \
      --germT gT.vcf.gz \
      --fpref ref.vcf.gz \
      --out results

Optional arguments:

    --filter PASS
    --ref reference.fa
    --reflen 3137300923

------------------------------------------------------------------------

# Outputs

Directory structure:

    results/
    ├── plots/
    ├── tables/
    ├── tmp/
    └── vcf/

## Tables

    <sample>.x.tsv
    <sample>.m.tsv
    <sample>.snv_xy.tsv
    <sample>.ind_xy.tsv
    <sample>.mix.tsv
    <sample>.burden.tsv

## Plots

    <sample>.snv_tp_vs_fp.pdf
    <sample>.ind_tp_vs_fp.pdf
    <sample>.mix_inferred.pdf
    <sample>.burden_snv.pdf
    <sample>.burden_ind.pdf

Temporary `bcftools isec` results are stored in:

    results/tmp/

------------------------------------------------------------------------


# License

For academic research use.
