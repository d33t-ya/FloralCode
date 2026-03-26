# FloralCode 🌸

A bioinformatics pipeline for differential gene expression analysis in plant stress response, built to process and visualize large-scale RNA-seq data across rose and rice samples.

## What It Does

FloralCode automates the full workflow from raw gene expression data to interpretable visualizations — cleaning, normalizing, and statistically analyzing 80,000+ gene expressions across 8 samples to identify which genes are differentially expressed under stress conditions.

## Data Sources

Real expression data obtained from publicly available biological databases including:
- [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/)
- [GSEA / MSigDB](https://www.gsea-msigdb.org/)

Organisms studied: *Rosa* (rose) and *Oryza sativa* (rice)

## Pipeline Overview

```
Raw Expression Data
        ↓
  Data Cleaning & Normalization
        ↓
  Differential Expression Analysis
  (log2 fold change, t-tests, FDR correction)
        ↓
  Visualization
  (Volcano Plots, PCA Plots, Heatmaps)
        ↓
  Stress-Response Pattern Identification
```

## Features

- **Data processing** — cleans and normalizes large biological datasets for consistent downstream analysis
- **Statistical analysis** — applies log2 fold change calculations, independent t-tests, and FDR (False Discovery Rate) correction to identify significantly differentially expressed genes
- **Visualization** — auto-generates volcano plots, PCA plots, and heatmaps to reveal expression patterns
- **Synthetic data generator** — built-in benchmarking tool to validate pipeline accuracy and ensure reproducibility
- **Scalable design** — architected to handle large datasets (80,000+ gene expressions) without manual intervention at each step

## Tech Stack

| Tool | Purpose |
|---|---|
| Python | Core pipeline logic |
| Kallisto | RNA-seq quantification |
| Matplotlib | Plot generation |
| Seaborn | Heatmap and statistical visualizations |
| Jupyter / Anaconda | Development environment |

## Output Visualizations

> *Figures folder in the outputfolder of the code*

Sample outputs include:
- **Volcano plots** — highlights statistically significant differentially expressed genes
- **PCA plots** — shows variance structure across the 8 samples
- **Heatmaps** — clusters genes by expression pattern across conditions

## Context

Built as a course project for a Bioinformatics class at the University of Washington Bothell. Given an open-ended prompt to create anything in the field, this pipeline was designed from scratch as a personal exploration of plant stress genomics and scalable data workflows.
