# Gallbladder Single-Cell Analysis Scripts

This repository stores R scripts for gallbladder cancer single-cell RNA-seq analysis.

The analysis is designed to run on the project Linux server. Raw 10x Genomics matrices and large generated results are intentionally excluded from GitHub.

## Repository Layout

```text
.
|-- README.md
|-- .gitignore
`-- scripts/
    |-- 01_tissue2_basic_seurat.R
    |-- 02_merged_basic_analysis.R
    `-- 03_tissue2_stk31_nk_analysis.R
```

## Data Policy

Do not commit raw sequencing data or large analysis artifacts.

Excluded examples:

- raw 10x data directories
- `gallbladder_cancer/`
- `results/**/*.rds`
- large PDF, PNG, HTML, and log outputs

The scripts assume raw 10x data is available on the analysis server under:

```text
/home/zhuweiyu/codex-r/gallbladder_cancer
```

## Script Order

1. `scripts/01_tissue2_basic_seurat.R`
   - Runs the basic Seurat workflow for tissue2.
   - Performs QC, normalization, PCA, clustering, UMAP, and saves the tissue2 Seurat object.

2. `scripts/02_merged_basic_analysis.R`
   - Runs the merged workflow across tissue1, tissue2, tissue4, and tissue5.
   - Includes the merged STK31/NK analysis, the former merged follow-up checks, CellChat/GO plotting, and focused STK31-high epithelial/NK CellChat figures.
   - Writes the core merged outputs to `/home/zhuweiyu/codex-r/results/merged_basic_seurat` and STK31/NK outputs under `/home/zhuweiyu/codex-r/results/merged_stk31_nk_analysis`.

3. `scripts/03_tissue2_stk31_nk_analysis.R`
   - Runs the tissue2 STK31/NK analysis from the tissue2 basic Seurat object.
   - Includes the former validation checks: cell identity markers, NK signatures, pathway panels, top DE exports, and optional GO analysis.
   - Writes the main outputs to `/home/zhuweiyu/codex-r/results/tissue2_stk31_nk_analysis` and validation outputs under its `validation/` subdirectory.

## Expected Runtime

These scripts are intended to run with R and Seurat on Linux. In this project, the expected server paths are:

```text
Working directory: /home/zhuweiyu/codex-r
Rscript: /usr/bin/Rscript
Scripts: /home/zhuweiyu/codex-r/scripts
Results: /home/zhuweiyu/codex-r/results
```

Example run pattern:

```bash
cd /home/zhuweiyu/codex-r
/usr/bin/Rscript scripts/01_tissue2_basic_seurat.R 2>&1 | tee logs/01_tissue2_basic_seurat.log
```

## Notes

The repository is for source scripts and lightweight project documentation. Reproducibility depends on the external data directory and the R package environment on the analysis server.
