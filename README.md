# Gallbladder RNA-seq / single-cell analysis scripts

This repository stores R scripts for gallbladder cancer single-cell RNA-seq analysis.

The analysis is designed to run on the remote Linux server used by this project, with input 10x Genomics matrices kept outside GitHub. Raw data and large generated results are intentionally excluded from the repository.

## Repository layout

```text
.
├── README.md
├── .gitignore
└── scripts/
    ├── 01_tissue2_basic_seurat.R
    ├── 02_merged_basic_analysis.R
    ├── 03_tissue2_stk31_nk_analysis.R
    └── 04_tissue2_stk31_nk_validation.R
```

## Data policy

Do not commit raw sequencing data or large analysis artifacts.

Excluded examples:

- `胆囊癌/`
- `gallbladder_cancer/`
- `results/**/*.rds`
- large PDF, PNG, HTML, and log outputs

The scripts assume the raw 10x data is available on the analysis server under:

```text
/home/zhuweiyu/codex-r/gallbladder_cancer
```

## Script order

1. `scripts/01_tissue2_basic_seurat.R`
   - Runs a basic Seurat workflow for tissue2.
   - Performs QC, normalization, PCA, clustering, UMAP, and saves a Seurat object.

2. `scripts/02_merged_basic_analysis.R`
   - Runs a merged/basic analysis workflow across available samples.

3. `scripts/03_tissue2_stk31_nk_analysis.R`
   - Focuses on STK31-high clusters and NK candidate populations.
   - Produces marker, annotation, UMAP, ligand-receptor, and pathway evidence outputs.

4. `scripts/04_tissue2_stk31_nk_validation.R`
   - Validates STK31/NK findings with cell identity markers, NK signatures, pathway panels, and optional GO analysis.

## Expected runtime

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
