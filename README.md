# scATAC-Seq Analysis Pipeline

End-to-end single-cell ATAC-seq analysis pipeline from fragment files to
motif activity and optional scRNA-seq label transfer. Implements rigorous
multi-metric QC, TF-IDF normalization, LSI dimensionality reduction,
MACS2 peak calling, gene activity scoring, and chromVAR motif deviation
analysis using the JASPAR 2020 TF database.

Designed to resolve cell-type-specific chromatin accessibility landscapes,
identify transcription factor activity differences across populations, and
characterise regulatory programmes underlying distinct cellular states - with
optional cross-modal integration via scRNA-seq label transfer.

This repository is a reproducible, dataset-agnostic pipeline template with
methodology validated against published best practices for single-cell ATAC-seq.

## Repository Structure

```
├── scatac_pipeline.R   # Complete end-to-end pipeline
└── README.md
```

## Pipeline Overview

| Step | Method | Tool |
|------|--------|------|
| Data loading | ChromatinAssay construction | Signac |
| QC | TSS enrichment, nucleosome signal, blacklist, FRiP | Signac |
| Normalization | TF-IDF (term frequency–inverse document frequency) | Signac |
| Dimensionality reduction | SVD / LSI (dim 1 excluded - depth correlated) | Signac |
| Clustering | KNN + Louvain (SLM algorithm) | Seurat |
| Peak calling | Per-cluster MACS2 re-calling | Signac / MACS2 |
| Gene activity | Fragment-based gene body + promoter scoring | Signac |
| Motif analysis | chromVAR deviation scores (JASPAR 2020) | chromVAR |
| DA peaks | Logistic regression with depth covariate | Signac |
| Label transfer | Anchor-based RNA to ATAC transfer (optional) | Seurat |

### Key Design Decisions

**LSI dim 1 exclusion**
The first LSI component almost always correlates strongly with total
sequencing depth per cell rather than biological variation. It is
excluded from UMAP and clustering (`LSI_DIMS_USE` defaults to `2:20`).
The `DepthCor` plot in outputs allows you to verify this empirically.

**Per-cluster MACS2 peak calling**
Cell Ranger ATAC calls peaks on all cells pooled together, which favours
peaks in abundant cell types. Re-calling peaks per cluster with MACS2
recovers cell-type-specific accessible regions and improves sensitivity
for rare populations.

**Logistic regression for DA peaks**
Logistic regression with `nCount_peaks` as a latent variable is the
recommended test for differential accessibility - it directly regresses
out sequencing depth confounding, which is critical for sparse ATAC data.

**chromVAR for TF activity**
chromVAR computes bias-corrected per-cell deviation scores for each
JASPAR motif. This is more informative than simply asking which peaks
are open - it quantifies the activity of specific transcription factors
across clusters, enabling mechanistic interpretation of chromatin states.

**GeneActivity normalization**
Gene activity scores are normalized using `median(nCount_peaks)` as the
scale factor rather than the default 10,000 - appropriate for ATAC sparse
count distributions where library-size assumptions differ from RNA-seq.

**UCSC chromosome style harmonization**
Ensembl-derived annotations use numeric chromosome names (1, 2, …) while
10x fragment files use UCSC-style names (chr1, chr2, …). `seqlevelsStyle`
is set to UCSC explicitly to prevent coordinate mismatch errors in TSS
enrichment and peak calling steps.

## Requirements

```r
# CRAN / GitHub
install.packages(c("Seurat", "Signac", "ggplot2", "patchwork",
                   "dplyr", "future", "Matrix"))

# Bioconductor
BiocManager::install(c(
  "GenomicRanges", "IRanges", "GenomeInfoDb",
  "chromVAR", "motifmatchr", "JASPAR2020", "TFBSTools",
  "BSgenome.Hsapiens.UCSC.hg38",   # or hg19/mm10
  "EnsDb.Hsapiens.v86"              # or EnsDb.Mmusculus.v79 for mouse
))
```

**MACS2** (required for peak calling):
```bash
pip install macs2
```

> **Note:** Requires Seurat v5+ and Signac >= 1.11.

## Usage

**Step 1 - Configure parameters**

Edit `SECTION 0` in `scatac_pipeline.R`:
```r
MATRIX_H5      <- "path/to/filtered_peak_bc_matrix.h5"
FRAGMENT_FILE  <- "path/to/fragments.tsv.gz"
FRAGMENT_INDEX <- "path/to/fragments.tsv.gz.tbi"
METADATA_CSV   <- "path/to/singlecell.csv"
GENOME         <- "hg38"
```

Calibrate QC thresholds after inspecting the VlnPlot and TSSPlot saved
to `scatac_outputs/plots/` before accepting defaults.

**Step 2 - Run pipeline**
```r
source("scatac_pipeline.R")
```

**Optional - Label transfer from scRNA-seq**

If you have a matched scRNA-seq Seurat object, set:
```r
RNASEQ_REF_RDS <- "path/to/seurat_annotated.rds"
RNA_LABEL_COL  <- "your_cell_type_column"
```

## Input Requirements

```
# 10x Cell Ranger ATAC output
filtered_peak_bc_matrix.h5    # peak-barcode count matrix (h5 format)
fragments.tsv.gz               # per-barcode fragment file
fragments.tsv.gz.tbi           # tabix index
singlecell.csv                 # per-barcode metadata
```

## Output Structure

```
scatac_outputs/
├── plots/
│   ├── qc_violin_before_filter.png
│   ├── tss_enrichment_plot.png
│   ├── nucleosome_histogram.png
│   ├── lsi_depth_correlation.png      # verify dim 1 excluded
│   ├── umap_clusters.png
│   ├── umap_macs2_peaks.png
│   ├── gene_activity_canonical_markers.png
│   └── label_transfer_predictions.png
├── tables/
│   ├── qc_summary.csv
│   ├── cluster_composition.csv
│   ├── chromvar_motif_markers.csv
│   └── differential_accessibility_peaks.csv
├── peaks/                             # MACS2 peak files
└── atac_final.rds
```

## References

- Stuart et al. 2021 - Signac (Nature Methods)
- Schep et al. 2017 - chromVAR (Nature Methods)
- Fornes et al. 2020 - JASPAR 2020 (Nucleic Acids Research)
- Zhang et al. 2008 - MACS2 (Genome Biology)
- Stuart et al. 2019 - Seurat anchor-based integration (Cell)
- Amemiya et al. 2019 - ENCODE Blacklist (Scientific Reports)

## Disclaimer

This repository contains reusable analysis code only. It is designed as a
reproducible, dataset-agnostic template applicable to any compatible 10x
Chromium ATAC or multiome dataset via user-supplied parameters.

Dataset-specific findings, result files, and analyses conducted as part of
ongoing or confidential research are not included and remain separate from
this codebase. Parameters require adjustment based on the target tissue
and dataset.

Shared for portfolio and demonstration purposes only. No license is granted
for reuse, modification, or redistribution - please reach out if any discussion
or collaboration.
