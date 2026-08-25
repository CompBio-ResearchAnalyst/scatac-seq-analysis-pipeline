# ==============================================================================
# scATAC-Seq Analysis Pipeline - Signac
# Author: Rahul V Sangoji
#
# Description:
#   End-to-end single-cell ATAC-seq analysis pipeline using Signac and Seurat.
#   Covers fragment file loading, chromatin assay construction, QC
#   (TSS enrichment, nucleosome signal, blacklist filtering), TF-IDF
#   normalization, LSI dimensionality reduction, UMAP/clustering, MACS2
#   peak calling, gene activity scoring, motif analysis via chromVAR,
#   and optional label transfer from a matched scRNA-seq object.
#
#   Designed as a dataset-agnostic template for any 10x Chromium ATAC
#   or multiome dataset. Supply parameters in SECTION 0.
#
# Requires: Signac >= 1.11, Seurat v5+, chromVAR, MACS2 (system install)
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION 0: CONFIGURATION - Edit these before running
# ------------------------------------------------------------------------------

# --- Paths ---
# MATRIX_H5: path to filtered_peak_bc_matrix.h5 from Cell Ranger ATAC output
# This is the count matrix file - NOT peaks.bed
MATRIX_H5       <- "YOUR_MATRIX_H5_PATH"         # e.g. "outs/filtered_peak_bc_matrix.h5"
FRAGMENT_FILE   <- "YOUR_FRAGMENTS_FILE_PATH"     # e.g. "outs/fragments.tsv.gz"
FRAGMENT_INDEX  <- "YOUR_FRAGMENTS_INDEX_PATH"    # e.g. "outs/fragments.tsv.gz.tbi"
METADATA_CSV    <- "YOUR_METADATA_CSV_PATH"       # e.g. "outs/singlecell.csv"
OUTPUT_DIR      <- "scatac_outputs"               # all outputs written here

# --- Genome ---
# Supported: "hg38", "hg19", "mm10"
# Install matching BSgenome + EnsDb packages (see Requirements in README)
GENOME          <- "hg38"

# --- QC thresholds - calibrate to your VlnPlot/TSSplot inspection ---
MIN_COUNTS      <- 1000    # minimum total ATAC fragments per cell - example
MAX_COUNTS      <- 50000   # maximum total ATAC fragments (likely doublets above) - example
MIN_TSS         <- 2       # TSS enrichment score cutoff - standard minimum; aim for > 3
MAX_NUCL_SIGNAL <- 4       # nucleosome signal cutoff - lower = better fragment quality
PEAK_PCT_MIN    <- 15      # minimum % reads in peaks per cell

# --- Dimensionality reduction ---
N_LSI_DIMS      <- 30      # total LSI dimensions to compute
LSI_DIMS_USE    <- 2:20    # dims to use for UMAP/clustering - EXCLUDE dim 1
                           # (dim 1 almost always correlates with sequencing depth)
CLUSTER_RES     <- 0.5     # Louvain resolution

# --- Peak calling ---
MACS2_GENOME    <- "hs"    # MACS2 genome size flag: "hs" (human) or "mm" (mouse)

# --- Optional: scRNA-seq label transfer ---
# Set to path of a saved Seurat scRNA-seq object (.rds) if available
# Leave as NULL to skip label transfer
RNASEQ_REF_RDS  <- NULL    # e.g. "seurat_annotated.rds"
RNA_LABEL_COL   <- "cell_type_major_revised"  # metadata column with cell type labels

set.seed(42)

# ------------------------------------------------------------------------------
# SECTION 1: PACKAGE SETUP
# ------------------------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c("Seurat", "Signac", "ggplot2", "patchwork",
               "dplyr", "future", "Matrix")
bioc_pkgs <- c("GenomicRanges", "IRanges", "GenomeInfoDb",
               "chromVAR", "motifmatchr", "JASPAR2020",
               "TFBSTools", "BSgenome.Hsapiens.UCSC.hg38",
               "EnsDb.Hsapiens.v86")

for (pkg in cran_pkgs) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
}
for (pkg in bioc_pkgs) {
  if (!require(pkg, character.only = TRUE))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
}

library(Signac);    library(Seurat);   library(ggplot2)
library(patchwork); library(dplyr);    library(Matrix)
library(GenomicRanges); library(IRanges)
library(chromVAR);  library(motifmatchr)
library(JASPAR2020); library(TFBSTools)
library(future);    plan("sequential")

options(future.globals.maxSize = 8000 * 1024^2)

dir.create(OUTPUT_DIR,                             showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "plots"),         showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "tables"),        showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "peaks"),         showWarnings = FALSE, recursive = TRUE)

# Select genome-specific objects
# Update these if using hg19 or mm10
genome_bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
genome_ensdb    <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
blacklist_gr    <- Signac::blacklist_hg38_unified  # built-in ENCODE blacklist

cat("Genome configured:", GENOME, "\n")

# ------------------------------------------------------------------------------
# SECTION 2: DATA LOADING - ChromatinAssay Construction
# ------------------------------------------------------------------------------
cat("\n=== SECTION 2: DATA LOADING ===\n")

# Load peak-barcode count matrix from Cell Ranger ATAC output
# MATRIX_H5 is the filtered_peak_bc_matrix.h5 file - the count matrix in h5 format
# Row names will be peak coordinates (e.g. "chr1:10000-10500")
counts <- Read10X_h5(MATRIX_H5)

# Read per-barcode metadata (singlecell.csv from Cell Ranger ATAC)
# Contains: peak_region_fragments, passed_filters, etc.
metadata <- read.csv(METADATA_CSV, row.names = 1)

# Align metadata to cells in count matrix
shared_cells <- intersect(colnames(counts), rownames(metadata))
counts       <- counts[, shared_cells]
metadata     <- metadata[shared_cells, ]

# Build ChromatinAssay with fragment file for downstream operations
# (TSS enrichment, nucleosome signal, peak calling all require fragment access)
chrom_assay <- CreateChromatinAssay(
  counts       = counts,
  sep          = c(":", "-"),
  fragments    = FRAGMENT_FILE,
  min.cells    = 10,
  min.features = 200
)

atac_obj <- CreateSeuratObject(
  counts    = chrom_assay,
  assay     = "peaks",
  meta.data = metadata
)

cat("Cells loaded:", ncol(atac_obj), "\n")
cat("Peaks loaded:", nrow(atac_obj), "\n")

# Add genome annotations - required for TSS enrichment calculation
annotations <- GetGRangesFromEnsDb(ensdb = genome_ensdb)
seqlevelsStyle(annotations) <- "UCSC"   # ensures chr prefix compatibility
Annotation(atac_obj) <- annotations

# ------------------------------------------------------------------------------
# SECTION 3: QUALITY CONTROL
# ------------------------------------------------------------------------------
cat("\n=== SECTION 3: QC ===\n")

# TSS enrichment: signal enrichment at transcription start sites - proxy for
# open chromatin quality. Higher = better. Cutoff typically >= 2, ideally >= 3.
atac_obj <- TSSEnrichment(atac_obj, fast = FALSE)

# Nucleosome signal: ratio of mono-nucleosomal to sub-nucleosomal fragments.
# Low signal = good library quality (sub-nucleosomal fragments dominate).
atac_obj <- NucleosomeSignal(atac_obj)

# Fraction of reads in peaks - calculated from Cell Ranger singlecell.csv metadata
atac_obj$pct_reads_in_peaks <- atac_obj$peak_region_fragments /
  atac_obj$passed_filters * 100

# Blacklist ratio - reads mapping to ENCODE blacklist signal artefacts
atac_obj$blacklist_ratio <- FractionCountsInRegion(
  object  = atac_obj,
  assay   = "peaks",
  regions = blacklist_gr
)

# Visualise QC distributions before filtering - inspect before setting cutoffs
p_qc <- VlnPlot(
  object   = atac_obj,
  features = c("nCount_peaks", "TSS.enrichment",
               "nucleosome_signal", "pct_reads_in_peaks",
               "blacklist_ratio"),
  pt.size  = 0,
  ncol     = 3
)
ggsave(file.path(OUTPUT_DIR, "plots", "qc_violin_before_filter.png"),
       p_qc, width = 18, height = 10, dpi = 300)

# TSS enrichment plot - bimodal distribution indicates mixed quality
p_tss <- TSSPlot(atac_obj, group.by = "high.tss") + NoLegend()
ggsave(file.path(OUTPUT_DIR, "plots", "tss_enrichment_plot.png"),
       p_tss, width = 8, height = 6, dpi = 300)

# Nucleosome banding pattern - should show clear mono/di-nucleosomal banding
p_nucl <- FragmentHistogram(atac_obj)
ggsave(file.path(OUTPUT_DIR, "plots", "nucleosome_histogram.png"),
       p_nucl, width = 8, height = 5, dpi = 300)

n_before <- ncol(atac_obj)

# Apply QC filters - adjust thresholds based on VlnPlot inspection above
atac_obj <- subset(
  atac_obj,
  subset = nCount_peaks       >= MIN_COUNTS        &
           nCount_peaks       <= MAX_COUNTS        &
           TSS.enrichment     >= MIN_TSS           &
           nucleosome_signal  <= MAX_NUCL_SIGNAL   &
           pct_reads_in_peaks >= PEAK_PCT_MIN      &
           blacklist_ratio    <= 0.05
)

n_after <- ncol(atac_obj)
cat("Cells before QC:", n_before, "\n")
cat("Cells after QC: ", n_after,  "\n")
cat("Cells removed:  ", n_before - n_after, "\n")

qc_summary <- data.frame(
  Metric = c("Cells before QC", "Cells after QC", "Cells removed",
             "nCount_peaks min", "nCount_peaks max",
             "TSS enrichment min", "Nucleosome signal max",
             "Min pct reads in peaks"),
  Value  = c(n_before, n_after, n_before - n_after,
             MIN_COUNTS, MAX_COUNTS, MIN_TSS,
             MAX_NUCL_SIGNAL, PEAK_PCT_MIN)
)
write.csv(qc_summary,
          file.path(OUTPUT_DIR, "tables", "qc_summary.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# SECTION 4: NORMALIZATION - TF-IDF
# ------------------------------------------------------------------------------
cat("\n=== SECTION 4: TF-IDF NORMALIZATION ===\n")

# TF-IDF: term frequency–inverse document frequency normalization
# Accounts for variable sequencing depth (TF) and downweights ubiquitous
# peaks open across many cells (IDF) - standard for sparse ATAC data.
atac_obj <- RunTFIDF(atac_obj)

# Feature selection: peaks with high variability across cells
# Excludes very low-frequency peaks (< 5% of cells) - likely noise
atac_obj <- FindTopFeatures(atac_obj, min.cutoff = "q5")

cat("Variable features selected:", length(VariableFeatures(atac_obj)), "\n")

# ------------------------------------------------------------------------------
# SECTION 5: DIMENSIONALITY REDUCTION - LSI + UMAP + CLUSTERING
# ------------------------------------------------------------------------------
cat("\n=== SECTION 5: DIMENSIONALITY REDUCTION + CLUSTERING ===\n")

# SVD / Latent Semantic Indexing - standard dim reduction for ATAC-seq
atac_obj <- RunSVD(atac_obj)

# CRITICAL: Check correlation of each LSI component with sequencing depth.
# Component 1 is almost always dominated by depth - must be excluded.
# Inspect the plot to confirm before proceeding.
depth_cor <- DepthCor(atac_obj)
ggsave(file.path(OUTPUT_DIR, "plots", "lsi_depth_correlation.png"),
       depth_cor, width = 8, height = 5, dpi = 300)
cat("LSI depth correlation plotted - verify dim 1 excluded in LSI_DIMS_USE\n")

# UMAP + clustering - using LSI_DIMS_USE (dim 1 excluded - see SECTION 0)
atac_obj <- RunUMAP(atac_obj, reduction = "lsi", dims = LSI_DIMS_USE)
atac_obj <- FindNeighbors(atac_obj, reduction = "lsi", dims = LSI_DIMS_USE)
atac_obj <- FindClusters(atac_obj, resolution = CLUSTER_RES,
                         algorithm = 3,    # SLM - preferred for ATAC data
                         verbose = FALSE)

cat("Clusters identified:", nlevels(atac_obj$seurat_clusters), "\n")

p_umap <- DimPlot(atac_obj, label = TRUE, pt.size = 0.5) +
          ggtitle("UMAP - Clusters")
ggsave(file.path(OUTPUT_DIR, "plots", "umap_clusters.png"),
       p_umap, width = 8, height = 7, dpi = 300)

saveRDS(atac_obj, file.path(OUTPUT_DIR, "atac_clustered.rds"))
cat("Clustering checkpoint saved.\n")

# ------------------------------------------------------------------------------
# SECTION 6: PEAK CALLING - MACS2
# ------------------------------------------------------------------------------
# Re-call peaks per cluster using MACS2 on the final filtered cell set.
# Cell Ranger peaks are called on all cells pooled; per-cluster calling
# recovers cell-type-specific accessible regions missed by pooled calling.
# Requires MACS2: pip install macs2
# ------------------------------------------------------------------------------
cat("\n=== SECTION 6: MACS2 PEAK CALLING ===\n")

tryCatch({
  peaks <- CallPeaks(
    atac_obj,
    group.by        = "seurat_clusters",
    macs2.path      = Sys.which("macs2"),
    outdir          = file.path(OUTPUT_DIR, "peaks"),
    additional.args = paste("-g", MACS2_GENOME)   # genome size flag for MACS2
  )

  # Remove blacklist peaks
  peaks <- keepSeqlevels(peaks, value = seqlevels(genome_bsgenome),
                         pruning.mode = "coarse")
  peaks <- peaks[!overlapsAny(peaks, blacklist_gr)]
  cat("Peaks after blacklist removal:", length(peaks), "\n")

  # Quantify fragments in new MACS2 peaks
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(atac_obj),
    features  = peaks,
    cells     = colnames(atac_obj)
  )

  # Add MACS2 peaks as a new assay
  atac_obj[["MACS2"]] <- CreateChromatinAssay(
    counts     = macs2_counts,
    fragments  = FRAGMENT_FILE,
    annotation = Annotation(atac_obj)
  )
  DefaultAssay(atac_obj) <- "MACS2"
  cat("MACS2 assay set as default.\n")

  # Re-normalize and re-cluster on MACS2 peaks
  atac_obj <- RunTFIDF(atac_obj)
  atac_obj <- FindTopFeatures(atac_obj, min.cutoff = "q5")
  atac_obj <- RunSVD(atac_obj)
  atac_obj <- RunUMAP(atac_obj, reduction = "lsi", dims = LSI_DIMS_USE,
                      reduction.name = "umap_macs2")

  p_macs2 <- DimPlot(atac_obj, reduction = "umap_macs2", label = TRUE) +
             ggtitle("UMAP - MACS2 Peaks")
  ggsave(file.path(OUTPUT_DIR, "plots", "umap_macs2_peaks.png"),
         p_macs2, width = 8, height = 7, dpi = 300)

}, error = function(e) {
  cat("MACS2 not available or failed:", e$message, "\n")
  cat("Proceeding with Cell Ranger peaks.\n")
})

# ------------------------------------------------------------------------------
# SECTION 7: GENE ACTIVITY SCORES
# ------------------------------------------------------------------------------
# Gene activity: sums ATAC signal in gene body + 2kb promoter region per cell.
# Used as a proxy for gene expression - enables comparison with scRNA-seq
# and label transfer between modalities.
# ------------------------------------------------------------------------------
cat("\n=== SECTION 7: GENE ACTIVITY SCORES ===\n")

gene_activities <- GeneActivity(atac_obj)

atac_obj[["GeneActivity"]] <- CreateAssayObject(counts = gene_activities)
atac_obj <- NormalizeData(
  atac_obj,
  assay                = "GeneActivity",
  normalization.method = "LogNormalize",
  scale.factor         = median(atac_obj$nCount_peaks)
)

cat("Gene activity scores computed for", nrow(gene_activities), "genes.\n")

# Visualise canonical accessibility markers
# Replace with markers relevant to your cell types of interest
canonical_markers <- c("CD14", "CD3D", "MS4A1", "GATA1", "PAX5")
canonical_present <- canonical_markers[
  canonical_markers %in% rownames(atac_obj[["GeneActivity"]])]

if (length(canonical_present) > 0) {
  p_ga <- FeaturePlot(
    atac_obj,
    features   = canonical_present,
    reduction  = "umap",
    max.cutoff = "q95",
    ncol       = 3
  )
  ggsave(file.path(OUTPUT_DIR, "plots", "gene_activity_canonical_markers.png"),
         p_ga, width = 15, height = 10, dpi = 300)
}

# ------------------------------------------------------------------------------
# SECTION 8: MOTIF ANALYSIS - chromVAR
# ------------------------------------------------------------------------------
# chromVAR computes per-cell deviation scores for JASPAR TF motifs.
# Captures TF activity differences across clusters - more informative than
# peak accessibility alone for biological interpretation.
# ------------------------------------------------------------------------------
cat("\n=== SECTION 8: MOTIF ANALYSIS (chromVAR) ===\n")

DefaultAssay(atac_obj) <- ifelse("MACS2" %in% Assays(atac_obj), "MACS2", "peaks")

tryCatch({
  # JASPAR 2020 motifs for human
  pfm <- getMatrixSet(
    JASPAR2020,
    opts = list(species = "Homo sapiens", all_versions = FALSE)
  )

  atac_obj <- AddMotifs(
    object = atac_obj,
    genome = genome_bsgenome,
    pfm    = pfm
  )

  atac_obj <- RunChromVAR(
    object = atac_obj,
    genome = genome_bsgenome
  )

  DefaultAssay(atac_obj) <- "chromvar"
  motif_markers <- FindAllMarkers(
    atac_obj,
    only.pos        = TRUE,
    min.pct         = 0.1,
    logfc.threshold = 0.1,
    verbose         = FALSE
  )

  write.csv(motif_markers,
            file.path(OUTPUT_DIR, "tables", "chromvar_motif_markers.csv"),
            row.names = FALSE)
  cat("chromVAR motif markers saved.\n")

  top_motifs <- motif_markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 3)
  cat("Top motifs per cluster:\n"); print(top_motifs)

}, error = function(e) {
  cat("chromVAR failed:", e$message, "\n")
  cat("Ensure BSgenome package matches GENOME parameter.\n")
})

# ------------------------------------------------------------------------------
# SECTION 9: DIFFERENTIAL ACCESSIBILITY
# ------------------------------------------------------------------------------
# Logistic regression with nCount_peaks as latent variable - recommended
# for ATAC-seq as it directly regresses out sequencing depth confounding.
# ------------------------------------------------------------------------------
cat("\n=== SECTION 9: DIFFERENTIAL ACCESSIBILITY ===\n")

DefaultAssay(atac_obj) <- ifelse("MACS2" %in% Assays(atac_obj), "MACS2", "peaks")
Idents(atac_obj) <- "seurat_clusters"

da_peaks <- FindAllMarkers(
  atac_obj,
  only.pos        = TRUE,
  min.pct         = 0.05,
  logfc.threshold = 0.1,
  verbose         = FALSE,
  test.use        = "LR",           # logistic regression - recommended for ATAC
  latent.vars     = "nCount_peaks"  # regress out depth differences
)

# Annotate DA peaks with nearest gene
closest_genes <- ClosestFeature(atac_obj, regions = rownames(da_peaks))
da_peaks$closest_gene  <- closest_genes$gene_name[
  match(rownames(da_peaks), closest_genes$query_region)]
da_peaks$gene_distance <- closest_genes$distance[
  match(rownames(da_peaks), closest_genes$query_region)]

write.csv(da_peaks,
          file.path(OUTPUT_DIR, "tables", "differential_accessibility_peaks.csv"))
cat("Differential accessibility results saved.\n")

# ------------------------------------------------------------------------------
# SECTION 10: OPTIONAL - Label Transfer from scRNA-seq
# ------------------------------------------------------------------------------
if (!is.null(RNASEQ_REF_RDS) && file.exists(RNASEQ_REF_RDS)) {
  cat("\n=== SECTION 10: LABEL TRANSFER FROM scRNA-seq ===\n")

  rna_obj <- readRDS(RNASEQ_REF_RDS)
  DefaultAssay(atac_obj) <- "GeneActivity"

  transfer_anchors <- FindTransferAnchors(
    reference       = rna_obj,
    query           = atac_obj,
    features        = VariableFeatures(rna_obj),
    reference.assay = "RNA",
    query.assay     = "GeneActivity",
    reduction       = "pcaproject",
    verbose         = FALSE
  )

  predicted_labels <- TransferData(
    anchorset        = transfer_anchors,
    refdata          = rna_obj[[RNA_LABEL_COL]][, 1],
    weight.reduction = atac_obj[["lsi"]],
    dims             = LSI_DIMS_USE,
    verbose          = FALSE
  )

  atac_obj <- AddMetaData(atac_obj, metadata = predicted_labels)

  p_transfer <- DimPlot(atac_obj, group.by = "predicted.id",
                        label = TRUE, repel = TRUE) +
               ggtitle("Label Transfer from scRNA-seq")
  ggsave(file.path(OUTPUT_DIR, "plots", "label_transfer_predictions.png"),
         p_transfer, width = 10, height = 8, dpi = 300)

  write.csv(predicted_labels,
            file.path(OUTPUT_DIR, "tables", "label_transfer_predictions.csv"))
  cat("Label transfer complete.\n")

} else {
  cat("\nSection 10 skipped - no scRNA-seq reference provided.\n")
  cat("Set RNASEQ_REF_RDS in SECTION 0 to enable label transfer.\n")
}

# ------------------------------------------------------------------------------
# SECTION 11: SAVE FINAL OBJECT + COMPOSITION TABLE
# ------------------------------------------------------------------------------
cat("\n=== SECTION 11: SAVING OUTPUTS ===\n")

DefaultAssay(atac_obj) <- ifelse("MACS2" %in% Assays(atac_obj), "MACS2", "peaks")
saveRDS(atac_obj, file.path(OUTPUT_DIR, "atac_final.rds"))

cluster_comp <- as.data.frame(table(atac_obj$seurat_clusters))
colnames(cluster_comp) <- c("Cluster", "N_Cells")
cluster_comp$Pct <- round(100 * cluster_comp$N_Cells / sum(cluster_comp$N_Cells), 1)
write.csv(cluster_comp,
          file.path(OUTPUT_DIR, "tables", "cluster_composition.csv"),
          row.names = FALSE)

cat("\n=== PIPELINE COMPLETE ===\n")
cat("Final cells:", ncol(atac_obj), "\n")
cat("Clusters:", nlevels(atac_obj$seurat_clusters), "\n")
cat("Outputs in:", OUTPUT_DIR, "\n")
