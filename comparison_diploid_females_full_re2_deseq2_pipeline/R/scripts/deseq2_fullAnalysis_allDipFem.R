###############################################################################
## DESeq2 — allDiploidFemales: asexual vs sexual, PSY_0145 EXCLUDED
##
## Input : inputs/allDiploidFemales_cluster_matrix.tsv  (cluster x sample, 18 samples)
##         inputs/allDiploidFemales_metadata.csv        (sample, sex, location)
## Output: outputs/  (all tables, figures, dds/res .rds, sessionInfo)
##
## Design ~ sex, sexual = reference. 4 sexual vs 13 asexual after dropping PSY_0145.
##
## COUNTS: the HPC matrix holds SUMMED PER-BASE COVERAGE per cluster. Every read is
## trimmed to exactly 140 bp and aligned end-to-end, so coverage / 140 = mapped reads
## (mates) per cluster. DESeq2's model assumes read counts, so the matrix is converted
## to read equivalents before anything else.
##
## FILTER: prevalence filter (DESeq2 vignette pattern). A cluster is kept only if it
## has >= minPerSample reads in >= minSamples samples, minSamples = smallest group (4).
## Clusters carried by 1-3 libraries cannot inform a group contrast (e.g. the
## PSY_0023-driven clusters: 2-3 non-zero samples, 0023 = 86-99% of counts).
##
## WHY PSY_0145 IS EXCLUDED:
##  - PRIMARY (a priori, taxonomic): barcode identification places PSY_0145 in
##    Pseudobankesia (P. alpestrella or P. vernella), NOT Dahlica triquetrella.
##  - CONSISTENT QC SIGNATURE (earlier 18-sample run, not shipped): lone outlier in the
##    sample-correlation heatmap (~0.4-0.5 with all others; others 0.85-1.0), and sole
##    driver of a high-LFC volcano arm (e.g. CL182: PSY_0145 = 1.65e10 normalized counts,
##    11/14 asexuals exactly 0).
##  - Taxonomic exclusion list, not the technical one (0139/0173/0174/0178-0181).
##  - NOTE: PSY_0145 reads were part of the RE2 clustering input, so the reference
##    contains some non-target clusters. The prevalence filter removes their effect
##    on the test; reclustering without 0145 is open work (see HANDOFF).
##
## HOW TO RUN: set the working directory to the allDiploidFemales/ folder
## (RStudio: Session > Set Working Directory > Choose Directory... > allDiploidFemales),
## open this script and press Source (runs top to bottom).
###############################################################################


################ 0. SETUP -----------------------------------------------------

cran_packages <- c("tidyverse", "cowplot", "RColorBrewer", "ggrepel", "pheatmap")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 missing. Install with: BiocManager::install('DESeq2')")
}

library(tidyverse)     # includes readr, dplyr, tidyr, ggplot2
library(DESeq2)
library(cowplot)
library(RColorBrewer)
library(ggrepel)
library(pheatmap)

options(stringsAsFactors = FALSE)


################ KNOBS --------------------------------------------------------

inDir        <- "inputs"
outDir       <- "outputs"
countsFile   <- file.path(inDir, "allDiploidFemales_cluster_matrix.tsv")   # 1st col Cluster, then samples
metadataFile <- file.path(inDir, "allDiploidFemales_metadata.csv")         # sample, sex, location
outPrefix    <- "allDiploidFemales_no0145"
contrastTag  <- "asexualVsSexual"

dropSamples  <- c("PSY_0145")   # Pseudobankesia (barcode), not D. triquetrella

readLength   <- 140    # cutadapt --length/--minimum-length: coverage / readLength = reads
minPerSample <- 10     # prevalence filter: cluster "present" in a sample at >= this many reads
minSamples   <- 4      # ...in >= this many samples (= smallest group: 4 sexuals)
padjCut      <- 0.05   # threshold for calling a cluster significant
lfcLabel     <- 2      # volcano: label clusters with |log2FoldChange| > this
padjLabel    <- 2      # volcano: label clusters with -log10(padj) > this

if (!file.exists(countsFile)) {
  stop("Input not found: ", countsFile, ". Set the working directory to allDiploidFemales/ (current: ", getwd(), ")")
}
dir.create(outDir, showWarnings = FALSE)


################ 1. LOAD COUNTS + CONVERT TO READ EQUIVALENTS ----------------

counts <- read_tsv(countsFile, show_col_types = FALSE) %>% as.data.frame()
rownames(counts) <- counts$Cluster
counts <- counts[, -1]                     # Cluster column now lives in rownames

counts <- round(counts / readLength)       # summed per-base coverage -> mapped reads

dim(counts)


################ 2. METADATA (colData) ----------------------------------------

coldata <- read_csv(metadataFile, show_col_types = FALSE) %>% as.data.frame()
rownames(coldata) <- coldata$sample
coldata <- coldata[, -1, drop = FALSE]     # drop redundant sample column

coldata$sex <- factor(coldata$sex, levels = c("sexual", "asexual"))  # sexual = reference


################ 2.1 DROP EXCLUDED SAMPLE(S) ----------------------------------

message("Dropping: ", paste(dropSamples, collapse = ", "))
coldata <- coldata[!(rownames(coldata) %in% dropSamples), , drop = FALSE]
counts  <- counts[, !(colnames(counts) %in% dropSamples)]
coldata$sex <- droplevels(coldata$sex)

## Matrix columns and metadata rows must match exactly, in the same order.
stopifnot(setequal(colnames(counts), rownames(coldata)))
coldata <- coldata[colnames(counts), , drop = FALSE]

table(coldata$sex)                         # expect 4 sexual vs 13 asexual


################ 3. ANALYSIS — asexual vs sexual ------------------------------

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = coldata,
                              design    = ~ sex)

## prevalence filter
keep <- rowSums(counts(dds) >= minPerSample) >= minSamples
dds  <- dds[keep, ]
message("Clusters kept after prevalence filter: ", sum(keep), " / ", length(keep))

dds <- DESeq(dds)

## Default alpha (0.1) for DESeq2's independent filtering; significance called at padjCut.
res <- results(dds, contrast = c("sex", "asexual", "sexual"))

cat("\n===== Summary: asexual vs sexual (diploid females, 0145 excluded) =====\n")
summary(res, alpha = 0.10)
summary(res, alpha = 0.05)


################ 3.0 EXPORT — single source of truth --------------------------

## res_df feeds the TSVs AND every plot below. Nothing downstream reads `res`.
res_df <- as.data.frame(res) %>%
  rownames_to_column("cluster") %>%
  arrange(padj)

res_sig <- res_df %>% filter(!is.na(padj), padj < padjCut)

cat("\nSig at padj <", padjCut, ": ", nrow(res_sig), " of ", sum(!is.na(res$padj)), " tested",
    " (up in asexual: ", sum(res_sig$log2FoldChange > 0),
    ", down: ", sum(res_sig$log2FoldChange < 0), ")\n", sep = "")

fullFile <- file.path(outDir, paste0("results_", outPrefix, "_", contrastTag, ".tsv"))
sigFile  <- file.path(outDir, paste0("results_", outPrefix, "_", contrastTag, "_sig.tsv"))
write_tsv(res_df,  fullFile)
write_tsv(res_sig, sigFile)

saveRDS(list(dds = dds, res = res),
        file.path(outDir, paste0(outPrefix, "_", contrastTag, "_dds_res.rds")))

## Read-back check: stop loudly if the files on disk don't match the live res.
stopifnot(
  nrow(read_tsv(fullFile, show_col_types = FALSE)) == nrow(res),
  nrow(read_tsv(sigFile,  show_col_types = FALSE)) == sum(res$padj < padjCut, na.rm = TRUE)
)
message("Exported ", nrow(res_df), " clusters; ", nrow(res_sig),
        " sig (padj < ", padjCut, ") -> ", outDir)


################ 3.1 rlog transform (computed once, reused by QC + PCA) --------

rld <- rlog(dds)


################ 3.2 QC — dispersion + sample correlation ---------------------

png(file.path(outDir, paste0(outPrefix, "_QC_dispersion.png")),
    width = 1800, height = 1400, res = 200)
plotDispEsts(dds)
dev.off()

## Should be uniformly high-correlation; watch for a NEW outlier.
sampleCor <- cor(assay(rld))
annot_col <- as.data.frame(colData(dds)[, "sex", drop = FALSE])
pheatmap(sampleCor,
         annotation_col = annot_col,
         main     = "Sample correlation (rlog) — diploid females (no 0145)",
         filename = file.path(outDir, paste0(outPrefix, "_QC_sampleCorr.png")),
         width = 8, height = 7)


################ 3.3 PCA — by sex --------------------------------------------

sex_palette <- c("sexual" = "#0072B2", "asexual" = "#D55E00")

pca_data <- plotPCA(rld, intgroup = "sex", returnData = TRUE)
pct      <- round(100 * attr(pca_data, "percentVar"))

pca_sex <- ggplot(pca_data, aes(x = PC1, y = PC2, color = sex, shape = sex)) +
  geom_point(size = 4) +
  scale_color_manual(values = sex_palette) +
  scale_shape_manual(values = c("sexual" = 16, "asexual" = 17)) +
  xlab(paste0("PC1: ", pct[1], "% variance")) +
  ylab(paste0("PC2: ", pct[2], "% variance")) +
  theme_cowplot(16) +
  ggtitle("PCA: diploid females (sexual vs asexual, no 0145)") +
  theme(panel.grid = element_blank())
print(pca_sex)
ggsave(file.path(outDir, paste0(outPrefix, "_PCA_bySex.png")), pca_sex,
       bg = "white", width = 10, height = 8, dpi = 300)


################ 3.4 VOLCANO (from res_df) ------------------------------------

volc_df <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(sig   = padj < padjCut,
         label = ifelse(abs(log2FoldChange) > lfcLabel | -log10(padj) > padjLabel,
                        cluster, NA))

volc <- ggplot(volc_df, aes(x = log2FoldChange, y = -log10(padj), colour = sig)) +
  geom_point(alpha = 0.6) +
  scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "#D55E00"),
                      name = paste0("padj < ", padjCut)) +
  geom_hline(yintercept = -log10(padjCut), linetype = "dashed", colour = "grey40") +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20,
                  colour = "black", na.rm = TRUE) +
  theme_cowplot(17) +
  xlab("sexual \u2192 asexual (log2 fold-change)") +
  ylab("-log10 adjusted p-value") +
  ggtitle(paste0("Volcano: asexual vs sexual (diploid females, no 0145) — ",
                 nrow(res_sig), " sig"))
print(volc)
ggsave(file.path(outDir, paste0(outPrefix, "_volcano.png")), volc,
       bg = "white", width = 10, height = 8, dpi = 300)


################ 3.5 LFC PLOT — significant clusters only (from res_sig) ------

lfc <- res_sig %>%
  mutate(cl_num = as.numeric(gsub("CL", "", cluster))) %>%
  ggplot(aes(x = cl_num, y = log2FoldChange)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_point() +
  geom_text_repel(aes(label = cluster), size = 3) +
  theme_cowplot(17) +
  xlab("Cluster") +
  ylab("log2 fold-change (asexual vs sexual)") +
  ggtitle(paste("LFC: asexual vs sexual, no 0145 (padj <", padjCut, ")"))
print(lfc)
ggsave(file.path(outDir, paste0(outPrefix, "_LFC.png")), lfc,
       bg = "white", width = 10, height = 8, dpi = 300)


################ 4. SESSION INFO ----------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(outDir, "sessionInfo.txt"))
message("Done. All outputs in ", normalizePath(outDir))