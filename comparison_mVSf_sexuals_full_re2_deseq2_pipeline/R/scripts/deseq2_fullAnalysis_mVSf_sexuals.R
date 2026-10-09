###############################################################################
## DESeq2 — mVSf_sexuals: sexual males + sexual females (29 samples)
##
## Input : inputs/mVSf_sexuals_cluster_matrix.tsv  (cluster x sample; mapped onto geo's RE2 reference)
##         inputs/mVSf_metadata.csv               (sample, population, sex, region)
## Output: outputs/  (all tables, figures, dds/res .rds, sessionInfo)
##
## Purpose: validate the female baseline. Sexual females should track sexual males,
## which justifies comparing females against the geo (male) reference.
##  1. Primary test: females only, Swiss vs Austrian (~ region, Austrian = reference).
##  2. Descriptive PCA of all 29 samples, coloured by sex / population / region.
##
## COUNTS: the HPC matrix holds SUMMED PER-BASE COVERAGE per cluster. Every read is
## trimmed to exactly 140 bp and aligned end-to-end, so coverage / 140 = mapped reads
## (mates) per cluster. DESeq2's model assumes read counts, so the matrix is converted
## to read equivalents before anything else.
##
## FILTER: prevalence filter (DESeq2 vignette pattern). A cluster is kept only if it
## has >= minPerSample reads in >= N samples, N = size of the smallest group in THAT
## test (computed per test).
##
## METADATA NOTE: the metadata lists PSY_0062, which has no column in the matrix;
## it is reported and ignored.
##
## HOW TO RUN: set the working directory to the mVSf_sexuals/ folder
## (RStudio: Session > Set Working Directory > Choose Directory... > mVSf_sexuals),
## open this script and press Source (runs top to bottom).
###############################################################################


################ 0. SETUP -----------------------------------------------------

cran_packages <- c("tidyverse", "cowplot", "RColorBrewer", "ggrepel")
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

options(stringsAsFactors = FALSE)


################ KNOBS --------------------------------------------------------

inDir        <- "inputs"
outDir       <- "outputs"
countsFile   <- file.path(inDir, "mVSf_sexuals_cluster_matrix.tsv")
metadataFile <- file.path(inDir, "mVSf_metadata.csv")
outPrefix    <- "mVSf"

austrianPops <- c("Linz", "Pfenningberg", "Luftenberg")   # shape split in the population PCA
readLength   <- 140    # cutadapt --length/--minimum-length: coverage / readLength = reads
minPerSample <- 10     # prevalence filter: cluster "present" in a sample at >= this many reads
padjCut      <- 0.05   # significance threshold
lfcLabel     <- 2      # volcano: label clusters with |log2FoldChange| > this
padjLabel    <- 30     # volcano: label clusters with -log10(padj) > this

if (!file.exists(countsFile)) {
  stop("Input not found: ", countsFile, ". Set the working directory to mVSf_sexuals/ (current: ", getwd(), ")")
}
dir.create(outDir, showWarnings = FALSE)

savePlot <- function(p, name, w = 10, h = 8) {
  ggsave(file.path(outDir, paste0(outPrefix, "_", name, ".png")), p,
         bg = "white", width = w, height = h, dpi = 300)
}

## prevalence filter: >= minPerSample reads in >= (smallest group size) samples
prevalenceFilter <- function(dds, group, label) {
  minSamples <- min(table(colData(dds)[[group]]))
  keep <- rowSums(counts(dds) >= minPerSample) >= minSamples
  message(label, ": prevalence filter (>= ", minPerSample, " reads in >= ", minSamples,
          " samples) kept ", sum(keep), " / ", length(keep), " clusters")
  dds[keep, ]
}

## colourblind-safe palettes
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#000000")
tol_9     <- c("#332288", "#117733", "#44AA99", "#88CCEE",
               "#DDCC77", "#CC6677", "#AA4499", "#882255", "#EE7733")


################ 1. COUNTS + CONVERT TO READ EQUIVALENTS ---------------------

counts <- read_tsv(countsFile, show_col_types = FALSE) %>% as.data.frame()
rownames(counts) <- counts$Cluster
counts <- counts[, -1]                     # Cluster column now lives in rownames

counts <- round(counts / readLength)       # summed per-base coverage -> mapped reads

dim(counts)                                # expect 29 sample columns


################ 2. METADATA (colData) ----------------------------------------

coldata <- read_csv(metadataFile, show_col_types = FALSE) %>% as.data.frame()
rownames(coldata) <- coldata$sample
coldata <- coldata[, -1, drop = FALSE]     # drop redundant sample column

coldata$population <- factor(coldata$population)
coldata$sex        <- factor(coldata$sex,    levels = c("male", "female"))
coldata$region     <- factor(coldata$region, levels = c("Austrian", "Swiss"))   # Austrian = reference

## every matrix sample must have metadata; extra metadata rows are allowed and ignored
stopifnot(all(colnames(counts) %in% rownames(coldata)))
extra <- setdiff(rownames(coldata), colnames(counts))
if (length(extra) > 0) message("Metadata samples not in the matrix (ignored): ", paste(extra, collapse = ", "))
coldata <- coldata[colnames(counts), , drop = FALSE]
coldata$population <- droplevels(coldata$population)

table(coldata$sex, coldata$region)


################ 3. PRIMARY — females only, Swiss vs Austrian -----------------

female_samples <- rownames(coldata)[coldata$sex == "female"]
counts_f  <- counts[, female_samples]
coldata_f <- coldata[female_samples, ]
coldata_f$population <- droplevels(coldata_f$population)
coldata_f$region     <- droplevels(coldata_f$region)

## sex is constant within females, so the model is region only
dds_f <- DESeqDataSetFromMatrix(countData = counts_f, colData = coldata_f, design = ~ region)
dds_f <- prevalenceFilter(dds_f, "region", "Females Swiss vs Austrian")
dds_f <- DESeq(dds_f)

res_f <- results(dds_f, contrast = c("region", "Swiss", "Austrian"))

cat("\n===== Summary: Swiss vs Austrian (females only) =====\n")
summary(res_f, alpha = 0.10)
summary(res_f, alpha = 0.05)

res_f_df  <- as.data.frame(res_f) %>% rownames_to_column("cluster") %>% arrange(padj)
res_f_sig <- res_f_df %>% filter(!is.na(padj), padj < padjCut)

write_tsv(res_f_df,  file.path(outDir, "results_females_SwissVsAustrian.tsv"))
write_tsv(res_f_sig, file.path(outDir, "results_females_SwissVsAustrian_sig.tsv"))
saveRDS(list(dds = dds_f, res = res_f), file.path(outDir, "mVSf_females_SwissVsAustrian_dds_res.rds"))

cat("Sig:", nrow(res_f_sig), "of", sum(!is.na(res_f$padj)), "tested\n")


################ 3.1 VOLCANO — females, outliers labelled ---------------------

volc_f_df <- res_f_df %>%
  filter(!is.na(padj)) %>%
  mutate(label = ifelse(abs(log2FoldChange) > lfcLabel | -log10(padj) > padjLabel,
                        cluster, NA))

volc_f <- ggplot(volc_f_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.6) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20, na.rm = TRUE) +
  theme_cowplot(17) +
  xlab("Austrian \u2192 Swiss (log2 fold-change)") +
  ylab("-log10 adjusted p-value") +
  ggtitle("Volcano: Swiss vs Austrian (females)")
print(volc_f)
savePlot(volc_f, "volcano_females_SwissVsAustrian")


################ 3.2 LFC — females, significant clusters only -----------------

lfc_f <- res_f_sig %>%
  mutate(cl_num = as.numeric(gsub("CL", "", cluster))) %>%
  ggplot(aes(x = cl_num, y = log2FoldChange)) +
  geom_point() +
  theme_cowplot(17) +
  xlab("Cluster") +
  ylab("log2 fold-change (Swiss vs Austrian)") +
  ggtitle(paste0("LFC: Swiss vs Austrian, females (padj < ", padjCut, ")"))
print(lfc_f)
savePlot(lfc_f, "LFC_females_SwissVsAustrian")


################ 3.3 PCA — females only ---------------------------------------

rld_f <- rlog(dds_f)

pca_f_data   <- plotPCA(rld_f, intgroup = "region", returnData = TRUE)
percentVar_f <- round(100 * attr(pca_f_data, "percentVar"))

pca_f <- ggplot(pca_f_data, aes(x = PC1, y = PC2, color = region, shape = region)) +
  geom_point(size = 4) +
  scale_color_manual(values = okabe_ito) +
  scale_shape_manual(values = c(17, 16)) +   # triangle / circle
  xlab(paste0("PC1: ", percentVar_f[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_f[2], "% variance")) +
  theme_cowplot(16) +
  ggtitle("PCA: females only (Swiss vs Austrian)") +
  theme(panel.grid = element_blank())
print(pca_f)
savePlot(pca_f, "PCA_females_SwissVsAustrian")


################ 4. DESCRIPTIVE PCA — all 29 samples --------------------------

## design ~ 1: no model, the dds is only used for the rlog transform.
## Filter uses the smallest region group across all 29 samples.
dds_all <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ 1)
dds_all <- prevalenceFilter(dds_all, "region", "All samples (descriptive PCA)")
rld_all <- rlog(dds_all)

## by sex
pca_sex_data   <- plotPCA(rld_all, intgroup = "sex", returnData = TRUE)
percentVar_all <- round(100 * attr(pca_sex_data, "percentVar"))

pca_sex <- ggplot(pca_sex_data, aes(x = PC1, y = PC2, color = sex)) +
  geom_point(size = 4) +
  scale_color_manual(values = okabe_ito) +
  xlab(paste0("PC1: ", percentVar_all[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_all[2], "% variance")) +
  theme_cowplot(16) +
  ggtitle("PCA: all samples (coloured by sex)")
print(pca_sex)
savePlot(pca_sex, "PCA_allSamples_bySex")

## by population, shape by country
pca_pop_data <- plotPCA(rld_all, intgroup = "population", returnData = TRUE)
pca_pop_data$country <- ifelse(pca_pop_data$population %in% austrianPops, "Austrian", "Swiss")

pca_pop <- ggplot(pca_pop_data, aes(x = PC1, y = PC2, color = population, shape = country)) +
  geom_point(size = 4) +
  scale_color_manual(values = tol_9) +
  scale_shape_manual(values = c("Austrian" = 17, "Swiss" = 16)) +
  xlab(paste0("PC1: ", percentVar_all[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_all[2], "% variance")) +
  theme_cowplot(16) +
  ggtitle("PCA: all samples (coloured by population)")
print(pca_pop)
savePlot(pca_pop, "PCA_allSamples_byPopulation")

## by region
pca_region_data <- plotPCA(rld_all, intgroup = "region", returnData = TRUE)

pca_region <- ggplot(pca_region_data, aes(x = PC1, y = PC2, color = region)) +
  geom_point(size = 4) +
  scale_color_manual(values = okabe_ito) +
  xlab(paste0("PC1: ", percentVar_all[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_all[2], "% variance")) +
  theme_cowplot(16) +
  ggtitle("PCA: all samples (coloured by region)")
print(pca_region)
savePlot(pca_region, "PCA_allSamples_byRegion")


################ 5. SESSION INFO ----------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(outDir, "sessionInfo.txt"))
message("Done. All outputs in ", normalizePath(outDir))