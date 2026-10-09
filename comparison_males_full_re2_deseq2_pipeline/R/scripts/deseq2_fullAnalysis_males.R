###############################################################################
## DESeq2 — geo (males): Swiss vs Austrian (Linz) repeat-cluster abundance
##
## Input : inputs/geo_cluster_matrix.tsv  (cluster x sample, from the HPC matrix step)
##         inputs/geo_metadata.csv        (sample, population)
## Output: outputs/  (all tables, figures, dds/res .rds, sessionInfo)
##
## Two analyses, design ~ geo_group, Linz = reference:
##  1. Overall: all Swiss males vs Linz.
##  2. Population-wise: each Swiss population (Nods, SwGuggis, RehetMid) vs Linz,
##     then clusters significant in all three / in at least two populations.
##
## COUNTS: the HPC matrix holds SUMMED PER-BASE COVERAGE per cluster. Every read is
## trimmed to exactly 140 bp and aligned end-to-end, so coverage / 140 = mapped reads
## (mates) per cluster. DESeq2's model assumes read counts, so the matrix is converted
## to read equivalents before anything else.
##
## FILTER: prevalence filter (DESeq2 vignette pattern). A cluster is kept only if it
## has >= minPerSample reads in >= N samples, N = size of the smallest group in THAT
## test (computed per test). Clusters carried by fewer libraries than the smallest
## group cannot inform a group contrast.
##
## HOW TO RUN: set the working directory to the geo/ folder
## (RStudio: Session > Set Working Directory > Choose Directory... > geo),
## open this script and press Source (runs top to bottom).
###############################################################################


################ 0. SETUP -----------------------------------------------------

cran_packages <- c("tidyverse", "cowplot", "RColorBrewer", "ggrepel", "ggVennDiagram")
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
library(ggVennDiagram)

options(stringsAsFactors = FALSE)


################ KNOBS --------------------------------------------------------

inDir        <- "inputs"
outDir       <- "outputs"
countsFile   <- file.path(inDir, "geo_cluster_matrix.tsv")
metadataFile <- file.path(inDir, "geo_metadata.csv")
outPrefix    <- "geo"

swissPops    <- c("Nods", "SwGuggis", "RehetMid")
readLength   <- 140    # cutadapt --length/--minimum-length: coverage / readLength = reads
minPerSample <- 10     # prevalence filter: cluster "present" in a sample at >= this many reads
padjCut      <- 0.05   # significance threshold
lfcLabel     <- 2      # volcano: label clusters with |log2FoldChange| > this
padjLabel    <- 30     # volcano: label clusters with -log10(padj) > this

if (!file.exists(countsFile)) {
  stop("Input not found: ", countsFile, ". Set the working directory to geo/ (current: ", getwd(), ")")
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


################ 1. COUNTS + CONVERT TO READ EQUIVALENTS ---------------------

counts <- read_tsv(countsFile, show_col_types = FALSE) %>% as.data.frame()
rownames(counts) <- counts$Cluster
counts <- counts[, -1]                     # Cluster column now lives in rownames

counts <- round(counts / readLength)       # summed per-base coverage -> mapped reads

dim(counts)


################ 2. METADATA (colData) ----------------------------------------

coldata <- read_csv(metadataFile, show_col_types = FALSE) %>%
  as.data.frame() %>%
  mutate(
    geo_group  = ifelse(population == "Linz", "Linz", "Swiss"),
    population = factor(population, levels = c("Linz", swissPops)),
    geo_group  = factor(geo_group,  levels = c("Linz", "Swiss"))   # Linz = reference
  )
rownames(coldata) <- coldata$sample
coldata <- coldata[, -1, drop = FALSE]     # drop redundant sample column

stopifnot(setequal(colnames(counts), rownames(coldata)))
coldata <- coldata[colnames(counts), , drop = FALSE]

table(coldata$population)
table(coldata$geo_group)


################ 3. OVERALL — all Swiss vs Linz -------------------------------

dds_all <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ geo_group)
dds_all <- prevalenceFilter(dds_all, "geo_group", "Overall")
dds_all <- DESeq(dds_all)

res_all <- results(dds_all, contrast = c("geo_group", "Swiss", "Linz"))

cat("\n===== Summary: all Swiss vs Linz =====\n")
summary(res_all, alpha = 0.10)
summary(res_all, alpha = 0.05)

res_all_df  <- as.data.frame(res_all) %>% rownames_to_column("cluster") %>% arrange(padj)
res_all_sig <- res_all_df %>% filter(!is.na(padj), padj < padjCut)

write_tsv(res_all_df,  file.path(outDir, "results_geoMales_SwissVsLinz.tsv"))
write_tsv(res_all_sig, file.path(outDir, "results_geoMales_SwissVsLinz_sig.tsv"))
saveRDS(list(dds = dds_all, res = res_all), file.path(outDir, "geo_SwissVsLinz_dds_res.rds"))

cat("Sig:", nrow(res_all_sig), "of", sum(!is.na(res_all$padj)), "tested\n")


################ 3.1 PCA — overall --------------------------------------------

rld_all <- rlog(dds_all)

country_palette <- c("Austrian" = "#CC6677", "Swiss" = "#332288")   # muted red / indigo

pca_data <- plotPCA(rld_all, intgroup = "population", returnData = TRUE)
pct      <- round(100 * attr(pca_data, "percentVar"))
pca_data$country <- ifelse(pca_data$population == "Linz", "Austrian", "Swiss")

pca_all <- ggplot(pca_data, aes(PC1, PC2, color = country, shape = country)) +
  geom_point(size = 3) +
  scale_color_manual(values = country_palette) +
  scale_shape_manual(values = c("Austrian" = 17, "Swiss" = 16)) +
  xlab(paste0("PC1: ", pct[1], "% variance")) +
  ylab(paste0("PC2: ", pct[2], "% variance")) +
  ggtitle("PCA: all samples (Swiss vs Austrian)") +
  theme_bw() +
  theme(panel.grid = element_blank())
print(pca_all)
savePlot(pca_all, "PCA_all", w = 7, h = 5)


################ 3.2 VOLCANO — overall, outliers labelled ---------------------

volc_all_df <- res_all_df %>%
  filter(!is.na(padj)) %>%
  mutate(label = ifelse(abs(log2FoldChange) > lfcLabel | -log10(padj) > padjLabel,
                        cluster, NA))

volc_all <- ggplot(volc_all_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.6) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20, na.rm = TRUE) +
  theme_cowplot(17) +
  xlab("Linz \u2192 Swiss (log2 fold-change)") +
  ylab("-log10 adjusted p-value") +
  ggtitle("Volcano: all Swiss vs Linz")
print(volc_all)
savePlot(volc_all, "volcano_all")


################ 4. POPULATION-WISE — each Swiss population vs Linz -----------

run_pop_DESeq <- function(swiss_pop) {

  pop_samples <- rownames(coldata)[coldata$population %in% c(swiss_pop, "Linz")]
  counts_tmp  <- counts[, pop_samples]
  coldata_tmp <- coldata[pop_samples, ]
  coldata_tmp$population <- droplevels(coldata_tmp$population)
  coldata_tmp$geo_group  <- droplevels(coldata_tmp$geo_group)

  dds_tmp <- DESeqDataSetFromMatrix(countData = counts_tmp, colData = coldata_tmp, design = ~ geo_group)
  dds_tmp <- prevalenceFilter(dds_tmp, "geo_group", paste(swiss_pop, "vs Linz"))
  dds_tmp <- DESeq(dds_tmp)

  res <- results(dds_tmp, contrast = c("geo_group", "Swiss", "Linz"))

  cat("\n===== Summary:", swiss_pop, "vs Linz =====\n")
  summary(res, alpha = 0.10)
  summary(res, alpha = 0.05)

  res_df <- as.data.frame(res) %>%
    rownames_to_column("clusterID") %>%
    mutate(cluster    = as.numeric(gsub("CL", "", clusterID)),
           population = swiss_pop)

  write_tsv(res_df, file.path(outDir, paste0("results_geo_", swiss_pop, "VsLinz.tsv")))

  list(res_df = res_df, dds = dds_tmp)
}

pop_out <- setNames(lapply(swissPops, run_pop_DESeq), swissPops)


################ 4.1 PCA + VOLCANO + LFC — per population ---------------------

make_pop_plots <- function(pop_name) {

  dds_tmp <- pop_out[[pop_name]]$dds
  res_df  <- pop_out[[pop_name]]$res_df

  pca_p <- plotPCA(rlog(dds_tmp), intgroup = "population") +
    ggtitle(paste("PCA:", pop_name, "vs Linz"))
  print(pca_p)
  savePlot(pca_p, paste0("PCA_", pop_name, "VsLinz"))

  volc_p <- res_df %>%
    filter(!is.na(padj)) %>%
    ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(alpha = 0.6) +
    theme_cowplot(17) +
    xlab("Linz \u2192 Swiss (log2 fold-change)") +
    ylab("-log10 adjusted p-value") +
    ggtitle(paste("Volcano:", pop_name, "vs Linz"))
  print(volc_p)
  savePlot(volc_p, paste0("volcano_", pop_name, "VsLinz"))

  lfc_p <- res_df %>%
    filter(!is.na(padj), padj < padjCut) %>%
    ggplot(aes(x = cluster, y = log2FoldChange)) +
    geom_point() +
    theme_cowplot(17) +
    xlab("Cluster") +
    ylab("log2 fold-change (Swiss vs Linz)") +
    ggtitle(paste0("LFC: ", pop_name, " vs Linz (padj < ", padjCut, ")"))
  print(lfc_p)
  savePlot(lfc_p, paste0("LFC_", pop_name, "VsLinz"))
}

invisible(lapply(swissPops, make_pop_plots))


################ 5. SHARED CLUSTERS ACROSS POPULATIONS ------------------------

pop_res_all <- bind_rows(lapply(pop_out, `[[`, "res_df"))

## 5.1 significant in ALL three populations
results_combined <- pop_res_all %>%
  filter(!is.na(padj), padj < padjCut) %>%
  pivot_wider(id_cols = "cluster", names_from = "population", values_from = "log2FoldChange") %>%
  filter(if_all(all_of(swissPops), ~ !is.na(.x)))
write_tsv(results_combined, file.path(outDir, "results_combined_geo.tsv"))

## 5.2 significant in AT LEAST two populations
results_shared2 <- pop_res_all %>%
  filter(!is.na(padj), padj < padjCut) %>%
  group_by(cluster) %>%
  summarise(n_pops   = n(),
            pops     = paste(population, collapse = "/"),
            mean_lfc = mean(log2FoldChange),
            .groups  = "drop") %>%
  filter(n_pops >= 2)
write_tsv(results_shared2, file.path(outDir, "results_shared2_geo.tsv"))
table(results_shared2$pops)

## 5.3 Venn of significant clusters per population
venn_list <- lapply(pop_out, function(x) x$res_df$cluster[!is.na(x$res_df$padj) & x$res_df$padj < padjCut])
venn_p <- ggVennDiagram(venn_list) +
  ggtitle("Shared significant clusters across Swiss populations (vs Linz)")
print(venn_p)
savePlot(venn_p, "venn_sharedClusters")

## 5.4 combined LFC plot (clusters significant in all three)
lfc_comb <- results_combined %>%
  pivot_longer(cols = all_of(swissPops), names_to = "population", values_to = "log2FC") %>%
  ggplot(aes(x = cluster, y = log2FC, color = population)) +
  geom_point() +
  theme_cowplot(17) +
  ylab("log2 fold-change (Swiss vs Linz)") +
  xlab("Cluster") +
  theme(legend.title = element_blank())
print(lfc_comb)
savePlot(lfc_comb, "LFC_combined")


################ 6. SESSION INFO ----------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(outDir, "sessionInfo.txt"))
message("Done. All outputs in ", normalizePath(outDir))