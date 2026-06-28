#!/usr/bin/env Rscript
# ════════════════════════════════════════════════════════════════
#  DESeq2 Analysis Script — 
# ════════════════════════════════════════════════════════════════

cat("\n")
cat("════════════════════════════════════════════════════════════\n")
cat("   DESeq2 Analysis — \n")
cat("════════════════════════════════════════════════════════════\n\n")

# ════════════════════════════════════════════════════════════════
# 0. SET R LIBRARY PATH FIRST
# ════════════════════════════════════════════════════════════════
cat("Step 0: Setting up R library path...\n")

user_lib <- file.path(Sys.getenv("HOME"), "R",
                      paste0("x86_64-pc-linux-gnu-library/", getRversion()))
dir.create(user_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

cat("  ✓ R library ready\n\n")

# ════════════════════════════════════════════════════════════════
# 1. ARGUMENT PARSING
# ════════════════════════════════════════════════════════════════
cat("Step 1: Parsing arguments...\n")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript deseq2_analysis.R <metadata_file> <output_dir>\n", call. = FALSE)
}

metadata_file <- args[1]
output_dir    <- args[2]

if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file, call. = FALSE)
}

cat("  Metadata file:", metadata_file, "\n")
cat("  Output dir   :", output_dir, "\n\n")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ════════════════════════════════════════════════════════════════
# 2. INSTALL & LOAD PACKAGES 
# ════════════════════════════════════════════════════════════════
cat("Step 2: Loading packages...\n")

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", quiet = TRUE)
}

# FIXED: Added tibble to the list
packages <- c("tximport", "readr", "DESeq2", "ggplot2", "ggrepel",
              "pheatmap", "RColorBrewer", "dplyr", "apeglm", "tidyr", "tibble")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing", pkg, "...\n")
    if (pkg %in% c("tximport", "DESeq2", "apeglm")) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
    } else {
      install.packages(pkg, quiet = TRUE)
    }
  }
}

suppressPackageStartupMessages({
  library(tximport)
  library(readr)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(apeglm)
  library(tidyr)
  library(tibble)  # FIXED: Added this
})

cat("  ✓ All packages loaded\n\n")

# ════════════════════════════════════════════════════════════════
# 3. READ METADATA
# ════════════════════════════════════════════════════════════════
cat("Step 3: Reading metadata...\n")

samples <- read.delim(metadata_file, stringsAsFactors = FALSE)

required_cols <- c("sample", "condition", "path")
missing_cols <- setdiff(required_cols, colnames(samples))
if (length(missing_cols) > 0) {
  stop("Missing columns in metadata: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

samples$condition <- factor(samples$condition)
samples$condition <- relevel(samples$condition, ref = levels(samples$condition)[1])

cat("  ✓ Read", nrow(samples), "samples\n")
cat("  Conditions:", paste(levels(samples$condition), collapse = ", "), "\n")
print(samples)
cat("\n")

# ════════════════════════════════════════════════════════════════
# 4. BUILD FILES LIST & VALIDATE
# ════════════════════════════════════════════════════════════════
cat("Step 4: Validating kallisto files...\n")

files <- samples$path
names(files) <- samples$sample

missing_files <- which(!file.exists(files))
if (length(missing_files) > 0) {
  cat("ERROR: Missing files:\n")
  print(files[missing_files])
  stop("Cannot proceed without all kallisto files", call. = FALSE)
}

cat("  ✓ All", length(files), "kallisto files found\n\n")

# ════════════════════════════════════════════════════════════════
# 5. CREATE tx2gene FROM FIRST KALLISTO FILE
# ════════════════════════════════════════════════════════════════
cat("Step 5: Creating tx2gene mapping from Kallisto data...\n")

first_file <- files[1]
cat("  Reading:", first_file, "\n")

abund_df <- read_tsv(first_file, show_col_types = FALSE)

if (!"target_id" %in% colnames(abund_df)) {
  stop("ERROR: 'target_id' column not found. Check Kallisto output format.", call. = FALSE)
}

cat("  Sample target_id entries (first 3):\n")
print(head(abund_df$target_id, 3))
cat("\n")

# ════════════════════════════════════════════════════════════════
# EXTRACT PIPE-DELIMITED KALLISTO FORMAT
# ════════════════════════════════════════════════════════════════
extract_field <- function(x, field_num = 1, sep = "|") {
  parts <- strsplit(x, sep, fixed = TRUE)[[1]]
  field <- parts[field_num]
  field <- sub("\\.[0-9]+$", "", field)
  return(field)
}

cat("  Parsing pipe-delimited transcript/gene IDs...\n")

tx_ids <- sapply(abund_df$target_id,
                 function(x) extract_field(x, field_num = 1, sep = "|"),
                 USE.NAMES = FALSE)

gene_ids <- sapply(abund_df$target_id,
                   function(x) extract_field(x, field_num = 2, sep = "|"),
                   USE.NAMES = FALSE)

tx2gene <- data.frame(
  TXNAME = tx_ids,
  GENEID = gene_ids,
  stringsAsFactors = FALSE
)

cat("  ✓ Parsed transcript and gene IDs\n")
cat("    Unique transcripts:", nrow(tx2gene), "\n")
cat("    Unique genes:", n_distinct(tx2gene$GENEID), "\n")
cat("    Sample mappings (first 3):\n")
print(head(tx2gene, 3))
cat("\n")

# ════════════════════════════════════════════════════════════════
# 6. IMPORT WITH TXIMPORT
# ════════════════════════════════════════════════════════════════
cat("Step 6: Importing kallisto counts with tximport...\n")

txi <- tximport(
  files,
  type = "kallisto",
  tx2gene = tx2gene,
  txOut = FALSE,
  ignoreTxVersion = TRUE,
  importer = function(x) read_tsv(x, show_col_types = FALSE)
)

cat("  ✓ Imported", nrow(txi$counts), "genes from", ncol(txi$counts), "samples\n")
cat("    Counts matrix dimensions:", nrow(txi$counts), "x", ncol(txi$counts), "\n\n")

# ════════════════════════════════════════════════════════════════
# 7. CREATE DESeq2 DATASET
# ════════════════════════════════════════════════════════════════
cat("Step 7: Creating DESeq2 dataset...\n")

dds <- DESeqDataSetFromTximport(
  txi,
  colData = samples,
  design = ~ condition
)

cat("  ✓ DESeq2 object created\n")
cat("    Genes:", nrow(dds), "\n")
cat("    Samples:", ncol(dds), "\n\n")

# ════════════════════════════════════════════════════════════════
# 8. PRE-FILTER LOW COUNTS 
# ════════════════════════════════════════════════════════════════
cat("Step 8: Pre-filtering low-count genes...\n")

keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

cat("  ✓ Filtered to", nrow(dds), "genes (removed", sum(!keep), "low-count genes)\n\n")

# ════════════════════════════════════════════════════════════════
# 9. RUN DESeq2 PIPELINE
# ════════════════════════════════════════════════════════════════
cat("Step 9: Running DESeq2 (estimating size factors, dispersion, testing)...\n")

dds <- DESeq(dds, quiet = FALSE)

cat("  ✓ DESeq2 pipeline complete\n\n")

# ════════════════════════════════════════════════════════════════
# 10. EXTRACT & SHRINK RESULTS
# ════════════════════════════════════════════════════════════════
cat("Step 10: Extracting and shrinking log fold changes...\n")

res <- results(dds, alpha = 0.05)

resultsNames(dds)
coef_name <- resultsNames(dds)[grep("condition", resultsNames(dds))][1]

cat("  Using coefficient:", coef_name, "\n")

res <- lfcShrink(
  dds,
  coef = coef_name,
  type = "apeglm",
  quiet = TRUE
)

cat("  ✓ LFC shrinkage complete\n\n")

# ════════════════════════════════════════════════════════════════
# 11. CONVERT TO DATA FRAME 
# ════════════════════════════════════════════════════════════════
cat("Step 11: Preparing results tables...\n")

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj) %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange >= 1 ~ "Up",
      padj < 0.05 & log2FoldChange <= -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

cat("  Significantly DE genes:\n")
cat("    Up-regulated  :", sum(res_df$significance == "Up", na.rm = TRUE), "\n")
cat("    Down-regulated:", sum(res_df$significance == "Down", na.rm = TRUE), "\n")
cat("    Not significant:", sum(res_df$significance == "NS", na.rm = TRUE), "\n\n")

# ════════════════════════════════════════════════════════════════
# 12. EXPORT TABLES
# ════════════════════════════════════════════════════════════════
cat("Step 12: Exporting CSV and text tables...\n")

write.csv(res_df,
          file.path(output_dir, "all_genes_DESeq2_results.csv"),
          row.names = FALSE)
cat("  ✓ all_genes_DESeq2_results.csv\n")

sig_df <- res_df %>% filter(padj < 0.05)
write.csv(sig_df,
          file.path(output_dir, "significant_genes.csv"),
          row.names = FALSE)
cat("  ✓ significant_genes.csv (", nrow(sig_df), " genes)\n", sep = "")

up_df <- res_df %>% filter(padj < 0.05 & log2FoldChange > 1)
write.csv(up_df,
          file.path(output_dir, "upregulated_genes.csv"),
          row.names = FALSE)
cat("  ✓ upregulated_genes.csv (", nrow(up_df), " genes)\n", sep = "")

down_df <- res_df %>% filter(padj < 0.05 & log2FoldChange < -1)
write.csv(down_df,
          file.path(output_dir, "downregulated_genes.csv"),
          row.names = FALSE)
cat("  ✓ downregulated_genes.csv (", nrow(down_df), " genes)\n", sep = "")

write.table(sig_df$gene_id,
            file.path(output_dir, "significant_gene_list.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("  ✓ significant_gene_list.txt\n")

write.table(up_df$gene_id,
            file.path(output_dir, "upregulated_gene_list.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("  ✓ upregulated_gene_list.txt\n")

write.table(down_df$gene_id,
            file.path(output_dir, "downregulated_gene_list.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("  ✓ downregulated_gene_list.txt\n\n")

# ════════════════════════════════════════════════════════════════
# 13. VARIANCE TRANSFORMATION
# ════════════════════════════════════════════════════════════════
cat("Step 13: Variance stabilizing transformation for plots...\n")

vsd <- vst(dds, blind = TRUE)

cat("  ✓ VST complete\n\n")

# ════════════════════════════════════════════════════════════════
# 14. PCA PLOT - GENE-LEVEL 
# ════════════════════════════════════════════════════════════════
cat("Step 14: Generating gene-level PCA plot...\n")

# Select top 500 variable genes for PCA
gene_var <- apply(assay(vsd), 1, var)
top_var_genes <- head(order(gene_var, decreasing = TRUE), 500)
gene_mat <- assay(vsd)[top_var_genes, ]

# Perform PCA on genes (genes as observations, samples as features)
pca_result <- prcomp(gene_mat, scale. = TRUE)
pca_data_genes <- as.data.frame(pca_result$x)
pca_data_genes$gene_id <- rownames(pca_data_genes)

# Join with significance information from results
pca_data_genes <- pca_data_genes %>%
  left_join(res_df[, c("gene_id", "significance")], by = "gene_id") %>%
  mutate(significance = ifelse(is.na(significance), "NS", significance))

percentVar <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)

cols <- c(Up = "#E63946", Down = "#457B9D", NS = "#CCCCCC")

pca_plot <- ggplot(pca_data_genes, aes(x = PC1, y = PC2, colour = significance)) +
  geom_point(size = 2, alpha = 0.6) +
  scale_colour_manual(values = cols) +
  labs(title = "PCA Plot - Genes (Top 500 variable)",
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       colour = "Significance") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

png(file.path(output_dir, "PCA_plot.png"), width = 1200, height = 1000, res = 100)
print(pca_plot)
dev.off()

cat("  ✓ PCA_plot.png (gene-level, top 500 variable genes)\n\n")

# ════════════════════════════════════════════════════════════════
# 15. SAMPLE DISTANCE HEATMAP
# ════════════════════════════════════════════════════════════════
cat("Step 15: Generating sample distance heatmap...\n")

sampleDists <- dist(t(assay(vsd)))

png(file.path(output_dir, "sample_distance_heatmap.png"),
    width = 1200, height = 1000, res = 100)
pheatmap(as.matrix(sampleDists),
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         color = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),
         main = "Sample-to-Sample Distance")
dev.off()

cat("  ✓ sample_distance_heatmap.png\n\n")

# ════════════════════════════════════════════════════════════════
# 16. MA PLOT
# ════════════════════════════════════════════════════════════════
cat("Step 16: Generating MA plot...\n")

png(file.path(output_dir, "MA_plot.png"), width = 1200, height = 1000, res = 100)
plotMA(res, ylim = c(-5, 5))
dev.off()

cat("  ✓ MA_plot.png\n\n")

# ════════════════════════════════════════════════════════════════
# 17. DISPERSION PLOT
# ════════════════════════════════════════════════════════════════
cat("Step 17: Generating dispersion plot...\n")

png(file.path(output_dir, "dispersion_plot.png"),
    width = 1200, height = 1000, res = 100)
plotDispEsts(dds)
dev.off()

cat("  ✓ dispersion_plot.png\n\n")

# ════════════════════════════════════════════════════════════════
# 18.  VOLCANO PLOT 
# ════════════════════════════════════════════════════════════════
cat("Step 18: Generating  volcano plot...\n")

cols <- c(Up = "#E63946", Down = "#457B9D", NS = "#CCCCCC")

volcano_plot <- ggplot(res_df,
                       aes(x = log2FoldChange, y = -log10(pvalue), colour = significance)) +
  geom_point(size = 2, alpha = 0.7, stroke = 0) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  scale_colour_manual(values = cols, name = "Significance") +
  labs(title = "Volcano Plot - Differential Expression Analysis",
       x = expression(log[2]~Fold~Change),
       y = expression(-log[10]~(p-value)),
       colour = "Significance") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5, margin = margin(b = 10)),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 11),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )

# Add top significant genes labels - improved spacing
top_sig <- res_df %>% 
  filter(significance != "NS") %>% 
  arrange(padj) %>% 
  slice_head(n = 10)

if (nrow(top_sig) > 0) {
  volcano_plot <- volcano_plot +
    geom_text_repel(
      data = top_sig, 
      aes(label = gene_id),
      size = 3.5, 
      colour = "black",
      fontface = "bold",
      box.padding = 0.5,
      point.padding = 0.3,
      segment.size = 0.4,
      max.overlaps = Inf
    )
}

png(file.path(output_dir, "Volcano_plot.png"), width = 1400, height = 1200, res = 100)
print(volcano_plot)
dev.off()

cat("  ✓ Volcano_plot.png (enhanced with gene labels)\n\n")

# ════════════════════════════════════════════════════════════════
# 19.  TOP 50 GENES HEATMAP - 
# ════════════════════════════════════════════════════════════════
cat("Step 19: Generating  top 50 DE genes heatmap...\n")

topgenes <- head(order(res$padj), min(50, nrow(res)))
mat <- assay(vsd)[topgenes, ]

png(file.path(output_dir, "Top50_heatmap.png"), width = 1300, height = 1400, res = 100)
pheatmap(
  mat, 
  scale = "row",
  main = "Top 50 Differentially Expressed Genes\n(VST, row-scaled)",
  color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(256),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = ifelse(nrow(mat) > 50, 7, 9),
  fontsize_col = 10,
  fontsize = 12,
  border_color = "white",
  cellwidth = 35,
  cellheight = ifelse(nrow(mat) > 50, 8, 12),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_method = "ward.D2",
  treeheight_row = 20,
  treeheight_col = 20
)
dev.off()

cat("  ✓ Top50_heatmap.png\n\n")

# ════════════════════════════════════════════════════════════════
# 20.  SCATTER PLOT - Gene Expression Correlation
# ════════════════════════════════════════════════════════════════
cat("Step 20: Generating gene expression scatter plot...\n")

# Extract normalized counts by condition - FIXED VERSION
norm_counts_matrix <- counts(dds, normalized = TRUE)
cond_levels <- levels(samples$condition)

# Calculate mean counts per condition
if (length(cond_levels) >= 2) {
  cond1_samples <- samples$sample[samples$condition == cond_levels[1]]
  cond2_samples <- samples$sample[samples$condition == cond_levels[2]]
  
  cond1_means <- rowMeans(norm_counts_matrix[, cond1_samples, drop = FALSE]) + 0.5
  cond2_means <- rowMeans(norm_counts_matrix[, cond2_samples, drop = FALSE]) + 0.5
  
  scatter_data_plot <- data.frame(
    gene_id = rownames(norm_counts_matrix),
    cond1_expr = cond1_means,
    cond2_expr = cond2_means,
    stringsAsFactors = FALSE
  )
  
  # Add significance information
  scatter_data_plot <- scatter_data_plot %>%
    left_join(res_df[, c("gene_id", "significance")], by = "gene_id") %>%
    mutate(significance = ifelse(is.na(significance), "NS", significance))
  
  scatter_plot <- ggplot(scatter_data_plot, 
                         aes(x = cond1_expr, y = cond2_expr, 
                             color = significance)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_colour_manual(values = cols, name = "Significance") +
    scale_x_log10() +
    scale_y_log10() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                color = "grey50", linewidth = 0.7) +
    labs(
      title = "Gene Expression Scatter Plot",
      x = paste(cond_levels[1], "- log10(Mean Normalized Counts)"),
      y = paste(cond_levels[2], "- log10(Mean Normalized Counts)"),
      color = "Significance"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      legend.position = "top",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
    )
  
  png(file.path(output_dir, "Gene_expression_scatter.png"), 
      width = 1200, height = 1100, res = 100)
  print(scatter_plot)
  dev.off()
  
  cat("  ✓ Gene_expression_scatter.png\n\n")
} else {
  cat("  ⚠ Cannot create scatter plot: fewer than 2 conditions\n\n")
}

# ════════════════════════════════════════════════════════════════
# 21.  HISTOGRAM - Gene Expression Distribution
# ════════════════════════════════════════════════════════════════
cat("Step 21: Generating gene expression distribution histogram...\n")

# Get all normalized counts and prepare for histogram - FIXED
norm_counts_matrix <- counts(dds, normalized = TRUE)
all_counts <- as.numeric(norm_counts_matrix)
log10_counts <- log10(all_counts + 0.5)

norm_counts_long <- data.frame(log10_counts = log10_counts)

histogram_plot <- ggplot(norm_counts_long, aes(x = log10_counts)) +
  geom_histogram(bins = 50, fill = "#3498db", alpha = 0.7, color = "black", linewidth = 0.4) +
  labs(
    title = "Gene Expression Distribution",
    x = "log10(Normalized Counts + 0.5)",
    y = "Frequency"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )

png(file.path(output_dir, "Gene_expression_histogram.png"), 
    width = 1200, height = 900, res = 100)
print(histogram_plot)
dev.off()

cat("  ✓ Gene_expression_histogram.png\n\n")

# ════════════════════════════════════════════════════════════════
# 22.  BOXPLOT - Gene Expression by Sample/Group
# ════════════════════════════════════════════════════════════════
cat("Step 22: Generating gene expression boxplot by group...\n")

# Prepare data for boxplot - VST counts by condition - OPTIMIZED
vst_matrix <- assay(vsd)

# Create long-format data frame efficiently
vst_counts <- data.frame(
  sample = rep(colnames(vst_matrix), each = nrow(vst_matrix)),
  vst_count = as.numeric(vst_matrix),
  stringsAsFactors = FALSE
)

# Merge with sample metadata
vst_counts <- vst_counts %>%
  left_join(samples[, c("sample", "condition")], by = "sample")

# Ensure proper factor ordering for plot
vst_counts$sample <- factor(vst_counts$sample, levels = samples$sample)
vst_counts$condition <- factor(vst_counts$condition)

boxplot <- ggplot(vst_counts, aes(x = sample, y = vst_count, fill = condition)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1.5, outlier.alpha = 0.6) +
  scale_fill_brewer(palette = "Set2", name = "Condition") +
  labs(
    title = "Gene Expression Distribution by Sample",
    x = "Sample",
    y = "VST-Transformed Counts",
    fill = "Condition"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )

png(file.path(output_dir, "Gene_expression_boxplot.png"), 
    width = 1300, height = 1000, res = 100)
print(boxplot)
dev.off()

cat("  ✓ Gene_expression_boxplot.png\n\n")

# ════════════════════════════════════════════════════════════════
# 23. SUMMARY STATISTICS
# ════════════════════════════════════════════════════════════════
cat("Step 23: Saving summary statistics...\n")

summary_stats <- data.frame(
  Metric = c("Total genes quantified", "Genes after filtering",
             "Significantly DE genes (padj < 0.05)",
             "Up-regulated (LFC > 1)", "Down-regulated (LFC < -1)",
             "Analysis date"),
  Value = c(nrow(res), nrow(dds), nrow(sig_df), nrow(up_df), nrow(down_df),
            as.character(Sys.Date()))
)

write.csv(summary_stats, file.path(output_dir, "summary_statistics.csv"),
          row.names = FALSE)

cat("  ✓ summary_statistics.csv\n")
cat("  Analysis summary:\n")
print(summary_stats)
cat("\n")

# ════════════════════════════════════════════════════════════════
# 24. SESSION INFO
# ════════════════════════════════════════════════════════════════
cat("Step 24: Saving session information...\n")

sink(file.path(output_dir, "session_info.txt"))
print(sessionInfo())
sink()

cat("  ✓ session_info.txt\n\n")

# ════════════════════════════════════════════════════════════════
# COMPLETION MESSAGE
# ════════════════════════════════════════════════════════════════
cat("════════════════════════════════════════════════════════════\n")
cat("   ✓ DESeq2 ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("════════════════════════════════════════════════════════════\n")
cat("\nOutput saved to:", output_dir, "\n\n")
cat("Generated files:\n")
cat("  CSV Tables:\n")
cat("    · all_genes_DESeq2_results.csv\n")
cat("    · significant_genes.csv\n")
cat("    · upregulated_genes.csv\n")
cat("    · downregulated_genes.csv\n")
cat("    · summary_statistics.csv\n\n")
cat("  Gene Lists (txt):\n")
cat("    · significant_gene_list.txt\n")
cat("    · upregulated_gene_list.txt\n")
cat("    · downregulated_gene_list.txt\n\n")
cat("  Standard Plots (PNG):\n")
cat("    · PCA_plot.png (gene-level)\n")
cat("    · sample_distance_heatmap.png\n")
cat("    · MA_plot.png\n")
cat("    · dispersion_plot.png\n\n")
cat("  Plots generated (PNG):\n")
cat("    · Volcano_plot.png\n")
cat("    · Top50_heatmap.png\n")
cat("    · Gene_expression_scatter.png\n")
cat("    · Gene_expression_histogram.png\n")
cat("    · Gene_expression_boxplot.png\n\n")
cat("════════════════════════════════════════════════════════════\n\n")
