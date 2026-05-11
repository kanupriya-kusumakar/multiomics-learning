library(SummarizedExperiment)
library(airway)
library(DESeq2)
BiocManager::install("AnnotationDbi")
library(AnnotationDbi)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

data(airway)

# Filter: keep genes with at least 10 reads across all samples combined
keep <- rowSums(assay(airway, "counts")) >= 10
airway_filtered <- airway[keep, ]

dim(airway_filtered) 

airway_removed <- sum(rowSums(assay(airway, "counts")) <10)

# Build the DESeqDataSet — specify which column is your experimental variable
dds <- DESeqDataSet(airway_filtered, design = ~dex)
dds

# coment - dds <- DESeqDataSet(airway_filtered, design = ~dex, condition = "untrt", "trt")

# Run the full DESeq2 pipeline (normalisation + dispersion + testing)
dds <- DESeq(dds)
# Extract results — by default: trt vs untrt
res <- results(dds)
res

summary(res)

# IMPORTANT: comparison is untrt vs trt
# Positive LFC = higher in untreated = SUPPRESSED by dexamethasone
# Negative LFC = higher in treated = INDUCED by dexamethasone

res <- results(dds, contrast = c("dex", "trt", "untrt"))
summary(res)

# Sort by adjusted p-value, view top results
res_sorted <- res[order(res$padj), ]
head(res_sorted, 10)

sig_genes <- res_sorted[!is.na(res_sorted$padj) &
                          res_sorted$padj < 0.05 &
                          abs(res_sorted$log2FoldChange) > 1, ]

nrow(sig_genes)

rownames(res_sorted)[1]

# Install if needed:
BiocManager::install("EnhancedVolcano")

library(EnhancedVolcano)

EnhancedVolcano(res,
                lab = rownames(res),
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "Dexamethasone effect — airway smooth muscle",
                subtitle = "trt vs untrt",
                pointSize = 2,
                labSize = 3)
gene_names <- rowData(airway_filtered)$symbol
res$symbol <- gene_names[match(rownames(res), rownames(airway_filtered))]
head(as.data.frame(res[, c("log2FoldChange", "padj", "symbol")]))

EnhancedVolcano(res,
                lab = res$symbol,          # <-- changed from rownames(res)
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "Dexamethasone effect — airway smooth muscle",
                subtitle = "trt vs untrt",
                pointSize = 2,
                labSize = 3)

ggsave("outputs/figures/day6_volcano.png", width = 10, height = 8, dpi = 150)

sig_genes$symbol <- gene_names[match(rownames(sig_genes), rownames(airway_filtered))]

# Plot all genes, label only significant ones
EnhancedVolcano(res,
                lab = res$symbol,
                selectLab = sig_genes$symbol,   # <-- only these get labels
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.05,
                FCcutoff = 1,
                drawConnectors = TRUE,          # lines from label to point
                title = "Dexamethasone effect — airway smooth muscle",
                subtitle = "trt vs untrt",
                pointSize = 2,
                labSize = 3)
ggsave("outputs/figures/day6_volcano_siglabel.png", width = 10, height = 8, dpi = 150)

# ── Day 7: LFC Shrinkage ──────────────────────────────────────
BiocManager::install("apeglm")
library(apeglm)

# Shrink LFC estimates for the trt vs untrt coefficient
# coef name must match what DESeq2 calls it internally
resultsNames(dds)   # run this first to see the exact name

# Shrink LFC estimates using apeglm method
res_shrunk <- lfcShrink(dds, coef = "dex_untrt_vs_trt", type = "apeglm")
summary(res_shrunk)

# Compare raw vs shrunk LFC for the same genes
comparison <- data.frame(
  raw_LFC    = res$log2FoldChange,
  shrunk_LFC = res_shrunk$log2FoldChange
)

# Look at the most extreme raw estimates
comparison[order(abs(comparison$raw_LFC), decreasing = TRUE), ] |> head(10)

# Scatter: raw vs shrunk (sign-flipped to match direction)
plot(comparison$raw_LFC, -comparison$shrunk_LFC,
     pch = 20, cex = 0.3,
     col = rgb(0, 0, 0, 0.1),
     xlab = "Raw LFC (trt vs untrt)",
     ylab = "Shrunk LFC (trt vs untrt)",
     main = "Effect of lfcShrink on fold change estimates")
abline(0, 1, col = "red", lty = 2)   # y = x line

# Volcano with shrunk LFC (flip sign to get trt vs untrt direction)
res_shrunk_df <- as.data.frame(res_shrunk)
res_shrunk_df$log2FoldChange <- -res_shrunk_df$log2FoldChange  # flip direction

# Add gene symbols
res_shrunk_df$symbol <- mapIds(org.Hs.eg.db,
                               keys = rownames(res_shrunk_df),
                               column = "SYMBOL",
                               keytype = "ENSEMBL",
                               multiVals = "first")

EnhancedVolcano(res_shrunk_df,
                lab = res_shrunk_df$symbol,
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.05,
                FCcutoff = 1,
                selectLab = res_shrunk_df$symbol[which(res_shrunk_df$padj < 0.05 &
                                                         abs(res_shrunk_df$log2FoldChange) > 1)],
                title = "Dex vs Untreated (shrunk LFC)",
                subtitle = "apeglm shrinkage")

ggsave("outputs/figures/day7_volcano_shrunk.png", width = 10, height = 8, dpi = 150)

# ── Day 7: VST + PCA ─────────────────────────────────────────
# Variance-stabilise counts for sample-level exploration
vsd <- vst(dds, blind = TRUE)
vsd

# PCA coloured by treatment
plotPCA(vsd, intgroup = "dex")
# Colour by cell line to see the batch structure
plotPCA(vsd, intgroup = "cell")

# Save both PCA plots
png("outputs/figures/day7_pca_treatment.png", width = 800, height = 600, res = 120)
plotPCA(vsd, intgroup = "dex")
dev.off()

png("outputs/figures/day7_pca_cellline.png", width = 800, height = 600, res = 120)
plotPCA(vsd, intgroup = "cell")
dev.off()