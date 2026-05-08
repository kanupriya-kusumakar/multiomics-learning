library(SummarizedExperiment)
library(airway)
library(DESeq2)

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

