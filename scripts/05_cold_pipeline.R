#---1. Install Packages and Load Libraries
library(tidyverse)
library(SummarizedExperiment)
library(ggplot2)
library(EnhancedVolcano)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DESeq2)
library(airway)
# ---2. Load 
data(airway) 

#---3. Filter the data
airway_filtered <- airway[rowSums(assay(airway, "counts")) >=10, ]
# Set reference level BEFORE building the model
airway_filtered$dex <- relevel(airway_filtered$dex, ref = "untrt")

#---4. Build DESeq2 dataset
dds <- DESeqDataSet(airway_filtered, design = ~ cell + dex)
dds

#---5.Run DESeq2 
dds <- DESeq(dds)
#--6. Extract results
res <- results(dds, contrast = c("dex", "trt", "untrt"))
summary(res)

#---7. lfc Shrink
library(apeglm)
resultsNames(dds) 
res_shrunk <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm")
summary(res_shrunk)

#---8. Add gene symbols
res_shrunk$symbol <- rowData(airway_filtered)$symbol[ match(rownames(res_shrunk), rownames(airway_filtered)) ]

#---9. Filter Significant DEGs ---
sig_genes <- res_shrunk [!is.na(res_shrunk$padj) & 
                           res_shrunk$padj <0.05 & 
                           abs(res_shrunk$log2FoldChange) >1, ]

#---10. Volcano Plot 
EnhancedVolcano(res_shrunk,
                lab = res_shrunk$symbol,
                selectLab = sig_genes$symbol,   # <-- only these get labels
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.05,
                FCcutoff = 1,
                drawConnectors = TRUE,         
                title = "Dexamethasone effect — airway smooth muscle",
                subtitle = "Treated vs Untreated",
                pointSize = 1.5,
                labSize = 3)

ggsave("outputs/figures/day8_volcano_shrunk_siggene.png", width = 10, height = 8, dpi = 150)

