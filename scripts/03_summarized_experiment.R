# Day 5: Bioconductor & SummarizedExperiment
# Kanupriya — multiomics-learning

library(SummarizedExperiment)
library(airway)

# Load dataset
data(airway)

# ── Inspect the object ──────────────────────────────────────────────────────
dim(airway)           # 63677 genes × 8 samples
colData(airway)       # sample metadata: cell line, treatment
head(rowData(airway)) # gene annotations: symbol, coordinates
assay(airway, "counts")[1:6, 1:4]  # count matrix corner

# ── Subsetting ──────────────────────────────────────────────────────────────
# Treated samples only
airway_treated <- airway[ , airway$dex == "trt"]
dim(airway_treated)   # 63677 × 4

# First 1000 genes only
airway_small <- airway[1:1000, ]
dim(airway_small)     # 1000 × 8

# Both: 1000 genes, treated samples
airway_subset <- airway[1:1000, airway$dex == "trt"]
dim(airway_subset)    # 1000 × 4

# ── Quick summary stats ─────────────────────────────────────────────────────
# Total counts per sample
colSums(assay(airway, "counts"))

# How many genes have zero counts across all samples?
zero_genes <- sum(rowSums(assay(airway, "counts")) == 0)
zero_genes