# =============================================================
# Day 3: ggplot2 for publication-quality figures
# Dataset: palmerpenguins + simulated DEG data
# =============================================================

library(tidyverse)
library(palmerpenguins)

# install.packages("patchwork")  # uncomment and run once if not installed
library(patchwork)

# --- Plot 1: Bill morphology scatterplot ---
p1 <- ggplot(penguins, aes(bill_length_mm, bill_depth_mm, colour = species)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  theme_classic() +
  scale_colour_manual(values = c("#378ADD", "#1D9E75", "#D85A30")) +
  labs(
    title  = "Bill morphology by species",
    x      = "Bill length (mm)",
    y      = "Bill depth (mm)",
    colour = "Species"
  )

p1 

p2 <- ggplot(penguins, aes(species, body_mass_g, fill = species)) +
  geom_boxplot(alpha= 0.4, outlier.shape = NA)+
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
  theme_classic() +
  theme(legend.position ="none") +
  scale_fill_manual(values = c("#378ADD", "#1D9E75", "#D85A30")) +
  labs( 
    title = "Species Body Mass" ,
    x = "Species" ,
    y = "Body Mass (g)", 
    )
p2

# Step 1 — simulate the data (copy this exactly, don't change it)
set.seed(42)
df <- data.frame(
  log2FC = rnorm(2000, 0, 2),
  pval   = runif(2000)
)
df$neglog10p <- -log10(df$pval)
df$sig       <- abs(df$log2FC) > 1 & df$pval < 0.05

p3 <- ggplot(df, aes(log2FC, neglog10p, colour = sig)) +
  geom_point(alpha = 0.5, size = 0.8) +
  theme_classic() +
  scale_color_manual(values = c("grey70", "#E24B4A")) +
  labs(
    title = "Significant Genes", 
    x= "Log2 fold change", 
    y= "-log10(p-value)", 
  )
p3

# Combine plots into a multipanel figure
combined <- (p1 | p2) / p3

combined

# Save individual plots
ggsave(
  filename = "outputs/figures/p1_scatterplot.png",
  plot     = p1,
  width    = 6,
  height   = 4,
  dpi      = 300
)

ggsave(
  filename = "outputs/figures/p2_boxplot.png",
  plot     = p2,
  width    = 5,
  height   = 4,
  dpi      = 300
)

ggsave(
  filename = "outputs/figures/p3_volcano.png",
  plot     = p3,
  width    = 6,
  height   = 5,
  dpi      = 300
)

# Save the combined figure
ggsave(
  filename = "outputs/figures/combined_day3.png",
  plot     = combined,
  width    = 12,
  height   = 8,
  dpi      = 300
)
