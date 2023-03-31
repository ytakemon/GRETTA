# Load library
library(tidyverse)
library(psych)

# Define paths
greta_data_dir <- "/projects/marralab/ytakemon_prj/DepMap/GRETA_data/20Q1/data/"
greta_output_dir <- "/projects/marralab/ytakemon_prj/DepMap/GRETA_troubleshooting/"

# Load data
load(paste0(greta_data_dir, "/gene_effect.rda"))

# Pan cancer -------------------------------------------
# Format as matrix
gene_effect_wide <- gene_effect

gene_effect_wide_mat <- gene_effect_wide %>%
  select(-DepMap_ID) %>%
  as.matrix()
rownames(gene_effect_wide_mat) <- gene_effect_wide$DepMap_ID

fit <- corr.test(
    gene_effect_wide_mat, 
    method = "pearson",
    adjust = "BH", 
    ci = FALSE)

#### SAVE pre-computed data for GRETA ! ----------------------------------------------------------
GRETA_dir_20Q1 <- "/projects/marralab/ytakemon_prj/DepMap/GRETA_data/20Q1/data/"

# save data individually for GRETA 
save(fit, file = paste0(GRETA_dir_20Q1, "pancan_coess_precomputed.rda"))