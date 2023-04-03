# Load library
library(tidyverse)
library(psych)

# Define paths
GRETTA_data_dir <- "/projects/marralab/ytakemon_prj/DepMap/GRETTA_data/20Q1/data/"
GRETTA_output_dir <- "/projects/marralab/ytakemon_prj/DepMap/GRETTA_troubleshooting/"

# Load data
load(paste0(GRETTA_data_dir, "/gene_effect.rda"))

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

#### SAVE pre-computed data for GRETTA ! ----------------------------------------------------------
GRETTA_dir_20Q1 <- "/projects/marralab/ytakemon_prj/DepMap/GRETTA_data/20Q1/data/"

# save data individually for GRETTA 
save(fit, file = paste0(GRETTA_dir_20Q1, "pancan_coess_precomputed.rda"))