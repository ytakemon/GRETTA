# Load library
library(tidyverse)
library(psych)

# Define paths
GRETTA_data_dir <- "./DepMap/23Q4/"
# Load data
load(paste0(GRETTA_data_dir, "/gene_effect.rda"))

# Pan cancer -------------------------------------------
# Format as matrix
gene_effect_wide <- gene_effect %>%
  pivot_wider(names_from = "GeneNameID", values_from = "Effect_score") %>%
  arrange(DepMap_ID)

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
# save data individually for GRETTA 
save(fit, file = paste0(GRETTA_data_dir, "pancan_coess_precomputed.rda"))