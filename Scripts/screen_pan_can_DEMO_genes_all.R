library(pacman)
p_load(tidyverse)
DepMap_dir <- "/projects/marralab/ytakemon_prj/DepMap/20Q1/"

# get list ------------------------------------------------------------------
gene_list <- c("ARID1A", "ARID1B")

# Generate batch submission for cluster -------------------------------------
create_cluster_batch(N_batches = length(gene_list), gene_list = gene_list, file = "screen_pan_can_DEMO_genes_all.sh",
  path = "~/GitHub/PhDwork/Genetic_interaction_projects/20Q1/COMPASS/KMT2D/Reverse_screen/",
  shell_script = "screen_pan_can_DEMO_genes_one.sh")
