# Reproduce tutorial found on https://github.com/ytakemon/GRETA 
message("This is script reproduces the tutorial that identified ARID1A genetic interactors and coessential genes.")

# Load libraries 
library(tidyverse)
library(GRETA)

# Set paths 
GRETA_data_dir <- "/opt/GRETA/data/20Q1/data/"
GRETA_output_dir <- paste0(getwd(),"/output")

# Create output directory if not already created
if(!dir.exists(GRETA_output_dir)){
  dir.create(GRETA_output_dir)
}

# Identify cell lines with ARID1A WT allele or LOF mutations
message("Assigning cancer cell lines into ARID1A control group and HomDel mutant group.")
ARID1A_groups <- select_cell_lines(Input_gene = "ARID1A", data_dir = GRETA_data_dir)
ARID1A_mutant_IDs <- ARID1A_groups %>% filter(Group %in% c("ARID1A_HomDel")) %>% pull(DepMap_ID)
ARID1A_control_IDs <- ARID1A_groups %>% filter(Group %in% c("Control")) %>% pull(DepMap_ID)

# Run genetic interaction screen ----------------------------------------------------------
# Total number of cores on machine:
all_cores <-  parallel::detectCores()

# This can take several hours depending on number of lines/cores used. Best to run this overnight.
message("Running genetic interaction screen using ", all_cores - 1, " threads.")
message("Depending on the number of threads used, this may take a while...")
screen_results <- GI_screen(
  control_IDs = ARID1A_control_IDs, 
  mutant_IDs = ARID1A_mutant_IDs,
  core_num = all_cores - 1, # leave one core free, just in case.
  output_dir = GRETA_output_dir, # Will save your results here as well as in the variable
  data_dir = GRETA_data_dir,
  test = FALSE) # use TRUE to run a short test to make sure all will run overnight.

message("Plotting ranked genetic interaction screen.")
# Plot ranked GI candidates 
GI_plot <- plot_screen(
    result_df = screen_results, 
    label_genes = TRUE, 
    label_n = 3)

# save plot 
pdf("./GRETA/output/Tutorial_ARID1A_GI_ranked_plot.pdf", width = 6, height = 4)
print(GI_plot)
dev.off()

# Run essentiality network analysis ---------------------------------------------------
# Pearson correlation
message("Running essentiality analysis.")
coess_df <- coessential_map(
  Input_gene = "ARID1A", 
  core_num = all_cores - 1, # leave one core free, just in case.
  output_dir = GRETA_output_dir, # Will save your results here as well as in the variable
  data_dir = GRETA_data_dir,
  test = FALSE)

# Calculate inflection points of positive and negative curve using co-essential gene results.
coess_inflection_df <- get_inflection_points(coess_df)

# Combine and annotate data frame containing co-essential genes
coess_annotated_df <- annotate_coessential_df(coess_df, coess_inflection_df)

# Plot ranked co-/anti-essential candidates 
message("Plotting ranked candidate co-/anti-essential genes.")
essential_plot <- plot_coessential_genes(
  result_df = coess_annotated_df, 
  inflection_df = coess_inflection_df,
  label_genes = TRUE, # Should gene names be labeled?
  label_n = 3) # Number of genes to display from each end

# save plot 
pdf("./GRETA/output/Tutorial_ARID1A_essentiality_ranked_plot.pdf", width = 6, height = 4)
print(essential_plot)
dev.off()

message("Results are save in `./GRETA/output/`")