# Reproduce tutorial found on https://github.com/ytakemon/GRETTA 
message("This is script reproduces the tutorial that identified ARID1A genetic interactors and coessential genes.")
num_threads <- as.numeric(commandArgs(trailingOnly = TRUE))[1]
data_dir <- as.character(commandArgs(trailingOnly = TRUE))[2]

message(commandArgs(trailingOnly = TRUE))
message(num_threads)
message(data_dir)
# Load libraries 
library(tidyverse)
library(GRETTA)

# Set paths 
gretta_data_dir <- data_dir
if(!dir.exists(gretta_data_dir)){
  stop(gretta_data_dir, "Doesn't exist.")
}

if(!dir.exists(paste0(getwd(),"/output/"))){
  dir.create(paste0(getwd(),"/output/"))
}
gretta_output_dir <- paste0(getwd(),"/output/")

# Create output directory if not already created
if(!dir.exists(gretta_output_dir)){
  dir.create(gretta_output_dir)
}

# Identify cell lines with ARID1A WT allele or LOF mutations
message("Assigning cancer cell lines into ARID1A control group and HomDel mutant group.")
ARID1A_groups <- select_cell_lines(input_gene = "ARID1A", data_dir = gretta_data_dir)
ARID1A_mutant_id <- ARID1A_groups %>% filter(Group %in% c("ARID1A_HomDel")) %>% pull(DepMap_ID)
ARID1A_control_id <- ARID1A_groups %>% filter(Group %in% c("Control")) %>% pull(DepMap_ID)

# Run genetic interaction screen ----------------------------------------------------------
# This can take several hours depending on number of lines/cores used. Best to run this overnight.
message("Running genetic interaction screen using ", num_threads, " threads.")
message("Depending on the number of threads used, this may take a while...")
screen_results <- GI_screen(
  control_id = ARID1A_control_id, 
  mutant_id = ARID1A_mutant_id,
  core_num = num_threads, # depends on how many cores you have  
  output_dir = gretta_output_dir, # Will save your results here as well as in the variable
  data_dir = gretta_data_dir,
  test = FALSE)
  
message("Plotting ranked genetic interaction screen.")
# Plot ranked GI candidates 

# save plot 
cairo_pdf(paste0(gretta_output_dir, "Tutorial_ARID1A_GI_ranked_plot.pdf"), width = 6, height = 4)
plot_screen(
    result_df = screen_results, 
    label_genes = TRUE, 
    label_n = 3)
dev.off()

# Run essentiality network analysis ---------------------------------------------------
# Pearson correlation
message("Running essentiality analysis.")
coess_df <- coessential_map(
  input_gene = "ARID1A", 
  data_dir = gretta_data_dir, 
  output_dir = gretta_output_dir) 

# Calculate inflection points of positive and negative curve using co-essential gene results.
coess_inflection_df <- get_inflection_points(coess_df)

# Combine and annotate data frame containing co-essential genes
coess_annotated_df <- annotate_df(coess_df, coess_inflection_df)

# Plot ranked co-/anti-essential candidates 
message("Plotting ranked candidate co-/anti-essential genes.")

# save plot 
cairo_pdf(paste0(gretta_output_dir, "Tutorial_ARID1A_essentiality_ranked_plot.pdf"), width = 6, height = 4)
plot_coess(
  result_df = coess_annotated_df, 
  inflection_df = coess_inflection_df,
  label_genes = TRUE, # Should gene names be labeled?
  label_n = 3)
dev.off()

message("Results are save in: ", gretta_output_dir)