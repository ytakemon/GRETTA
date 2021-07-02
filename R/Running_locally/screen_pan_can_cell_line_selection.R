# To render rmarkdown to html report for cell line selection:
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rmarkdown)
Local_dir <- "/where/ever/you/decide/"
# This might only be a thing for me....
# Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio-server/bin/pandoc")

# Example for only 1 gene
gene_group <- c("ARID1A")

for(gene in gene_group){
  print(gene)

  render(
    input = paste0(Local_dir,"Script/20Q1_select_cell_lines.Rmd"),
    output_format = 'html_document',
    output_file = paste0(Local_dir,"/Results/Pan_Cancer",gene,"/20Q1_select_cell_lines.html"),
    clean = TRUE,
    params = list(Target_gene = gene),
    quiet = FALSE)
}
