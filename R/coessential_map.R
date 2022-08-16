#' @title Perform co-essentially mapping
#' 
#' @description Performs multiple correltion coefficent analysese and determines cut to identify most likely co-essential genes.
#' 
#' @param Input_gene string, A vector containing one Hugo Symbol, Default: NULL
#' @param Input_disease string, A vector one or more disease contexts, Will perform pan-cancer analyses 
#' (all cell lines) by default, Default: NULL
#' @param Input_cell_lines string, A vector DepMap_IDs for which co-essentiality mapping will be performed on. 
#' Will perform pan-cancer analyses (all cell lines) by default, Default: NULL
#' @param core_num integer, Number of cores to run analysis, Default: NULL
#' @param output_dir string, Full path to where output file should be saved, Default: NULL
#' @param data_dir string Path to GINIR_data
#' @param test logical, TRUE/FALSE whether you want to run only a small subset (first 10 genes) to ensure function will run properly 
#' priort to running all 18,333 genes. Default: FALSE.
#'
#' @return A data frame containing Pearson correlation coefficients. A copy is also saved to the 
#' directory defined in `output_dir`.
#' 
#' @details Description of output data frame
#' * `GeneNameID_A` - Hugo symbol with NCBI gene ID of query gene.
#' * `GeneNameID_B` - Hugo symbol with NCBI gene ID of all genes targeted in the DepMap KO screen.
#' * `estimate` - Correlation coefficient output from `?cor.test`.
#' * `statistic` - Pearson's correlation statistic. Output from `?cor.test`.
#' * `p.value` - P-value from Pearson's correlation statistic. Output from `?cor.test`.
#' * `parameter` - Degrees of freedom. Output from `?cor.test`.
#' * `conf.low` - Confidence interval low end. Output from `?cor.test`.
#' * `conf.high` - Confidence interval high end. Output from `?cor.test`.
#' * `method` - Type of cor.test used. Output from `?cor.test`.
#' * `alternative` - Alternative hypothesis selected. Output from `?cor.test`.
#' * `Rank` - Rank by correlation coefficient. 
#' * `Padj_BH` - Benjamini-Hochberg adjusted p-value.
#' @md
#' 
#' @examples 
#' \dontrun{
#' 
#' Screen_results <- GINI_screen(
#' Input_gene = "ARID1A", 
#' output_dir = "~/Desktop/GINI_test_dir/",
#' data_dir = "/path/to/DepMap_data/",
#' test = TRUE) # turn on for shorter test runs
#' 
#' }
#' 
#' @rdname coessential_map
#' @export 
#' @importFrom parallel detectCores
#' @importFrom doMC registerDoMC
#' @importFrom dplyr mutate filter pull rename arrange case_when
#' @importFrom foreach `%dopar%` foreach
#' @importFrom rcompanion cliffDelta
#' @importFrom diptest dip.test
#' @importFrom broom tidy
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect everything
#' @importFrom readr write_csv
#' @importFrom stats cor.test

coessential_map <- function(Input_gene = NULL, Input_disease = NULL, Input_cell_lines = NULL, core_num = NULL, output_dir = NULL, data_dir = NULL, test = FALSE){
  
  # Check that essential inputs are given:
  if(is.null(Input_gene)){
    stop("No control IDs detected")
  }
  if(is.null(core_num)){
    cores_detected <-  parallel::detectCores()
    print("No cores specified")
    print(paste0("Detected: ", cores_detected," cores"))
    print(paste0("Using: ", cores_detected/2," cores"))
    doMC::registerDoMC(cores_detected/2)
  }
  if(is.null(output_dir)){
    output_dir <- paste0(getwd(),"/GINIR_",Sys.Date())
    print(paste0("No output directory specified. Creating: ", output_dir))
    dir.create(output_dir)
  }
  if(!dir.exists(output_dir)){
    stop("Output directory does not exist. Please provide full path to directory.")
  }
  if(is.null(data_dir)){
    stop(paste0("No directory to data was specified. Please provide path to DepMap data."))
  }
  if(!dir.exists(data_dir)){
    stop(paste0("DepMap data directory does not exists. Please check again and provide the full path to the DepMap data directory."))
  }
  
  # Set cores:
  if(!is.null(core_num)){
    doMC::registerDoMC(core_num)
  }
  
  # Load necessary data
  gene_effect <- sample_annot <- NULL # see: https://support.bioconductor.org/p/24756/
  load(paste0(data_dir, "/gene_effect.rda"), envir = environment())
  load(paste0(data_dir, "/sample_annot.rda"), envir = environment())
  
  if(!is.null(Input_disease)){
    selected_cell_lines <- sample_annot %>%
      dplyr::filter(
        .data$DepMap_ID %in% gene_effect$DepMap_ID,
        .data$disease %in% Input_disease) %>%
      dplyr::pull(.data$DepMap_ID)
  } else if(!is.null(Input_cell_lines)){
    selected_cell_lines <- sample_annot %>%
      dplyr::filter(
        .data$DepMap_ID %in% gene_effect$DepMap_ID,
        .data$DepMap_ID %in% Input_cell_lines) %>%
      dplyr::pull(.data$DepMap_ID)
  } else {
    selected_cell_lines <- sample_annot %>%
      dplyr::filter(.data$DepMap_ID %in% gene_effect$DepMap_ID) %>%
      dplyr::pull(.data$DepMap_ID)
  }
  
  # account for differences between DepMap versions
  if(ncol(gene_effect) == 3){
    AllGenes <- unique(gene_effect$GeneNameID)
    gene_effect_long <- gene_effect %>%
      dplyr::filter(.data$DepMap_ID %in% selected_cell_lines)
    
  } else if(ncol(gene_effect) > 3) {
    AllGenes <- colnames(gene_effect)[-1] # removes DepMap_ID column  
    gene_effect_long <- gene_effect %>%
      tidyr::pivot_longer(
        cols = matches("\\d"),
        names_to = "GeneNameID",
        values_to = "Effect_score") %>%
      dplyr::filter(.data$DepMap_ID %in% selected_cell_lines)
    
  }
  
  Gene_A_GeneNameID <-  get_GeneNameID(Input_gene, data_dir = data_dir)
  Gene_A_effect <- gene_effect_long %>% dplyr::filter(.data$GeneNameID == Gene_A_GeneNameID)
  
  # Need to define function. A fix for a strange bug:
  `%dopar%` <- foreach::`%dopar%`
  
  # Begin loop
  print("This may take a few mins... Consider running with a higher core numbers to speed up the analysis.")
  if(test == TRUE){
    run <- 10
  } else {
    run <- length(unique(AllGenes))
  }
  res <- each <- NULL
  res <- foreach::foreach(each = 1:run, .combine = bind_rows) %dopar% {
    if(each == 1){
      print(paste0("Processing ", each, " of ", length(AllGenes),"\n"))
    } else if(each == length(AllGenes)){
      print(paste0("Processing ", each, " of ", length(AllGenes),"\n"))
    } else if(each%%1000 == 0){
      print(paste0("Processing ", each, " of ", length(AllGenes),"\n"))
    }
    
    Gene_B_effect <- gene_effect_long %>% dplyr::filter(.data$GeneNameID == AllGenes[each])
    
    res_pearson <- cor.test(Gene_A_effect$Effect_score, Gene_B_effect$Effect_score,
                            alternative = "two.sided",
                            method = "pearson",
                            na.action = "na.omit") %>% broom::tidy() %>%
      dplyr::mutate(GeneNameID_A = Gene_A_GeneNameID,
                    GeneNameID_B = AllGenes[each]) %>%
      dplyr::select(.data$GeneNameID_A, .data$GeneNameID_B, tidyr::everything())
    
    res_pearson
  }
  
  # save and return output
  output <- res %>%
    dplyr::arrange(-.data$estimate, .data$p.value) %>%
    dplyr::mutate(
      Rank = order(-.data$estimate, decreasing = F),
      Padj_BH = p.adjust(.data$p.value, method = "BH", n = (length(.data$p.value)))) %>%
    readr::write_csv(file = paste0(output_dir,"/GINI_coessentiality_network_results.csv"))
  
  print(paste0("Coessentiality mapping finished. Outputs were written to: ", output_dir,"/GINI_coessentiality_network_results.csv"))
  return(output)
}
