#' @title Extract RNA expression data for given cell lines
#' 
#' @description FUNCTION_DESCRIPTION
#' 
#' @param input_samples string A vector of DepMap_ID(s) must be provided, Default: NULL
#' @param input_genes string Optional Hugo Symbol(s), Default: NULL
#' @param data_dir string Path to GINIR_data
#' @return Data frame containing RNA expression (TPM) for sample provided in the input. 
#' If no genes were specified, the function will return a data frame of all genes profiled in DepMap
#' 
#' @details See also `extract_prot` to extract proteomics profile data
#' 
#' @examples 
#' gretta_data_dir <- './GRETTA_example/'
#' gretta_output_dir <- './GRETTA_example_output/'
#' 
#' if(!dir.exists(gretta_data_dir)){
#'   download_example_data(".")
#' }
#' 
#' extract_rna(
#' input_samples = c('ACH-001642','ACH-000688'), 
#' input_genes = c('ARID1A'),
#' data_dir = gretta_data_dir)
#' 
#' @rdname extract_rna
#' @export 
#' @importFrom dplyr filter select

extract_rna <- function(input_samples = NULL, input_genes = NULL,
                        data_dir = NULL) {
  
  # Print and check to see input was provided
  if (is.null(input_samples)) {
    stop("No samples given. Please input sample DepMap_ID")
  }
  if (is.null(data_dir)) {
    stop("No directory to data was specified. Please provide path to DepMap data.")
  }
  if (!dir.exists(data_dir)) {
    stop("DepMap data directory does not exists.",
         "Please check again and provide the full path to the DepMap data directory.")
  }
  
  # Load necessary data
  CCLE_exp <- CCLE_exp_annot <- sample_annot <- NULL  # see: https://support.bioconductor.org/p/24756/
  load(paste0(data_dir, "/CCLE_exp.rda"), envir = environment())
  load(paste0(data_dir, "/CCLE_exp_annot.rda"), envir = environment())
  load(paste0(data_dir, "/sample_annot.rda"), envir = environment())
  
  # Check if inputs are recognized
  if (!all(input_samples %in% sample_annot$DepMap_ID)) {
    stop(input_samples[!input_samples %in% sample_annot$DepMap_ID],
         ", not recognized as a valid sample")
  }
  if (!all(input_genes %in% CCLE_exp_annot$GeneNames)) {
    stop(input_genes[!input_genes %in% CCLE_exp_annot$GeneNames],
         ", not recognized. Please check spelling or remove gene name from input")
  }
  
  # If no input gene is given, give full expr
  # table
  if (is.null(input_genes)) {
    res <- CCLE_exp %>%
      dplyr::filter(.data$DepMap_ID %in% input_samples)
    
    return(res)
  }
  
  # Otherwise, provide only expr of genes of
  # interst
  gene_colname <- CCLE_exp_annot %>% 
    dplyr::filter(.data$GeneNames %in% input_genes) %>% 
    dplyr::pull(.data$GeneNameID)
  
  res <- CCLE_exp %>%
    dplyr::select(.data$DepMap_ID, all_of(gene_colname)) %>%
    dplyr::filter(.data$DepMap_ID %in% input_samples) %>%
    ungroup()
  
  # Notify if some samples do not have
  # expression data
  if (!all(input_samples %in% res$DepMap_ID)) {
    GRETTA_says <- paste0("Following sample did not contain RNA data: ",
                          paste0(input_samples[!input_samples %in%
                                                 res$DepMap_ID], collapse = ", "))
    message(GRETTA_says)
    return(res)
    
  } else {
    return(res)
  }
}