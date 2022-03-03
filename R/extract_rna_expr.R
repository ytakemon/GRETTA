#' @title Extract RNA expression data for given cell lines
#' 
#' @description FUNCTION_DESCRIPTION
#' 
#' @param Input_samples string A vector of DepMap_ID(s) must be provided, Default: NULL
#' @param Input_genes string Optional Hugo Symbol(s), Default: NULL
#' @param data_dir string Path to GINIR_data
#' @return Data frame containing RNA expression (TPM) for sample provided in the input. 
#' If no genes were specified, the function will return a data frame of all genes profiled in DepMap
#' 
#' @details See also `extract_protein_expr` to extract proteomics profile data
#' 
#' @examples 
#' \dontrun{
#' extract_rna_expr(Input_samples = c("ACH-001642","ACH-000688"), Input_genes = c("TP53","ARID1A"))
#' }
#' @rdname extract_rna_expr
#' @export 
#' @importFrom dplyr filter select

extract_rna_expr <- function(Input_samples = NULL, Input_genes = NULL, data_dir = NULL){

   # Print and check to see input was provided
  if(is.null(Input_samples)){
    stop("No samples given. Please input sample DepMap_ID")
  } 
  
  # Load necessary data
  CCLE_exp <- CCLE_exp_annot <- sample_annot <- NULL # see: https://support.bioconductor.org/p/24756/
  load(paste0(data_dir,"/CCLE_exp.rda"), envir = environment())
  load(paste0(data_dir,"/CCLE_exp_annot.rda"), envir = environment())
  load(paste0(data_dir,"/sample_annot.rda"), envir = environment())
  
  # Check if inputs are recognized 
  if(!all(Input_samples %in% sample_annot$DepMap_ID)){
    stop(paste0(Input_samples[!Input_samples %in% sample_annot$DepMap_ID],", not recognized as a valid sample"))
  }
  if(!all(Input_genes %in% CCLE_exp_annot$GeneNames)){
    stop(paste0(Input_genes[!Input_genes %in% CCLE_exp_annot$GeneNames],", not recognized. Please check spelling or remove gene name from input"))
  }
  
  # If no input gene is given, give full expr table
  if(is.null(Input_genes)){
    res <- CCLE_exp %>%
      dplyr::filter(.data$DepMap_ID %in% Input_samples) 
    
    return(res)
  }
  
  # Otherwise, provide only expr of genes of interst 
  res <- CCLE_exp %>%
    dplyr::select(.data$DepMap_ID, get_GeneNameID(Input_genes, data_dir = data_dir)) %>%
    dplyr::filter(.data$DepMap_ID %in% Input_samples)
  
  # Notify if some samples do not have expression data 
  if(!all(Input_samples %in% res$DepMap_ID)){
    print(paste0("Following sample did not contain profile data: ", 
                 paste0(Input_samples[!Input_samples %in% res$DepMap_ID], collapse = ", ")))
    return(res)
  } else {
    return(res)
  }
}