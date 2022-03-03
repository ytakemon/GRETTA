#' Get unique DepMap-compatible gene IDs
#' 
#' @description 
#' `get_GeneNameID()` and `get_DepMapID()` provide tools for converting gene symbols/ids and 
#' common cell line aliases to DepMap's unique gene name and cell line identifiers, respectively.
#' 
#' @param GeneName string containing a vector of either Hugo gene symbol or numeric NCBI ID
#' @param data_dir string Path to GINIR_data
#' @return string
#' 
#' @import rlang
#' @import dplyr
#' @import utils
#' 
#' @export
#' @examples
#' \dontrun{
#' get_GeneNameID("A1CF", data_dir = "/path/")
#' }
get_GeneNameID <- function(GeneName, data_dir){
  # Load necessary data
  dep_annot <- CCLE_exp_annot <- NULL # see: https://support.bioconductor.org/p/24756/
  load(paste0(data_dir, "/dep_annot.rda"), envir = environment())
  load(paste0(data_dir, "/CCLE_exp_annot.rda"), envir = environment())
  
  # For Hugo symbols/NCBI IDs - from dependency prob.
  if(any(dep_annot$GeneNames %in% GeneName | dep_annot$GeneID %in% GeneName)){
    res <- dep_annot %>%
      dplyr::filter(.data$GeneNames %in% GeneName | .data$GeneID %in% GeneName) %>%
      dplyr::pull(.data$GeneNameID)
    
  # For Hugo IDs - from CCLE RNA-seq expr.
  } else if(any(CCLE_exp_annot$GeneNames %in% GeneName | CCLE_exp_annot$GeneID %in% GeneName)){
    res <- CCLE_exp_annot %>%
      dplyr::filter(.data$GeneNames %in% GeneName | .data$GeneID %in% GeneName) %>%
      dplyr::pull(.data$GeneNameID)
    
  # If no matches are found
  } else {
    stop("Cannot find gene name. Please check the spelling. Gene names should be a valid Hugo Symbol and in all caps!")
  }
  return(res)
}

#' @describeIn get_GeneNameID Get unique DepMap-compatible sample IDs
#' @param sample_name string containing a vector of unique sample_id used in proteomics data or common cell line names
#' @param data_dir string Path to GINIR_data
#' @return string
#' 
#' @import rlang
#' @import dplyr
#' @import utils
#' 
#' @export
#' @examples
#' \dontrun{
#' get_DepMapID("JURKAT", data_dir = "/path/")
#' }
get_DepMapID <- function(sample_name, data_dir){
  # Load necessary data
  protein_annot <- sample_annot <- NULL # see: https://support.bioconductor.org/p/24756/
  load(paste0(data_dir, "/protein_annot.rda"), envir = environment())
  load(paste0(data_dir, "/sample_annot.rda"), envir = environment())
  
  # For sample names from proteomics data 
  if(any(protein_annot$GygiNames %in% sample_name)){
    res <- protein_annot %>% 
      dplyr::filter(.data$GygiNames %in% sample_name) %>% 
      dplyr::pull(.data$DepMap_ID)
  
  # For sample names that are commonly used in the wild
  } else if(any(sample_annot$stripped_cell_line_name %in% sample_name)){
    res <- sample_annot %>% 
      dplyr::filter(.data$stripped_cell_line_name %in% sample_name) %>% 
      dplyr::pull(.data$DepMap_ID)
    
  # If no matches are found
  } else {
    stop("Cannot find sample. Please check the spelling. Common sample names should be all caps without spaces or special symbols!")
  }
  return(res)
}
