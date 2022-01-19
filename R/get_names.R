#' Get unique DepMap-compatible names
#' 
#' @description 
#' `get_GeneNameID()` and `get_DepMapID()` provide tools for converting gene symbols/ids and 
#' common cell line aliases to DepMap's unique gene name and cell line identifiers, respectively.
#' 
#' @param GeneName string containing a vector of either Hugo gene symbol or numeric NCBI ID
#' @return string
#' @export
#' @examples
#' get_GeneNameID("A1CF")
get_GeneNameID <- function(GeneName){
    # For Hugo symbols/NCBI IDs - from dependency prob.
  if(any(dep_annot$GeneNames %in% GeneName | dep_annot$GeneID %in% GeneName)){
    res <- dep_annot %>%
      dplyr::filter(GeneNames %in% GeneName | GeneID %in% GeneName) %>%
      dplyr::pull(GeneNameID)
    
    # For Hugo IDs - from CCLE RNA-seq expr.
  } else if(any(CCLE_exp_annot$GeneNames %in% GeneName | CCLE_exp_annot$GeneID %in% GeneName)){
    res <- CCLE_exp_annot %>%
      dplyr::filter(GeneNames %in% GeneName | GeneID %in% GeneName) %>%
      dplyr::pull(GeneNameID)
    
    # If no matches are found
  } else {
    stop("Cannot find gene name. Please check the spelling. Gene names should be a valid Hugo Symbol and in all caps!")
  }
  return(res)
}

#' @export
#' @rdname get_GeneNameID
#' @param sample_name string containing a vector of unique sample_id used in proteomics data or common cell line names
#' @return string
#' @examples
#' get_DepMapID("JURKAT")
get_DepMapID <- function(sample_name){
  if(any(protein_annot$GygiNames %in% sample_name)){
    res <- protein_annot %>% 
      dplyr::filter(GygiNames %in% x) %>% 
      dplyr::pull(DepMap_ID)
  } else if(any(sample_annot$stripped_cell_line_name %in% sample_name)){
    res <- sample_annot %>% 
      dplyr::filter(stripped_cell_line_name %in% sample_name) %>% 
      dplyr::pull(DepMap_ID)
    
    # If no matches are found
  } else {
    stop("Cannot find sample. Please check the spelling. Common sample names should be all caps without spaces or special symbols!")
  }
  return(res)
}
