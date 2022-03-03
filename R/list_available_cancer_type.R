#' List cancer types that are available
#' 
#' @description 
#' `list_available_cancer_types()` and `list_available_cancer_subtypes()` provide tools for identifying cancer (sub)types that are available in DepMap.
#' 
#' @return string A vector containing unique cancer types available
#' 
#' @import rlang
#' @import dplyr
#' @import utils
#' 
#' @export
#' @examples
#' list_available_cancer_types()
list_available_cancer_types <- function(){
  # Load necessary data
  sample_annot <- NULL # see: https://support.bioconductor.org/p/24756/
  load(paste0(system.file(package = "GINIR"), "/data/sample_annot.rda"), envir = environment())
  
  # Main
  sample_annot %>% 
    dplyr::pull(.data$disease) %>% unique
}

#' @describeIn list_available_cancer_types List cancer subtypes that are available
#' 
#' @param input_disease string A vector of unique with one or more cancer types listed in `list_available_cancer_types()`
#' @importFrom rlang .data
#' @export
#' @examples
#' list_available_cancer_subtypes("Lung Cancer")
list_available_cancer_subtypes <- function(input_disease){
  # Load necessary data
  sample_annot <- NULL # see: https://support.bioconductor.org/p/24756/
  load(paste0(system.file(package = "GINIR"), "/data/sample_annot.rda"), envir = environment())
  
  # Main
  sample_annot %>% 
    dplyr::filter(.data$disease %in% input_disease) %>%
    dplyr::pull(.data$disease_subtype) %>% unique
}
