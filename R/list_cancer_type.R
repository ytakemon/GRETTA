#' List cancer types that are available
#' 
#' @description 
#' `list_cancer_types()` and `list_cancer_subtypes()` provide tools for identifying cancer (sub)types that are available in DepMap.
#' 
#' @return string A vector containing unique cancer types available

#' @param data_dir string Path to GINIR_data
#' @import rlang
#' @import dplyr
#' @import utils
#' 
#' @export
#' @examples
#' gretta_data_dir <- './GRETTA_example/'
#' gretta_output_dir <- './GRETTA_example_output/'
#' 
#' if(!dir.exists(gretta_data_dir)){
#'   download_example_data(".")
#' }
#' 
#' list_cancer_types(data_dir = gretta_data_dir)
#' 
list_cancer_types <- function(data_dir) {
  # Load necessary data
  sample_annot <- NULL  # see: https://support.bioconductor.org/p/24756/
  load(paste0(data_dir, "/sample_annot.rda"), envir = environment())

  # Main
  sample_annot %>%
    dplyr::pull(.data$disease) %>%
    unique
}

#' @describeIn list_cancer_types List cancer subtypes that are available
#' 
#' @param input_disease string A vector of unique with one or more cancer types listed in `list_cancer_types()`
#' @param data_dir string Path to GINIR_data
#' @importFrom rlang .data
#' @export
#' @examples
#' gretta_data_dir <- './GRETTA_example/'
#' gretta_output_dir <- './GRETTA_example_output/'
#' 
#' list_cancer_subtypes('Lung Cancer', data_dir = gretta_data_dir)
#' 
list_cancer_subtypes <- function(input_disease, data_dir) {
  if (is.null(data_dir)) {
    stop("No directory to data was specified. ",
         "Please provide path to DepMap data.")
  }
  if (!dir.exists(data_dir)) {
    stop("DepMap data directory does not exists. ",
         "Please check again and provide the full path to the DepMap data directory.")
  }
  # Load necessary data
  sample_annot <- NULL  # see: https://support.bioconductor.org/p/24756/
  load(paste0(data_dir, "/sample_annot.rda"), envir = environment())
  
  # Main
  sample_annot %>%
    dplyr::filter(.data$disease %in% input_disease) %>%
    dplyr::pull(.data$disease_subtype) %>%
    unique
  
  # Message:
  message("See sample_annot.rda in your data directory for a complist list of cell line annotation.")
}
