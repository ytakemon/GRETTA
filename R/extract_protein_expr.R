#' @title Extract protein expression data for given cell lines
#' 
#' @description FUNCTION_DESCRIPTION
#' 
#' @param input_samples string A vector of DepMap_ID(s) must be provided, Default: NULL
#' @param input_genes string Optional Hugo Symbol(s) encoding proteins of interest, Default: NULL
#' @param data_dir string Path to GINIR_data
#' @return Data frame containing protein expression for samples provided in the input. 
#' If no genes were specified, the function will return a data frame of proteins profiled in DepMap
#' 
#' @details See also `extract_protein_rna` to extract proteomics profile data
#' 
#' @examples 
#' \dontrun{
#' extract_protein_expr(
#' input_samples = c('ACH-000004', 'ACH-000146'), 
#' input_genes = c('ATM','TOP1'),
#' data_dir = '/path/to/DepMap_data/')
#' }
#' @rdname extract_protein_expr
#' @export 
#' @importFrom dplyr select contains left_join filter
#' @importFrom tidyr pivot_longer

extract_protein_expr <- function(input_samples = NULL, input_genes = NULL, data_dir = NULL) {
  
  # Print and check to see input was provided
  if (is.null(input_samples)) {
    stop("No samples given. Please input sample DepMap_ID")
  }
  if (is.null(data_dir)) {
    stop(
      paste0("No directory to data was specified. Please provide path to DepMap data.")
    )
  }
  if (!dir.exists(data_dir)) {
    stop(
      paste0(
        "DepMap data directory does not exists. Please check again and provide the full path to the DepMap data directory."
      )
    )
  }
  
  # Load necessary data
  protein_annot <- protein_nodup <- sample_annot <- NULL  # see: https://support.bioconductor.org/p/24756/
  load(
    paste0(data_dir, "/sample_annot.rda"),
    envir = environment()
  )
  load(
    paste0(data_dir, "/protein_nodup.rda"),
    envir = environment()
  )
  load(
    paste0(data_dir, "/protein_annot.rda"),
    envir = environment()
  )
  
  # Check if inputs are recognized
  if (!all(input_samples %in% sample_annot$DepMap_ID)) {
    stop(
      paste0(
        input_samples[!input_samples %in% sample_annot$DepMap_ID], ", not recognized as a valid sample"
      )
    )
  }
  if (!all(input_genes %in% protein_nodup$Gene_Symbol)) {
    stop(
      paste0(
        input_genes[!input_genes %in% protein_nodup$Gene_Symbol], ", not recognized or protein expression is not available. Please check spelling or remove gene name from input"
      )
    )
  }
  
  # If no input gene is given, give full expr table
  if (is.null(input_genes)) {
    res <- protein_nodup %>%
      dplyr::select(
        .data$Gene_Symbol, .data$Description, .data$Uniprot, .data$Uniprot_Acc,
        dplyr::contains("_TenPx")
      ) %>%
      tidyr::pivot_longer(
        -c(.data$Gene_Symbol, .data$Description, .data$Uniprot, .data$Uniprot_Acc),
        names_to = "Gygi_ID", values_to = "protein_expr"
      ) %>%
      dplyr::left_join(protein_annot, by = c(Gygi_ID = "GygiNames")) %>%
      dplyr::select(.data$DepMap_ID, dplyr::everything())
    
    return(res)
  }
  
  # Otherwise, provide only expr of genes of interst
  res <- protein_nodup %>%
    dplyr::filter(.data$Gene_Symbol %in% input_genes) %>%
    dplyr::select(
      .data$Gene_Symbol, .data$Description, .data$Uniprot, .data$Uniprot_Acc,
      dplyr::contains("_TenPx")
    ) %>%
    tidyr::pivot_longer(
      -c(.data$Gene_Symbol, .data$Description, .data$Uniprot, .data$Uniprot_Acc),
      names_to = "Gygi_ID", values_to = "protein_expr"
    ) %>%
    dplyr::left_join(protein_annot, by = c(Gygi_ID = "GygiNames")) %>%
    dplyr::filter(.data$DepMap_ID %in% input_samples,) %>%
    dplyr::select(.data$DepMap_ID, dplyr::everything())
  
  # Notify if some samples do not have expression data
  if (!all(input_samples %in% res$DepMap_ID)) {
    print(
      paste0(
        "Following sample did not contain profile data: ", paste0(input_samples[!input_samples %in% res$DepMap_ID], collapse = ", ")
      )
    )
    return(res)
  } else {
    return(res)
  }
}

