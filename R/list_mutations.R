#' @title Search mutations detected in DepMap cancer cell lines 
#' 
#' @description Quick method to look up mutations of interest prior to using `select_cell_lines()`
#' 
#' @param gene string, Hugo Symbol, Default: NULL
#' @param chr string, Default: NULL
#' @param start_bp integer, Default: NULL
#' @param end_bp integer, Default: NULL
#' @param is_hotspot logical, TCGA or COSMIC hotspot, Default: NULL
#' @param is_damaging logical, Default: NULL 
#' @param variant_classification string, Select from the following: 3'UTR, 5'Flank, 5'UTR, 
#' De_novo_Start_OutOfFrame, Frame_Shift_Del, Frame_Shift_Ins, IGR, In_Frame_Del, In_Frame_Ins, 
#' Intron, Missense_Mutation, Nonsense_Mutation, Nonstop_Mutation, Silent, Splice_Site, 
#' Start_Codon_Del, Start_Codon_Ins, Start_Codon_SNP, Stop_Codon_Del, Stop_Codon_Ins, Default: NULL
#' @param data_dir string Path to GINIR_data
#' 
#' @return A data frame containing mutations matching criteria of input arguments
#' 
#' @examples 
#' gretta_data_dir <- './GRETTA_example/'
#' gretta_output_dir <- './GRETTA_example_output/'
#' 
#' if(!dir.exists(gretta_data_dir)){
#'   download_example_data(".")
#' }
#' 
#' list_mutations(
#' gene = 'ARID1A', 
#' is_damaging = TRUE,
#' data_dir = gretta_data_dir)
#' 
#' @rdname list_mutations
#' @export 
#' @importFrom dplyr filter select arrange distinct
#' 
list_mutations <- function(gene = NULL, chr = NULL,
                           start_bp = NULL, end_bp = NULL, is_hotspot = NULL,
                           is_damaging = NULL, variant_classification = NULL,
                           data_dir = NULL) {
  
  # Print and check to see input
  if (is.null(data_dir)) {
    stop("No directory to data was specified. Please provide path to DepMap data.")
  }
  if (!dir.exists(data_dir)) {
    stop("DepMap data directory does not exists. ",
         "Please check again and provide the full path to the DepMap data directory.")
  }
  if (is.null(c(gene, chr))) {
    stop("No input gene or region given. Please prvide a Hugo gene symbol.")
  }
  
  # Load necessary data
  mut_calls <- NULL  # see: https://support.bioconductor.org/p/24756/
  load(paste0(data_dir, "/mut_calls.rda"), envir = environment())
  
  # If gene is provided look for mutations:
  if (!is.null(gene))
  {
    # Check if input gene mutations exist
    if (!any(mut_calls$Hugo_Symbol %in% gene)) {
      stop("No mutations were found for: ",
           gene, ". Please check spelling and for valid Hugo Symbols")
    }
    
    # Get ALL Mutations, and apply
    # additional filters if they exist
    target_mut <- mut_calls %>%
      dplyr::filter(.data$Hugo_Symbol %in%
                      gene) %>%
      dplyr::select(.data$DepMap_ID, .data$Hugo_Symbol:.data$Annotation_Transcript,
                    .data$cDNA_Change:.data$COSMIChsCnt,
                    .data$Variant_annotation) %>%
      dplyr::arrange(.data$Start_position) %>%
      dplyr::distinct()
    
    if (!is.null(chr)) {
      target_mut <- target_mut %>%
        dplyr::filter(.data$Chromosome %in%
                        as.character(chr))
    }
    if (!is.null(start_bp)) {
      target_mut <- target_mut %>%
        dplyr::filter(.data$Start_position >=
                        start_bp)
    }
    if (!is.null(end_bp)) {
      target_mut <- target_mut %>%
        dplyr::filter(.data$End_position <=
                        end_bp)
    }
    if (!is.null(is_hotspot)) {
      target_mut <- target_mut %>%
        dplyr::filter(.data$isTCGAhotspot ==
                        is_hotspot | .data$isCOSMIChotspot ==
                        is_hotspot)
    }
    if (!is.null(is_damaging)) {
      target_mut <- target_mut %>%
        dplyr::filter(.data$isDeleterious ==
                        is_damaging)
    }
    if (!is.null(variant_classification)) {
      target_mut <- target_mut %>%
        dplyr::filter(.data$Variant_Classification %in%
                        variant_classification)
    }
  }  # End of: If gene is provided look for mutations:
  
  # If no gene is provided, but chr is provided
  if (is.null(gene) & !is.null(chr))
  {
    
    # Check if input chr exists
    if (!any(mut_calls$Chromosome %in% as.character(chr))) {
      stop("No mutations were found for: ",
           gene, ". Please check spelling and for valid Hugo Symbols")
    }
    
    # Get ALL Mutations, and apply
    # additional filters if they exist
    target_mut <- mut_calls %>%
      dplyr::filter(.data$Chromosome %in%
                      as.character(chr)) %>%
      dplyr::select(.data$DepMap_ID, .data$Hugo_Symbol:.data$Annotation_Transcript,
                    .data$cDNA_Change:.data$COSMIChsCnt,
                    .data$Variant_annotation) %>%
      dplyr::arrange(.data$Start_position) %>%
      dplyr::distinct()
    
    if (!is.null(start_bp)) {
      target_mut <- target_mut %>%
        dplyr::filter(.data$Start_position >=
                        start_bp)
    }
    if (!is.null(end_bp)) {
      target_mut <- target_mut %>%
        dplyr::filter(.data$End_position <=
                        end_bp)
    }
    if (!is.null(is_hotspot)) {
      target_mut <- target_mut %>%
        dplyr::filter(.data$isTCGAhotspot ==
                        is_hotspot | .data$isCOSMIChotspot ==
                        is_hotspot)
    }
    if (!is.null(is_damaging)) {
      target_mut <- target_mut %>%
        dplyr::filter(.data$isDeleterious ==
                        is_damaging)
    }
    if (!is.null(variant_classification)) {
      target_mut <- target_mut %>%
        dplyr::filter(.data$Variant_Classification %in%
                        variant_classification)
    }
  }  # End of: If no gene is provided, but chr is provided
  
  return(target_mut)
}
