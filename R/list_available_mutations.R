#' @title Search mutations detected in DepMap cancer cell lines 
#' 
#' @description Quick method to look up mutations of interest prior to using `select_cell_lines()`
#' 
#' @param Gene string, Hugo Symbol, Default: NULL
#' @param Chr string, Default: NULL
#' @param Start_bp integer, Default: NULL
#' @param End_bp integer, Default: NULL
#' @param Is_hotspot logical, TCGA or COSMIC hotspot, Default: NULL
#' @param Is_damaging logical, Default: NULL 
#' @param Variant_classification string, Select from the following: 3'UTR, 5'Flank, 5'UTR, 
#' De_novo_Start_OutOfFrame, Frame_Shift_Del, Frame_Shift_Ins, IGR, In_Frame_Del, In_Frame_Ins, 
#' Intron, Missense_Mutation, Nonsense_Mutation, Nonstop_Mutation, Silent, Splice_Site, 
#' Start_Codon_Del, Start_Codon_Ins, Start_Codon_SNP, Stop_Codon_Del, Stop_Codon_Ins, Default: NULL
#' 
#' @return A data frame containing mutations matching criteria of input arguments
#' 
#' @examples 
#' \dontrun{
#' list_available_mutations(Gene = "TP53")
#' list_available_mutations(Chr = 12)
#' }
#' 
#' @rdname list_available_mutations
#' @export 
#' @importFrom dplyr filter select arrange distinct
list_available_mutations <- function(Gene = NULL, 
                             Chr = NULL, Start_bp = NULL, End_bp = NULL,
                             Is_hotspot = NULL, 
                             Is_damaging = NULL,
                             Variant_classification = NULL){
  
  # Print and check to see input
  if(is.null(c(Gene, Chr))){
    stop("No input gene or region given. Please prvide a Hugo gene symbol.")
    
  }
  
  # Load necessary data
  mut_calls <- NULL # see: https://support.bioconductor.org/p/24756/
  load("data/mut_calls.rda", envir = environment())
  
  # If Gene is provided look for mutations:
  if(!is.null(Gene)){
    # Check if input gene mutations exist
    if(!any(mut_calls$Hugo_Symbol %in% Gene)){
      stop(paste0("No mutations were found for: ", Gene,". Please check spelling and for valid Hugo Symbols"))
    }
    
    # Get ALL Mutations, and apply additional filters if they exist
    target_mut <- mut_calls %>%
      dplyr::filter(.data$Hugo_Symbol %in% Gene) %>%
      dplyr::select(.data$Hugo_Symbol:.data$Annotation_Transcript, 
                    .data$cDNA_Change:.data$COSMIChsCnt, 
                    .data$Variant_annotation) %>%
      dplyr::arrange(.data$Start_position) %>%
      dplyr::distinct()
    
    if(!is.null(Chr)){
      target_mut <- target_mut %>% dplyr::filter(.data$Chromosome %in% as.character(Chr))
    }
    if(!is.null(Start_bp)){
      target_mut <- target_mut %>% dplyr::filter(.data$Start_position >= Start_bp)
    }
    if(!is.null(End_bp)){
      target_mut <- target_mut %>% dplyr::filter(.data$End_position <= End_bp)
    }
    if(!is.null(Is_hotspot)){
      target_mut <- target_mut %>% dplyr::filter(.data$isTCGAhotspot == Is_hotspot |
                                                   .data$isCOSMIChotspot == Is_hotspot)
    }
    if(!is.null(Is_damaging)){
      target_mut <- target_mut %>% dplyr::filter(.data$isDeleterious == Is_damaging)
    }
    if(!is.null(Variant_classification)){
      target_mut <- target_mut %>% dplyr::filter(.data$Variant_Classification %in% Variant_classification)
    }
  } # End of: If Gene is provided look for mutations:
  
  # If no gene is provided, but Chr is provided
  if(is.null(Gene) & !is.null(Chr)){
    
    # Check if input chr exists
    if(!any(mut_calls$Chromosome %in% as.character(Chr))){
      stop(paste0("No mutations were found for: ", Gene,". Please check spelling and for valid Hugo Symbols"))
    }
    
    # Get ALL Mutations, and apply additional filters if they exist
    target_mut <- mut_calls %>%
      dplyr::filter(.data$Chromosome %in% as.character(Chr)) %>%
      dplyr::select(.data$Hugo_Symbol:.data$Annotation_Transcript, 
                    .data$cDNA_Change:.data$COSMIChsCnt, 
                    .data$Variant_annotation) %>%
      dplyr::arrange(.data$Start_position) %>%
      dplyr::distinct()
    
    if(!is.null(Start_bp)){
      target_mut <- target_mut %>% dplyr::filter(.data$Start_position >= Start_bp)
    }
    if(!is.null(End_bp)){
      target_mut <- target_mut %>% dplyr::filter(.data$End_position <= End_bp)
    }
    if(!is.null(Is_hotspot)){
      target_mut <- target_mut %>% dplyr::filter(.data$isTCGAhotspot == Is_hotspot |
                                                   .data$isCOSMIChotspot == Is_hotspot)
    }
    if(!is.null(Is_damaging)){
      target_mut <- target_mut %>% dplyr::filter(.data$isDeleterious == Is_damaging)
    }
    if(!is.null(Variant_classification)){
      target_mut <- target_mut %>% dplyr::filter(.data$Variant_Classification %in% Variant_classification)
    }
  } # End of: If no gene is provided, but Chr is provided
  
  return(target_mut)
}
