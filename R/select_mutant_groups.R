#' Select mutant groups based on input gene of interest
#' 
#' @description
#' `select_mutant_groups()` assigns cancer cell lines to either `Control` groups or one of the following mutant groups: `HomDel`,
#' `T-HetDel`, `HetDel`, `Amplified`, or `Others` (see below for more details).
#' 
#' @param Input_gene string Input Hugo Symbol
#' @return data frame containing a summary of mutations found in cell lines and their control and mutant group assignments.
#' @import rlang
#' @import dplyr
#' @import utils
#' 
#' @export
#' @details Mutant groups in more detail:
#' * Control` cell lines do not harbor any single nucleotide variations (SNVs) or insertions and deletions (InDels) with a neutral copy number (CN).
#' * `HomDel` cell lines harbor one or more homozygous deleterious SNVs or have deep CN loss.
#' * `T-HetDel` cell lines harbor two or more heterozygous deleterious SNVs/InDels with neutral or CN loss.
#' * `HetDel` cell lines harbor one heterozygous deleterious SNV/InDel with neutral CN, or no SNV/InDel with CN loss. 
#' * `Amplified` cell lines harbor no SNVs/InDels with increased CN.
#' * `Others` cell lines harbor deleterious SNVs with increased CN.
#' 
#' @examples
#' A1CF_groups <- select_mutant_groups("A1CF")

select_mutant_groups <- function(Input_gene){
  # Check to see input is given
  if(missing(Input_gene)){
    stop("Input gene is missing")
  } else {
    # Print input
    print(paste0("Selecting mutant groups for: ", Input_gene))
  }
  
  # Load necessary data
  mut_calls <- copy_num_annot <- copy_num <- dep <- sample_annot <- NULL # see: https://support.bioconductor.org/p/24756/
  data(list = list("mut_calls", "copy_num_annot", "copy_num", "dep", "sample_annot"), envir = environment())
  
  # Check if input gene mutations exist
  if(!any(mut_calls$Hugo_Symbol %in% Input_gene)|!any(copy_num_annot$GeneNames %in% Input_gene)){
    stop(paste0("No mutations were found for: ", Input_gene,"\nPlease check spelling and for valid Hugo Symbols"))
  }
  # Convert to unique geneID
  Input_geneID <- get_GeneNameID(Input_gene)

  # Get copy number
  if(any(names(copy_num) %in% Input_geneID)){
    target_copy_num <- copy_num %>%
      dplyr::select(.data$DepMap_ID, dplyr::all_of(Input_geneID)) %>%
      dplyr::filter(.data$DepMap_ID %in% dep$DepMap_ID) %>%
      dplyr::arrange(.data$DepMap_ID) %>%
      dplyr::mutate( Status = dplyr::case_when(
        !!as.name(Input_geneID) <= 0.25 ~ "Deep_del",
        !!as.name(Input_geneID) > 0.25 & !!as.name(Input_geneID) < 0.75 ~ "Loss",
        !!as.name(Input_geneID) >= 0.75 & !!as.name(Input_geneID) < 1.25 ~ "Neutral",
        !!as.name(Input_geneID) >= 1.25 ~ "Amplified",
        TRUE ~ "Other"))
    
  } else {
    target_copy_num <- dep %>% 
      dplyr::select(.data$DepMap_ID) %>%
      dplyr::arrange(.data$DepMap_ID) %>%
      dplyr::mutate(!!as.name(Input_geneID) := NA,
             Status = "Unknown")
  }
  
  # Get ALL Mutations
  target_mut <- mut_calls %>%
    dplyr::filter((.data$DepMap_ID %in% dep$DepMap_ID) & (.data$Hugo_Symbol %in% Input_gene)) %>%
    dplyr::mutate(AC_combined = dplyr::coalesce(.data$CGA_WES_AC, .data$SangerRecalibWES_AC, .data$SangerWES_AC, .data$RNAseq_AC, .data$HC_AC, .data$RD_AC, .data$WGS_AC), #(Alt:REF)
           AC_ref_NULL = grepl(":0", .data$AC_combined)) %>%
    dplyr::mutate(AC_Variant = dplyr::case_when(
      .data$AC_ref_NULL == "TRUE" ~ "Hom_Mut",
      TRUE ~ "Het_Mut"),
      AC_Variant = paste0(.data$Variant_Classification," ", .data$AC_Variant)) %>% # are there any with 0 contribution from reference?
    dplyr::select(.data$Hugo_Symbol, .data$Chromosome:.data$Annotation_Transcript, .data$cDNA_Change:.data$COSMIChsCnt, .data$Variant_annotation:.data$AC_Variant) %>%
    dplyr::arrange(.data$Start_position) %>%
    dplyr::select(.data$DepMap_ID, everything())
  # Count number of mutations found per sample
  all_mutations_count_by_sample <- target_mut %>% dplyr::count(.data$DepMap_ID)
  
  # Get deleterious/damaging mutations
  mut_dels <- target_mut %>% dplyr::filter(.data$Variant_annotation == "damaging")
  # Count number of deleterious mutations found per sample
  del_mutations_count_by_sample <- mut_dels %>% dplyr::count(.data$DepMap_ID)
  # Find samples with homozygous deleterious mutations (HomDels)
  hom_del_muts_by_sample <- mut_dels %>% dplyr::select(.data$DepMap_ID, .data$AC_ref_NULL) %>%
    dplyr::distinct() %>%
    dplyr::filter(.data$AC_ref_NULL == TRUE)
  
  # Find samples with multiple heterozygous deleterious mutations. Trans-heterozygous mutants (T-HetDels)
  multi_mut_dels <- mut_dels %>%
    dplyr::filter(.data$AC_ref_NULL == FALSE) %>%
    dplyr::add_count(.data$DepMap_ID) %>% 
    dplyr::filter(.data$n > 1)
  
  # Summarize mutations for all samples 
  summary <- sample_annot %>% 
    dplyr::filter(.data$DepMap_ID %in% dep$DepMap_ID) %>%
    dplyr::select(.data$DepMap_ID, .data$stripped_cell_line_name, .data$disease, .data$lineage_subtype, .data$primary_or_metastasis) %>%
    dplyr::left_join(all_mutations_count_by_sample, by = "DepMap_ID") %>%
    dplyr::rename(Total_mutations = .data$n) %>%
    dplyr::mutate(Total_mutations = dplyr::case_when(
      is.na(.data$Total_mutations) ~ as.double(0),
      TRUE ~  as.double(.data$Total_mutations))) %>%
    dplyr::left_join(del_mutations_count_by_sample, by = "DepMap_ID") %>%
    dplyr::rename(Del_mutations = .data$n) %>%
    dplyr::mutate(Del_mutations = dplyr::case_when(
      is.na(.data$Del_mutations) ~ as.double(0),
      TRUE ~ as.double(.data$Del_mutations)),
      Del_hom_mut = dplyr::case_when(
        .data$DepMap_ID %in% hom_del_muts_by_sample$DepMap_ID ~ TRUE,
        TRUE ~ FALSE)) %>%
    dplyr::left_join(target_copy_num %>% 
                       dplyr::select(.data$DepMap_ID, .data$Status), by = "DepMap_ID") %>%
    dplyr::rename(CN_status = .data$Status) %>% 
    dplyr::arrange(-.data$Del_mutations, -.data$Total_mutations, .data$CN_status)
  
  # Annotate mutant group types based on conditions
  if(!all(summary$CN_status == "Unknown")){
    Groups <- summary %>% 
      dplyr::mutate(
        Group = dplyr::case_when(
        ((.data$CN_status == "Deep_del") | (.data$Del_hom_mut == TRUE)) ~ paste0(Input_gene,"_HomDel"),
        ((.data$Del_mutations > 1) & (.data$CN_status == "Neutral")) |
          ((.data$Del_mutations == 1) & (.data$CN_status == "Loss")) ~ paste0(Input_gene,"_T-HetDel"),
        ((.data$Del_mutations == 1) & (.data$CN_status == "Neutral")) |
          ((.data$Del_mutations == 0) & (.data$CN_status == "Loss")) ~ paste0(Input_gene,"_HetDel"),
        ((.data$Total_mutations == 0) & (.data$CN_status == "Neutral")) ~ "Control",
        ((.data$Total_mutations == 0) & (.data$CN_status == "Amplified")) ~ "Amplified",
        ((.data$Del_mutations == 1) & (.data$CN_status == "Amplified")) ~ "Others",
        TRUE ~ "Others")) %>%
      dplyr::mutate(GeneNameID = Input_geneID,
             GeneName = Input_gene)
  } else {
    Groups <- summary %>% 
      dplyr::mutate(
        Group = dplyr::case_when(
        (.data$Del_hom_mut == TRUE) ~ paste0(Input_gene,"_HomDel"),
        (.data$Del_mutations > 1) ~ paste0(Input_gene,"_T-HetDel"),
        (.data$Del_mutations == 1) ~ paste0(Input_gene,"_HetDel"),
        (.data$Total_mutations == 0) ~ "Control",
        TRUE ~ "Others")) %>%
      dplyr::mutate(GeneNameID = Input_geneID,
             GeneName = Input_gene)
  }
  
  # Print quick summary of mutants found
  if(any(stringr::str_detect(Groups$Group, "HomDel|HetDel"))){
    print("Mutants found! Summary of mutant cell lines:")
    print(table(Groups$Group))
  } else {
    print(paste0("Mutants found!\nSummary of mutant cell lines: ", table(Groups$Group)))
  }
  return(Groups)
}
