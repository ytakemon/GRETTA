#' Select mutant groups based on input gene of interest
#' 
#' @description
#' `select_mutant_groups()` assigns cancer cell lines to either `Control` groups or one of the following mutant groups: `HomDel`,
#' `T-HetDel`, `HetDel`, `Amplified`, or `Others` (see below for more details).
#' \describe{
#'  \item `Control` cell lines do not harbor any single nucleotide variations (SNVs) or insertions and deletions (InDels) with a neutral copy number (CN).
#'  \item `HomDel` cell lines harbor one or more homozygous deleterious SNVs or have deep CN loss.
#'  \item `T-HetDel` cell lines harbor two or more heterozygous deleterious SNVs/InDels with neutral or CN loss.
#'  \item `HetDel` cell lines harbor one heterozygous deleterious SNV/InDel with neutral CN, or no SNV/InDel with CN loss. 
#'  \item `Amplified` cell lines harbor no SNVs/InDels with increased CN.
#'  \item `Others` cell lines harbor deleterious SNVs with increased CN.
#' }
#' 
#' @param Input_gene string Input Hugo Symbol
#' @return data frame containing a summary of mutations found in cell lines and their control and mutant group assignments.
#' 
#' @export
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
  
  # Check if input gene mutations exist
  if(!any(mut_calls$Hugo_Symbol %in% Input_gene)|!any(copy_num_annot$GeneNames %in% Input_gene)){
    stop(paste0("No mutations were found for: ", Input_gene,"\nPlease check spelling and for valid Hugo Symbols"))
  }
  # Convert to unique geneID
  Input_geneID <- get_GeneNameID(Input_gene)

  # Get copy number
  if(any(names(copy_num) %in% Input_geneID)){
    target_copy_num <- copy_num %>%
      dplyr::select(DepMap_ID, dplyr::all_of(Input_geneID)) %>%
      dplyr::filter(DepMap_ID %in% dep$DepMap_ID) %>%
      dplyr::arrange(DepMap_ID) %>%
      dplyr::mutate( Status = dplyr::case_when(
        !!as.name(Input_geneID) <= 0.25 ~ "Deep_del",
        !!as.name(Input_geneID) > 0.25 & !!as.name(Input_geneID) < 0.75 ~ "Loss",
        !!as.name(Input_geneID) >= 0.75 & !!as.name(Input_geneID) < 1.25 ~ "Neutral",
        !!as.name(Input_geneID) >= 1.25 ~ "Amplified",
        TRUE ~ "Other"))
    
  } else {
    target_copy_num <- dep %>% 
      dplyr::select(DepMap_ID) %>%
      dplyr::arrange(DepMap_ID) %>%
      dplyr::mutate(!!as.name(Input_geneID) := NA,
             Status = "Unknown")
  }
  
  # Get ALL Mutations
  target_mut <- mut_calls %>%
    dplyr::filter((DepMap_ID %in% dep$DepMap_ID) & (Hugo_Symbol %in% Input_gene)) %>%
    dplyr::mutate(AC_combined = coalesce(CGA_WES_AC, SangerRecalibWES_AC, SangerWES_AC, RNAseq_AC,HC_AC, RD_AC, WGS_AC), #(Alt:REF)
           AC_ref_NULL = grepl(":0", AC_combined)) %>%
    dplyr::mutate(AC_Variant = case_when(
      .$AC_ref_NULL == "TRUE" ~ "Hom_Mut",
      TRUE ~ "Het_Mut"),
      AC_Variant = paste0(Variant_Classification," ",AC_Variant)) %>% # are there any with 0 contribution from reference?
    dplyr::select(Hugo_Symbol, Chromosome:Annotation_Transcript, cDNA_Change:COSMIChsCnt, Variant_annotation:AC_Variant) %>%
    dplyr::arrange(Start_position) %>%
    dplyr::select(DepMap_ID, everything())
  # Count number of mutations found per sample
  all_mutations_count_by_sample <- target_mut %>% dplyr::count(DepMap_ID)
  
  # Get deleterious/damaging mutations
  mut_dels <- target_mut %>% dplyr::filter(Variant_annotation == "damaging")
  # Count number of deleterious mutations found per sample
  del_mutations_count_by_sample <- mut_dels %>% dplyr::count(DepMap_ID)
  # Find samples with homozygous deleterious mutations (HomDels)
  hom_del_muts_by_sample <- mut_dels %>% dplyr::select(DepMap_ID, AC_ref_NULL) %>%
    dplyr::distinct() %>%
    dplyr::filter(AC_ref_NULL == TRUE)
  
  # Find samples with multiple heterozygous deleterious mutations. Trans-heterozygous mutants (T-HetDels)
  multi_mut_dels <- mut_dels %>%
    dplyr::filter(AC_ref_NULL == FALSE) %>%
    dplyr::add_count(DepMap_ID) %>% 
    dplyr::filter(n > 1)
  
  # Summarize mutations for all samples 
  summary <- sample_annot %>% 
    dplyr::filter(DepMap_ID %in% dep$DepMap_ID) %>%
    dplyr::select(DepMap_ID, stripped_cell_line_name, disease, lineage_subtype, primary_or_metastasis) %>%
    dplyr::left_join(. , all_mutations_count_by_sample, by = "DepMap_ID") %>%
    dplyr::rename(Total_mutations = n) %>%
    dplyr::mutate(Total_mutations = case_when(
      is.na(Total_mutations) ~ as.double(0),
      TRUE ~  as.double(Total_mutations))) %>%
    dplyr::left_join(., del_mutations_count_by_sample, by = "DepMap_ID") %>%
    dplyr::rename(Del_mutations = n) %>%
    dplyr::mutate(Del_mutations = case_when(
      is.na(Del_mutations) ~ as.double(0),
      TRUE ~ as.double(Del_mutations)),
      Del_hom_mut = case_when(
        DepMap_ID %in% hom_del_muts_by_sample$DepMap_ID ~ TRUE,
        TRUE ~ FALSE)) %>%
    dplyr::left_join(., target_copy_num %>% select(DepMap_ID, Status), by = "DepMap_ID") %>%
    dplyr::rename(CN_status = Status) %>% arrange(-Del_mutations, -Total_mutations, CN_status)
  
  # Annotate mutant group types based on conditions
  if(!all(summary$CN_status == "Unknown")){
    Groups <- summary %>% 
      dplyr::mutate(
        Group = case_when(
        ((CN_status == "Deep_del") | (Del_hom_mut == TRUE)) ~ paste0(Input_gene,"_HomDel"),
        ((Del_mutations > 1) & (CN_status == "Neutral")) |
          ((Del_mutations == 1) & (CN_status == "Loss")) ~ paste0(Input_gene,"_T-HetDel"),
        ((Del_mutations == 1) & (CN_status == "Neutral")) |
          ((Del_mutations == 0) & (CN_status == "Loss")) ~ paste0(Input_gene,"_HetDel"),
        ((Total_mutations == 0) & (CN_status == "Neutral")) ~ "Control",
        ((Total_mutations == 0) & (CN_status == "Amplified")) ~ "Amplified",
        ((Del_mutations == 1) & (CN_status == "Amplified")) ~ "Others",
        TRUE ~ "Others")) %>%
      dplyr::mutate(GeneNameID = Input_geneID,
             GeneName = Input_gene)
  } else {
    Groups <- summary %>% 
      dplyr::mutate(
        Group = case_when(
        (Del_hom_mut == TRUE) ~ paste0(Input_gene,"_HomDel"),
        (Del_mutations > 1) ~ paste0(Input_gene,"_T-HetDel"),
        (Del_mutations == 1) ~ paste0(Input_gene,"_HetDel"),
        (Total_mutations == 0) ~ "Control",
        TRUE ~ "Others")) %>%
      dplyr::mutate(GeneNameID = Input_geneID,
             GeneName = Input_gene)
  }
  
  # Print quick summary of mutants found
  if(any(str_detect(Groups$Group, "HomDel|HetDel"))){
    print("Mutants found! Summary of mutant cell lines:")
    print(table(Groups$Group))
  } else {
    print(paste0("Mutants found!\nSummary of mutant cell lines: ", table(Groups$Group)))
  }
  return(Groups)
}
