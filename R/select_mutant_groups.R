select_mutant_groups <- function(Input_gene){
  print(paste0("Selecting mutant groups for: ", Input_gene))
  
  
}








start_time <- Sys.time()
print(start_time)
# Run: Rscript4.0.2 Gather_all_gene_groups_pan_cancer.R > Rout/Gather_all_gene_groups_pan_cancer.Rout 2>&1 &
#Create a master list of all groups
library(pacman)
p_load(tidyverse, foreach, doParallel, doMC, broom)
registerDoMC(5)
DepMap_dir <- "/projects/marralab/ytakemon_prj/DepMap/20Q1/"

## Load files ------------------------------------------------------------------
load(paste0(DepMap_dir,"/Data/RData/20Q1_full.RData"))

AllGenes <- dep_annot$GeneNames[-1] #remove X1

res <- NULL
combine_by <- function(x, y){
  bind_rows(x, y)
}

# 1:length(AllGenes)
res <- foreach(i = 1:length(AllGenes), .combine = combine_by) %dopar% {
  
  print(paste0("Processing ", i, " of ",length(AllGenes)))
  
  Target_gene <- AllGenes[i]
  Target_ID <- dep_annot %>% filter(GeneNames == Target_gene) %>% pull(GeneNameID)
  
  #Copy number
  if(any(names(copy_num) %in% Target_ID)){
    target_copy_num <- copy_num %>%
      select(DepMap_ID, all_of(Target_ID)) %>%
      filter(DepMap_ID %in% dep$DepMap_ID) %>%
      arrange(DepMap_ID) %>%
      mutate( Status = case_when(
        !!as.name(Target_ID) <= 0.25 ~ "Deep_del",
        !!as.name(Target_ID) > 0.25 & !!as.name(Target_ID) < 0.75 ~ "Loss",
        !!as.name(Target_ID) >= 0.75 & !!as.name(Target_ID) < 1.25 ~ "Neutral",
        !!as.name(Target_ID) >= 1.25 ~ "Amplified",
        TRUE ~ "Other"))
    
  } else {
    target_copy_num <- dep %>% select(DepMap_ID) %>%
      arrange(DepMap_ID) %>%
      mutate(!!as.name(Target_ID) := NA,
             Status = "Unknown")
  }
  
  # ALL Mutations
  target_mut <- mut_calls %>%
    filter((DepMap_ID %in% dep$DepMap_ID) & (Hugo_Symbol %in% Target_gene)) %>%
    mutate(AC_combined = coalesce(CGA_WES_AC, SangerRecalibWES_AC, SangerWES_AC, RNAseq_AC,HC_AC, RD_AC, WGS_AC), #(Alt:REF)
           AC_ref_NULL = grepl(":0", AC_combined)) %>%
    mutate(AC_Variant = case_when(
      .$AC_ref_NULL == "TRUE" ~ "Hom_Mut",
      TRUE ~ "Het_Mut"),
      AC_Variant = paste0(Variant_Classification," ",AC_Variant)) %>% # are there any with 0 contribution from reference?
    select(Hugo_Symbol, Chromosome:Annotation_Transcript, cDNA_Change:COSMIChsCnt, Variant_annotation:AC_Variant) %>%
    arrange(Start_position) %>%
    select(DepMap_ID, everything())
  
  all_mutations_count_by_sample <- target_mut %>% count(DepMap_ID)
  
  # Del mutations
  mut_dels <- target_mut %>% filter(Variant_annotation == "damaging")
  num_del_samples <- mut_dels %>% select(DepMap_ID) %>% distinct()
  del_mutations_count_by_sample <- mut_dels %>% count(DepMap_ID)
  hom_del_muts_by_sample <- mut_dels %>% select(DepMap_ID, AC_ref_NULL) %>%
    distinct() %>%
    filter(AC_ref_NULL == TRUE)
  
  # Multi muts
  multi_mut_dels <- mut_dels %>%
    filter(AC_ref_NULL == FALSE) %>%
    add_count(DepMap_ID) %>% filter(n > 1)
  count_sample_multi_mut_dels <- multi_mut_dels %>% select(DepMap_ID) %>% distinct()
  
  # Summarise results
  summary <- sample_annot %>% filter(DepMap_ID %in% dep$DepMap_ID) %>%
    select(DepMap_ID, stripped_cell_line_name, disease, lineage_subtype, primary_or_metastasis) %>%
    left_join(. , all_mutations_count_by_sample, by = "DepMap_ID") %>%
    rename(Total_mutations = n) %>%
    mutate(Total_mutations = case_when(
      is.na(Total_mutations) ~ as.double(0),
      TRUE ~  as.double(Total_mutations))) %>%
    left_join(., del_mutations_count_by_sample, by = "DepMap_ID") %>%
    rename(Del_mutations = n) %>%
    mutate(Del_mutations = case_when(
      is.na(Del_mutations) ~ as.double(0),
      TRUE ~ as.double(Del_mutations)),
      Del_hom_mut = case_when(
        DepMap_ID %in% hom_del_muts_by_sample$DepMap_ID ~ TRUE,
        TRUE ~ FALSE)) %>%
    left_join(., target_copy_num %>% select(DepMap_ID, Status), by = "DepMap_ID") %>%
    rename(CN_status = Status) %>% arrange(-Del_mutations, -Total_mutations, CN_status)
  
  if(!all(summary$CN_status == "Unknown")){
    Groups <- summary %>% mutate(
      Group = case_when(
        ((CN_status == "Deep_del") | (Del_hom_mut == TRUE)) ~ paste0(Target_gene,"_mut_1"),
        ((Del_mutations > 1) & (CN_status == "Neutral")) |
          ((Del_mutations == 1) & (CN_status == "Loss")) ~ paste0(Target_gene,"_mut_2"),
        ((Del_mutations == 1) & (CN_status == "Neutral")) |
          ((Del_mutations == 0) & (CN_status == "Loss")) ~ paste0(Target_gene,"_mut_3"),
        ((Total_mutations == 0) & (CN_status == "Neutral")) ~ "Control",
        ((Total_mutations == 0) & (CN_status == "Amplified")) ~ "Amplified",
        ((Del_mutations == 1) & (CN_status == "Amplified")) ~ "Others",
        TRUE ~ "Others")) %>%
      mutate(Group = fct_relevel(Group, "Control", "Amplified", paste0(Target_gene,"_mut_1"),  paste0(Target_gene,"_mut_2"), paste0(Target_gene,"_mut_3"), "Others")) %>%
      mutate(GeneNameID = Target_ID,
             GeneName = Target_gene)
  } else {
    Groups <- summary %>% mutate(
      Group = case_when(
        (Del_hom_mut == TRUE) ~ paste0(Target_gene,"_mut_1"),
        (Del_mutations > 1) ~ paste0(Target_gene,"_mut_2"),
        (Del_mutations == 1) ~ paste0(Target_gene,"_mut_3"),
        (Total_mutations == 0) ~ "Control",
        TRUE ~ "Others")) %>%
      mutate(Group = fct_relevel(Group, "Control", paste0(Target_gene,"_mut_1"),  paste0(Target_gene,"_mut_2"), paste0(Target_gene,"_mut_3"), "Others")) %>%
      mutate(GeneNameID = Target_ID,
             GeneName = Target_gene)
  }
}

res %>%
  relocate(c(GeneNameID, GeneName)) %>%
  mutate(disease = str_replace_all(disease, " |/","_")) %>%
  write_csv(file = paste0(DepMap_dir, "Analysis/Automate_single_screen/Pan_Cancer/AllGene_Groups_pan_cancer.csv"))

# Tally pan cancer
res %>%
  mutate(disease = str_replace_all(disease, " |/","_")) %>%
  group_by(GeneNameID, GeneName, Group) %>% tally %>%
  write_csv(.,
            file = paste0(DepMap_dir, "Analysis/Automate_single_screen/Pan_Cancer/AllGene_Groups_tally_pan_cancer.csv"))

# Tally by disease type ---------------------------------------------------------
# Unique disease types
sample_annot_used <- sample_annot %>%
  filter(DepMap_ID %in% dep$DepMap_ID) %>%
  mutate(disease = str_replace_all(disease, " |/","_"))

# less than 3: insufficient group size
disease_insuff <- sample_annot_used %>% group_by(disease) %>% tally %>%
  arrange(disease) %>% filter(n < 3) %>% pull(disease)

# greater or equal to 3: sufficient group size
disease_suff <- sample_annot_used %>% group_by(disease) %>% tally %>%
  arrange(disease) %>% filter(n >= 3) %>% pull(disease)

# Create directory for disease if not already available
for(disease_type in disease_suff){
  if(dir.exists(paste0(DepMap_dir,"/Analysis/Automate_single_screen/Per_Disease/",disease_type)) == FALSE){
    dir.create(paste0(DepMap_dir,"/Analysis/Automate_single_screen/Per_Disease/",disease_type))
  } else {
    print(paste("Directory exists for:", disease_type))
  }
}

res %>%
  mutate(disease = str_replace_all(disease, " |/","_")) %>%
  group_by(GeneNameID, GeneName, disease, Group) %>% tally %>%
  write_csv(.,
            file = paste0(DepMap_dir, "Analysis/Automate_single_screen/Per_Disease/AllGene_Groups_tally_per_disease.csv"))

end_time <- Sys.time()
print(paste("Start time", start_time))
print(paste("End time", end_time))
print(end_time - start_time)
