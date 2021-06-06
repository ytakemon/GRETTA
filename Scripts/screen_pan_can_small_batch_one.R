# Batch screening adapted from: 20Q1_par_MWU_genetic_screen.R
batch_list <- commandArgs(trailingOnly = TRUE)
cat(paste("Scanning", length(batch_list), "genes", "\n"))
print(batch_list)
start_time <- Sys.time()

# Set up------------------------------------------------------
library(pacman)
p_load(tidyverse, rcompanion, doMC, broom, diptest)
registerDoMC(10)
DepMap_dir <- "/projects/marralab/ytakemon_prj/DepMap/20Q1/"

## Load files ---------------------------------------------------------
load(paste0(DepMap_dir,"/Data/RData/20Q1_minimum_screening_data.RData"))

Groups <- read_csv(paste0(DepMap_dir, "Analysis/Automate_single_screen/Pan_Cancer/AllGene_Groups_pan_cancer.csv")) %>%
  filter(GeneName %in% batch_list)

Groups_tally <- read_csv(paste0(DepMap_dir, "Analysis/Automate_single_screen/Pan_Cancer/AllGene_Groups_tally_pan_cancer.csv")) %>%
    filter(GeneName %in% batch_list)

# Perpare Data -----------------------------------------------------------
# Remove essential data
dep <- dep %>%
  pivot_longer(cols = matches("\\d"), names_to = "GeneNameID", values_to = "DepProb") %>%
  filter(!(GeneNameID %in% essential_genes$GeneNameID),
    !(GeneNameID %in% nonessential_genes$GeneNameID))

# function for forloop
combine_all_by <- function(x, y){
  bind_rows(x, y)
}

# Begin single gene analysis -----------------------------------------------
batch_list <- as.character(batch_list)
for(i in 1:length(batch_list)){
  Target_gene <- batch_list[i]
  # Target_gene <- "CIC"
  if(!any(Target_gene %in% Groups_tally$GeneName)){
    # if this gene doesn't exist
    cat(paste0("Warning: ", Target_gene," not in DepMap. Is the spelling correct?\n"))
    next
  }

  # Create Directory if not present ----------------------------------------------
  if(dir.exists(paste0(DepMap_dir, "Analysis/Automate_single_screen/Pan_Cancer", "/", Target_gene))){
    print(paste0("Directory for ", Target_gene, " exists"))
  } else{
    print(paste0("Directory for ", Target_gene, " does not exist. Creating new directory ... "))
    dir.create(paste0(Report_dir, "/", Target_gene))
  }

  Target_ID <- Groups_tally %>% select(GeneNameID, GeneName) %>%
    distinct() %>%
    filter(GeneName == !!Target_gene) %>%
    pull(GeneNameID)

  Control_group_avail <- Groups_tally %>%
    filter(GeneName == !!Target_gene,
      str_detect(Group, "Control"),
      n >= 2) %>% pull(Group)

  Mutant_groups_avail <- Groups_tally %>%
    filter(GeneName == !!Target_gene,
      str_detect(Group, "_mut_|Amp"),
      n >= 2) %>% pull(Group)

  # Check point --------------------------------------------------------
  if(length(Control_group_avail) == 0){
    cat(paste0("No control group found: ", Target_gene,"\n"))
    next
  } else {
    Control_group <- Control_group_avail
    Control_samples <- Groups %>%
      filter(GeneName == !!Target_gene,
        Group == !!Control_group) %>%
      pull(DepMap_ID)
  }

  if(length(Mutant_groups_avail) == 0){
    cat(paste0("No mutant groups found: ", Target_gene,"\n"))
    next # if no mutants are found go to next gene
  }

  # Begin nested loop
  # 1:length(unique(dep$GeneNameID))
  All_res <- NULL
  All_res <- foreach(Mutant_group = Mutant_groups_avail, .combine = bind_rows) %dopar% {
    foreach(each = 1:length(unique(dep$GeneNameID)), .combine = bind_rows) %dopar% {
      if(each == 1){
        cat(paste0("Processing ", each, " of ", length(unique(dep$GeneNameID))),"\n")
        cat(paste0("Running ", Target_gene, " mutant group: ", Mutant_group, "\n"))
      } else if(each == length(unique(dep$GeneNameID))){
        cat(paste0("Processing ", each, " of ", length(unique(dep$GeneNameID))),"\n")
      } else if(each%%1000 == 0){
        cat(paste0("Processing ", each, " of ", length(unique(dep$GeneNameID))),"\n")
      }

      # Get samples in group
      Mutant_samples <- Groups %>%
        filter(GeneName == !!Target_gene,
          Group == !!Mutant_group) %>%
        pull(DepMap_ID)

      select_dep <- dep %>%
        mutate(CellType = case_when(
          DepMap_ID %in% Mutant_samples ~ "Mutant",
          DepMap_ID %in% Control_samples ~ "Control",
          TRUE ~ "Others")) %>%
        filter(CellType != "Others") %>%
        mutate(CellType = fct_relevel(CellType, "Control", "Mutant"))

      # Get each gene
      geneID <- unique(select_dep$GeneNameID)[each]
      df <- select_dep %>% filter(GeneNameID == geneID) %>%
        filter(!is.na(DepProb))

      # Begin analysis
      if(all(df$DepProb == 0)){
        populate <- rep(0,11)
      # } else if(all(df$DepProb < 0.01)) {
      #   populate <- rep(0,11)
      } else {
        # # MWU doesn't handle na or zero's well so
        # # FOR NOW remove zeros.
        # df <- df %>% filter(!is.na(DepProb)) %>%
        #   filter(DepProb != 0)

        stats <- df %>%
          group_by(CellType) %>%
          summarize(Median = median(DepProb, na.rm = TRUE),
                    Mean = mean(DepProb, na.rm = TRUE),
                    SD = sd(DepProb, na.rm = TRUE),
                    IQR = IQR(DepProb, na.rm = TRUE),
                   .groups = "drop")

        if((any(is.na(stats)) != TRUE) & (nrow(stats) == 2)){

          fit_pval <- wilcox.test(DepProb ~ CellType, df,
                             paired = F,
                             alternative = "two.sided",
                             conf.int = T,
                             na.action = "na.omit")$p.value

          # If group size is < 3 cliffDelta will have error:
          # missing value where TRUE/FALSE needed
          # Important note: because celltype has a specific factor order specified for select_dep, a delta > 0 indicates an effect score greater in the Control (first level) and a delta < 0 means an effect score greater in the Mutant( second level)

          CliffDelta <- cliffDelta(DepProb ~ CellType, df)

          # Add diptest for uni/multi modality
          # null hypothesis if p > 0.05 the data is unimodal
          # alternative hyp if p < 0.05 the data is multimodal
          dip_pval <- df %>%
            filter(CellType == "Mutant") %>%
            pull(DepProb) %>%
            dip.test(.) %>% tidy %>%
            pull(p.value)
        } else if((any(is.na(stats)) == TRUE) & (nrow(stats) == 2)){
          populate <- rep(0,11)
        }

        if((any(is.na(stats)) != TRUE) & (nrow(stats) == 2)){
          populate <- as.numeric(c(unlist(stats)[-c(1,2)], fit_pval, -CliffDelta[[1]], dip_pval))
        } else {
          populate <- rep(0,11)
        }
      }

      tibble(Result = c("Control_median", "Mutant_median", "Control_mean", "Mutant_mean","Control_sd", "Mutant_sd", "Control_iqr","Mutant_iqr", "Pval",
        "CliffDelta", "dip_pval")) %>%
          mutate(!!sym(geneID) := populate) %>%
          pivot_longer(-Result) %>%
          pivot_wider(names_from = Result, values_from = value) %>%
          rename(GeneNameID = name) %>%
          mutate(Mutant_group = !!Mutant_group,
                 Control_group = !!Control_group)
    }} # End of for loop

    # Add mutant group name
  output <- All_res %>%
    mutate(log2FC_by_median = log2(Mutant_median / Control_median),
           log2FC_by_mean = log2(Mutant_mean / Control_mean)) %>%
    left_join(., dep_annot %>% select(GeneNameID, GeneNames),
        by = "GeneNameID") %>%
    select(Control_group, Mutant_group, GeneNameID, GeneNames, Control_median:Pval, log2FC_by_median, log2FC_by_mean, everything())

  # output
  write_csv(output,
      file = paste0(DepMap_dir, "Analysis/Automate_single_screen/Pan_Cancer/",Target_gene, "/", Target_gene,"_MWU_genetic_screening_results.csv"))
} # End of this Target_gene

end_time <- Sys.time()
cat(paste("Start time", start_time,"\n"))
cat(paste("End time", end_time,"\n"))
cat(end_time - start_time,"\n")

# did you fix 1:length(unique(dep$GeneNameID)) ?
