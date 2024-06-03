options(tibble.width = Inf)

# Clean up data and package them in .RData according to usage type (ie dep comparison, etc.)
library(pacman)
p_load(tidyverse, janitor)
DepMap_dir <- "./DepMap/23Q4/"

# Essential genes -----------------------------------------------------------------------
# Load files
Achilles_common_essentials <- read_csv(paste0(DepMap_dir,"/raw_data/AchillesCommonEssentialControls.csv")) # Essential in all DepMap cell lines
Inferred_common_essentials <- read_csv(paste0(DepMap_dir,"/raw_data/CRISPRInferredCommonEssentials.csv")) %>% select(Gene = Essentials) # Inferred from all models

# Join list
essential_genes <- full_join(Achilles_common_essentials, Inferred_common_essentials) %>% distinct

dim(Achilles_common_essentials)
# [1] 1247    1
dim(Inferred_common_essentials)
# [1] 1552    1
dim(essential_genes)
# [1] 1831    1

# Total of 2089 essential genes
essential_genes <- essential_genes %>% mutate(
    GeneNames = str_split_fixed(Gene, " ", n = 2)[,1],
    GeneID = str_split_fixed(Gene, "\\(|\\)", n = 3)[,2],
    GeneNameID = paste(GeneNames, GeneID, sep = "_")) %>%
  select(GeneNameID, GeneNames, GeneID) %>%
  arrange(GeneNameID)

# Non-Essential genes -----------------------------------------------------------------------
nonessentials <- read_csv(paste0(DepMap_dir,"/raw_data/AchillesNonessentialControls.csv"))

nonessential_genes <- nonessentials %>% mutate(
    GeneNames = str_split_fixed(Gene, " ", n = 2)[,1],
    GeneID = str_split_fixed(Gene, "\\(|\\)", n = 3)[,2],
    GeneNameID = paste(GeneNames, GeneID, sep = "_")) %>%
  select(GeneNameID, GeneNames, GeneID) %>%
  arrange(GeneNameID) %>% distinct

dim(nonessential_genes)
# A tibble: 781 x 3

# Sample Annotation ----------------------------------------------------------------------
sample_annot <- read_csv(paste0(DepMap_dir,"/raw_data/Model.csv"), guess_max = 2000) %>% # sample anot
  rename(
    DepMap_ID = ModelID,
    stripped_cell_line_name = StrippedCellLineName,
    primary_or_metastasis = PrimaryOrMetastasis) %>% 
  mutate(
    disease = OncotreePrimaryDisease,
    disease_subtype = OncotreeSubtype) %>%
  select(DepMap_ID:stripped_cell_line_name, disease, disease_subtype, everything())

# Mut calls ----------------------------------------------------------------------
omics_sample_key <- read_csv(paste0(DepMap_dir,"/raw_data/OmicsProfiles.csv"),
  guess_max = 2000)

mut_calls <- read_csv(paste0(DepMap_dir,"/raw_data/OmicsSomaticMutationsProfile.csv"),
  guess_max = Inf) %>%
  left_join((omics_sample_key %>% select(ProfileID, DepMap_ID = ModelID))) %>%
  select(DepMap_ID, everything()) %>%
  rename(Hugo_Symbol = HugoSymbol) %>%
  mutate(chr = as.character(str_split_fixed(Chrom, "chr", 2)[,2]))

# This version denotes genotype differently from previous version. Now Requires a change in zygosity determination. GT == "1|1", indicates homozygous alteration. The others are all hets. 

# Copy number --------------------------------------------------------------------------
copy_num <- read_csv(paste0(DepMap_dir,"/raw_data/OmicsCNGene.csv")) %>%  # copy number
  rename(DepMap_ID = "...1")
copy_num_annot <- tibble(names = colnames(copy_num)) %>%
  mutate(GeneNames = str_split_fixed(names, " ", n = 2)[,1],
         GeneID = str_split_fixed(names, "\\(|\\)", n = 3)[,2],
         GeneNameID = paste(GeneNames, GeneID, sep = "_")) %>%
         .[-1,]

# Dependency ----------------------------------------------------------------------
dep <- read_csv(paste0(DepMap_dir,"/raw_data/CRISPRGeneDependency.csv")) # BROAD + SANGER

# Fix column names 
dep_col_names_raw <- colnames(dep)[-1] %>%
  str_split_fixed(., regex("\ |\\(|\\)"), 4) %>%
  as_tibble(.name_repair = "unique") %>%
  mutate(new_names = paste0(...1,"_",...3))
colnames(dep) <- c("DepMap_ID", dep_col_names_raw$new_names)

# Create annotation for matching
dep_annot <- tibble(names = colnames(dep)) %>%
  mutate(GeneNames = str_split_fixed(names, "_", n = 2)[,1],
         GeneID = str_split_fixed(names, "_", n = 2)[,2],
         GeneNameID = names) %>%
         .[-1,]

# Gene effects --------------------------------------------------------------------
gene_effect <- read_csv(paste0(DepMap_dir,"/raw_data/CRISPRGeneEffect.csv"))

# Fix column names 
gene_effect_col_names_raw <- colnames(gene_effect)[-1] %>%
  str_split_fixed(., regex("\ |\\(|\\)"), 4) %>%
  as_tibble(.name_repair = "unique") %>%
  mutate(new_names = paste0(...1,"_",...3))
colnames(gene_effect) <- c("DepMap_ID", gene_effect_col_names_raw$new_names)

# pivoting
gene_effect <- gene_effect %>%
  pivot_longer(-DepMap_ID, names_to = "GeneNameID", values_to = "Effect_score")

# GeneExp ---------------------------------------------------------------------------------
CCLE_exp_raw <- read_csv(paste0(DepMap_dir,"/raw_data/OmicsExpressionAllGenesTPMLogp1Profile.csv")) # gene expression protein coding genes only , log2 transformed + 1

# Fix column names 
CCLE_exp_col_names_raw <- colnames(CCLE_exp_raw)[-1] %>%
  str_split_fixed(., regex("\ |\\(|\\)"), 4) %>%
  as_tibble(.name_repair = "unique") %>%
  mutate(new_names = paste0(...1,"_",...3))
colnames(CCLE_exp_raw) <- c("ProfileID", CCLE_exp_col_names_raw$new_names)

CCLE_exp_temp <- CCLE_exp_raw %>% 
  left_join(omics_sample_key %>% select(ProfileID, DepMap_ID = ModelID)) %>%
  select(DepMap_ID, -ProfileID, everything())

# Check to see how different the expressions are 
test <- CCLE_exp_temp %>% filter(DepMap_ID == "ACH-000029") %>% 
  pivot_longer(-c(DepMap_ID, ProfileID), names_to = "Gene", values_to = "TPM") %>%
  select(ProfileID, Gene, TPM) %>%
  pivot_wider(names_from = ProfileID, values_from = TPM)

ggplot(test, aes(x = test$"PR-pOBrMJ", y = test$"PR-lqUArB"))+
  geom_point()

# t = 1391.5, df = 53959, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9861203 0.9865779
# sample estimates:
#      cor
# 0.986351

# Take mean for cell lines with multiple sequencing results
get_count <- CCLE_exp_temp %>% count(DepMap_ID)
single <- get_count %>% filter(n == 1)
multi <- get_count %>% filter(n > 1)

CCLE_exp_single <- CCLE_exp_temp %>% 
  filter(DepMap_ID %in% single$DepMap_ID) %>%
  pivot_longer(-c(DepMap_ID, ProfileID), names_to = "Gene", values_to = "TPM") %>%
  group_by(DepMap_ID, Gene) %>%
  summarise(mean_TPM = mean(TPM))

CCLE_exp_multi <- CCLE_exp_temp %>% 
  filter(DepMap_ID %in% multi$DepMap_ID) %>%
  pivot_longer(-c(DepMap_ID, ProfileID), names_to = "Gene", values_to = "TPM") %>%
  group_by(DepMap_ID, Gene) %>%
  summarise(mean_TPM = mean(TPM))

CCLE_exp <- bind_rows(CCLE_exp_single, CCLE_exp_multi) %>% 
  pivot_wider(names_from = Gene, values_from = mean_TPM)

# Create annotation for matching
CCLE_exp_annot <- tibble(names = colnames(CCLE_exp)) %>%
  mutate(GeneNames = str_split_fixed(names, "_", n = 2)[,1],
         GeneID = str_split_fixed(names, "_", n = 2)[,2],
         GeneNameID = names) %>%
         .[-1,]

# Proteomics --------------------------------------------------------------------------------
# This dataset has not been updated since 22Q2, just copy and paste it over. 

#### SAVE Data for GINI! ----------------------------------------------------------
GRETTA_dir_23Q4 <- "./DepMap/GRETTA_data/23Q4/data/"

# save data individually for GRETTA 
save(essential_genes, file = paste0(GRETTA_dir_23Q4, "essential_genes.rda"))
save(nonessential_genes, file = paste0(GRETTA_dir_23Q4, "nonessential_genes.rda"))
save(sample_annot, file = paste0(GRETTA_dir_23Q4, "sample_annot.rda"))
save(mut_calls, file = paste0(GRETTA_dir_23Q4, "mut_calls.rda"))
save(copy_num_annot, file = paste0(GRETTA_dir_23Q4, "copy_num_annot.rda"))
save(copy_num, file = paste0(GRETTA_dir_23Q4, "copy_num.rda"))
save(dep_annot, file = paste0(GRETTA_dir_23Q4, "dep_annot.rda"))
save(dep, file = paste0(GRETTA_dir_23Q4, "dep.rda"))
save(gene_effect, file = paste0(GRETTA_dir_23Q4, "gene_effect.rda"))
save(CCLE_exp_annot, file = paste0(GRETTA_dir_23Q4, "CCLE_exp_annot.rda"))
save(CCLE_exp, file = paste0(GRETTA_dir_23Q4, "CCLE_exp.rda"))
# Protein data: This dataset has not been updated since 22Q2, just copy it over. 