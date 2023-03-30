options(tibble.width = Inf)

# Clean up data and package them in .RData according to usage type (ie dep comparison, etc.)
library(pacman)
p_load(tidyverse, janitor)
DepMap_dir <- "/projects/marralab/ytakemon_prj/DepMap/21Q4/"

# Essential genes -----------------------------------------------------------------------
# Load files
Achilles_common_essentials <- read_csv(paste0(DepMap_dir,"/Data/Achilles_common_essentials.csv")) # Essential in all DepMap cell lines
common_essentials <- read_csv(paste0(DepMap_dir,"/Data/common_essentials.csv")) # essential genes defined by Biomen (2014) and Hart (2015)

# Join list
essential_genes <- full_join(Achilles_common_essentials, common_essentials) %>% distinct

dim(Achilles_common_essentials)
#[1] 1868   1
dim(common_essentials)
#[1] 1247    1
dim(essential_genes)
#[1] 2089    1

# Total of 2089 essential genes
essential_genes <- essential_genes %>% mutate(
    GeneNames = str_split_fixed(gene, " ", n = 2)[,1],
    GeneID = str_split_fixed(gene, "\\(|\\)", n = 3)[,2],
    GeneNameID = paste(GeneNames, GeneID, sep = "_")) %>%
  select(GeneNameID, GeneNames, GeneID) %>%
  arrange(GeneNameID)

# Non-Essential genes -----------------------------------------------------------------------
nonessentials <- read_csv(paste0(DepMap_dir,"/Data/nonessentials.csv"))

nonessential_genes <- nonessentials %>% mutate(
    GeneNames = str_split_fixed(gene, " ", n = 2)[,1],
    GeneID = str_split_fixed(gene, "\\(|\\)", n = 3)[,2],
    GeneNameID = paste(GeneNames, GeneID, sep = "_")) %>%
  select(GeneNameID, GeneNames, GeneID) %>%
  arrange(GeneNameID)

dim(nonessential_genes)
# A tibble: 781 x 3

# Sample Annotation ----------------------------------------------------------------------
sample_annot <- read_csv(paste0(DepMap_dir,"/Data/sample_info.csv"), guess_max = 2000) %>% # sample anot
  rename(
    disease = primary_disease,
    disease_subtype = Subtype
  )

# Mut calls ----------------------------------------------------------------------
# load(paste0(DepMap_dir,"/Data/RData/20Q1_full.RData"))
mut_calls <- read_csv(paste0(DepMap_dir,"/Data/CCLE_mutations.csv"),
  guess_max = 2000)

# Copy number --------------------------------------------------------------------------
copy_num <- read_csv(paste0(DepMap_dir,"/Data/CCLE_gene_cn.csv")) %>%  # copy number
  rename(DepMap_ID = "...1")
copy_num_annot <- tibble(names = colnames(copy_num)) %>%
  mutate(GeneNames = str_split_fixed(names, " ", n = 2)[,1],
         GeneID = str_split_fixed(names, "\\(|\\)", n = 3)[,2],
         GeneNameID = paste(GeneNames, GeneID, sep = "_"))

# Dependency ----------------------------------------------------------------------
# https://forum.depmap.org/t/how-comparable-are-the-achilles-and-crispr-gene-dependency-files/1410/7
# current recommendation is to use the integrated dependency data
dep <- read_csv(paste0(DepMap_dir,"/Data/CRISPR_gene_dependency.csv")) # BROAD + SANGER

# Fix column names 
dep_col_names_raw <- colnames(dep)[-1] %>%
  str_split_fixed(., regex("\ |\\(|\\)"), 4) %>%
  as_tibble() %>%
  mutate(new_names = paste0(V1,"_",V3))
colnames(dep) <- c("DepMap_ID", dep_col_names_raw$new_names)

# Create annotation for matching
dep_annot <- tibble(names = colnames(dep)) %>%
  mutate(GeneNames = str_split_fixed(names, "_", n = 2)[,1],
         GeneID = str_split_fixed(names, "_", n = 2)[,2],
         GeneNameID = names)

dim(dep)
dim(dep_annot)

# Gene effects --------------------------------------------------------------------
gene_effect <- read_csv(paste0(DepMap_dir,"Data/CRISPR_gene_effect.csv"))

# Fix column names 
gene_effect_col_names_raw <- colnames(gene_effect)[-1] %>%
  str_split_fixed(., regex("\ |\\(|\\)"), 4) %>%
  as_tibble() %>%
  mutate(new_names = paste0(V1,"_",V3))

# pivoting
colnames(gene_effect) <- c("DepMap_ID", gene_effect_col_names_raw$new_names)
gene_effect <- gene_effect %>%
  pivot_longer(-DepMap_ID, names_to = "GeneNameID", values_to = "Effect_score")

# GeneExp ---------------------------------------------------------------------------------
CCLE_exp <- read_csv(paste0(DepMap_dir,"/Data/CCLE_expression.csv")) # gene expression protein coding genes only , log2 transformed + 1

# Fix column names 
CCLE_exp_col_names_raw <- colnames(CCLE_exp)[-1] %>%
  str_split_fixed(., regex("\ |\\(|\\)"), 4) %>%
  as_tibble() %>%
  mutate(new_names = paste0(V1,"_",V3))
colnames(CCLE_exp) <- c("DepMap_ID", CCLE_exp_col_names_raw$new_names)

# Create annotation for matching
CCLE_exp_annot <- tibble(names = colnames(CCLE_exp)) %>%
  mutate(GeneNames = str_split_fixed(names, "_", n = 2)[,1],
         GeneID = str_split_fixed(names, "_", n = 2)[,2],
         GeneNameID = names)

# Proteomics --------------------------------------------------------------------------------
protein <-  read_csv(paste0(DepMap_dir,"../20Q1/Data/protein_quant_current_normalized.csv")) # Proteomics
# need to make a conversion table to connect Gygi nomenclature to DepMap_ID

# Remove duplicates
protein_nodup <- protein %>%
  filter(!is.na(Gene_Symbol)) %>%
  mutate(Gene_human = paste0(Gene_Symbol,"_HUMAN")) %>%
  filter(Uniprot == Gene_human) %>%
  select(-Gene_human)

GetID <- function(x){
  if(x %in% sample_annot$CCLE_Name){
    y <- filter(sample_annot, CCLE_Name %in% x) %>% pull(DepMap_ID)
    return(y)
  } else {
    return(x)
  }
}

protein_annot <- tibble(GygiNames = colnames(protein)) %>%
  mutate(GygiNames_celllines = str_split_fixed(GygiNames, "_TenPx", n = 2)[,1],
         DepMap_ID = map_chr(GygiNames_celllines, GetID)) %>%
  select(GygiNames, DepMap_ID)

#### SAVE Data for GINI! ----------------------------------------------------------
GRETA_dir_21Q4 <- "/projects/marralab/ytakemon_prj/DepMap/GRETA_data/21Q4/data/"

# save data individually for GRETA 
save(CCLE_exp_annot, file = paste0(GRETA_dir_21Q4, "CCLE_exp_annot.rda"))
save(CCLE_exp, file = paste0(GRETA_dir_21Q4, "CCLE_exp.rda"))
save(copy_num_annot, file = paste0(GRETA_dir_21Q4, "copy_num_annot.rda"))
save(copy_num, file = paste0(GRETA_dir_21Q4, "copy_num.rda"))
save(dep_annot, file = paste0(GRETA_dir_21Q4, "dep_annot.rda"))
save(dep, file = paste0(GRETA_dir_21Q4, "dep.rda"))
save(essential_genes, file = paste0(GRETA_dir_21Q4, "essential_genes.rda"))
save(gene_effect, file = paste0(GRETA_dir_21Q4, "gene_effect.rda"))
save(mut_calls, file = paste0(GRETA_dir_21Q4, "mut_calls.rda"))
save(nonessential_genes, file = paste0(GRETA_dir_21Q4, "nonessential_genes.rda"))
save(protein_annot, file = paste0(GRETA_dir_21Q4, "protein_annot.rda"))
save(protein_nodup, file = paste0(GRETA_dir_21Q4, "protein_nodup.rda"))
save(protein, file = paste0(GRETA_dir_21Q4, "protein.rda"))
save(sample_annot, file = paste0(GRETA_dir_21Q4, "sample_annot.rda"))