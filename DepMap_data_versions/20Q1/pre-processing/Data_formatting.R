options(tibble.width = Inf)

# Clean up data and package them in .RData according to usage type (ie dep comparison, etc.)
library(tidyverse)
DepMap_dir <- "/projects/marralab/ytakemon_prj/DepMap/20Q1/"

# Essential genes -----------------------------------------------------------------------
# Load files
Achilles_common_essentials <- read_csv(paste0(DepMap_dir,"/Data/Achilles_common_essentials.csv")) # Essential in all DepMap cell lines
common_essentials <- read_csv(paste0(DepMap_dir,"/Data/common_essentials.csv")) # essential genes defined by Biomen (2014) and Hart (2015)

# Join list
essential_genes <- full_join(Achilles_common_essentials, common_essentials)

#> dim(Achilles_common_essentials)
#[1] 2148    1
#> dim(common_essentials)
#[1] 1246    1
#> dim(essential_genes)
#[1] 2290    1

# Total of 2290 essential genes
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
# A tibble: 758 x 3

# Sample Annotation ----------------------------------------------------------------------
sample_annot <- read_csv(paste0(DepMap_dir,"/Data/sample_info.csv"), guess_max = 2000) # sample anot

# Mut calls ----------------------------------------------------------------------
# load(paste0(DepMap_dir,"/Data/RData/20Q1_full.RData"))
mut_calls <- read_csv(paste0(DepMap_dir,"/Data/CCLE_mutations.csv"),
  col_types = "dcddcddccccccccccccclldldcccccccccc") %>%
  select(-X1) # mutation called

# # test - Add SIFT
# pacman::p_load(PolyPhen.Hsapiens.dbSNP131, SIFT.Hsapiens.dbSNP132)
# cic_muts <- mut_calls %>% filter(Hugo_Symbol == "CIC")
# cic_mut_wrsid <- cic_muts %>% filter(!is.na(dbSNP_RS))
#
#
# sil <- "rs371965355"
# dam <- "rs571341549"
#
# rsids <- cic_mut_wrsid$dbSNP_RS
# select(SIFT.Hsapiens.dbSNP132, keys=rsids)
#
#
# select(PolyPhen.Hsapiens.dbSNP131, keys=rsids)


# Copy number --------------------------------------------------------------------------
copy_num <- read_csv(paste0(DepMap_dir,"/Data/CCLE_gene_cn.csv")) # copy number, better
copy_num_annot <- tibble(names = colnames(copy_num)) %>%
  mutate(GeneNames = str_split_fixed(names, " ", n = 2)[,1],
         GeneID = str_split_fixed(names, "\\(|\\)", n = 3)[,2],
         GeneNameID = paste(GeneNames, GeneID, sep = "_")) %>%
  mutate(GeneNameID = case_when(.$GeneNameID == "X1_" ~ "DepMap_ID",
                              TRUE ~ .$GeneNameID))
colnames(copy_num) <- copy_num_annot$GeneNameID

# Dependency ----------------------------------------------------------------------
dep <- read_csv(paste0(DepMap_dir,"/Data/Achilles_gene_dependency.csv")) # dependency
dep_annot <- tibble(names = colnames(dep)) %>%
  mutate(GeneNames = str_split_fixed(names, " ", n = 2)[,1],
         GeneID = str_split_fixed(names, "\\(|\\)", n = 3)[,2],
         GeneNameID = paste(GeneNames, GeneID, sep = "_")) %>%
  mutate(GeneNameID = case_when(.$GeneNameID == "X1_" ~ "DepMap_ID",
                              TRUE ~ .$GeneNameID))
colnames(dep) <- dep_annot$GeneNameID

# GeneExp ---------------------------------------------------------------------------------
CCLE_exp <- read_csv(paste0(DepMap_dir,"/Data/CCLE_expression.csv")) # gene expression protein coding genes only , log2 transformed + 1
#CCLE_exp_full <- read_csv(paste0(DepMap_dir,"/Data/CCLE_expression_full.csv")) # gene expression
CCLE_exp_annot <- tibble(names = colnames(CCLE_exp)) %>%
  mutate(GeneNames = str_split_fixed(names, " ", n = 2)[,1],
         GeneID = str_split_fixed(names, "\\(|\\)", n = 3)[,2],
         GeneNameID = paste(GeneNames, GeneID, sep = "_")) %>%
  mutate(GeneNameID = case_when(.$GeneNameID == "X1_" ~ "DepMap_ID",
                              TRUE ~ .$GeneNameID))
colnames(CCLE_exp) <- CCLE_exp_annot$GeneNameID

# Proteomics --------------------------------------------------------------------------------
protein <-  read_csv(paste0(DepMap_dir,"/Data/protein_quant_current_normalized.csv")) # Proteomics
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

# Gene effects --------------------------------------------------------------------
gene_effect <- read_csv(paste0(DepMap_dir,"Data/Achilles_gene_effect.csv"))

gene_effect <- gene_effect %>%
  rename(DepMap_ID = X1) %>%
  pivot_longer(-DepMap_ID, names_to = "GeneNames_origin", values_to = "DepProb") %>%
  mutate(GeneNames_ID = str_replace_all(GeneNames_origin, " \\(", "_"),
         GeneNames_ID = str_replace_all(GeneNames_ID, "\\)", "")) %>%
  select(-GeneNames_origin) %>%
  pivot_wider(
    id_cols = DepMap_ID,
    names_from = GeneNames_ID,
    values_from = DepProb)


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