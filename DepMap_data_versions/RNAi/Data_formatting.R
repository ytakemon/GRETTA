options(tibble.width = Inf)

# For using RNAi data in GRETTA --------------------------------------------------
library(pacman)
p_load(tidyverse, janitor)
RNAi_dir <- "./DepMap/RNAi/"
DepMap_data_dir <- "./DepMap/23Q4/data/" # Eg. using 23Q4

# Load files
load(paste0(DepMap_data_dir,"dep_annot.rda")) 
load(paste0(DepMap_data_dir,"sample_annot.rda"))  
rnai <- read_csv(paste0(RNAi_dir,"raw_data/D2_combined_gene_dep_scores.csv"))
rnai_sample <- read_csv(paste0(RNAi_dir,"raw_data/sample_info.csv")) %>% rename(CCLEName = "CCLE_ID")

# Wrangle and format --------------------------------------------------
# Match CCLE ID to DepMap_ID 
rnai_sample_match <- rnai_sample %>% 
  left_join(sample_annot %>% select(CCLEName, DepMap_ID), by = "CCLEName")

# Check how many are matched and rescue unmatched lines
sample_annot <- sample_annot %>% mutate( CCLEName_split = str_split_fixed(CCLEName, "_", 2)[,1])
matched <- rnai_sample_match %>% filter(!is.na(DepMap_ID))
unmatched <- rnai_sample_match %>% filter(is.na(DepMap_ID)) %>% select(-DepMap_ID) %>% distinct %>%
  mutate(CCLEName_split = str_split_fixed(CCLEName, "_", 2)[,1]) %>% 
  left_join(sample_annot %>% select(CCLEName_split, DepMap_ID), by = "CCLEName_split") %>%
  filter(!is.na(DepMap_ID)) # 5 lines unmatchable.
rnai_sample_match <- bind_rows(matched, unmatched) %>% select(DepMap_ID, CCLEName)

# Wrangle RNAi data into similar format at dep data 
rnai_long <- rnai %>% rename(GeneID = 1) %>%
  pivot_longer(-GeneID, names_to = "CCLEName", values_to = "values") %>%
  left_join(rnai_sample_match, by = "CCLEName") %>% distinct %>% 
  filter(!is.na(DepMap_ID)) %>% # remove unmatchable lines
  mutate(GeneNames = str_split_fixed(GeneID, "\\ \\(", n = 2)[,1]) %>%
  left_join(dep_annot %>% select(GeneNames, GeneNameID), by = "GeneNames") %>% distinct %>%
  select(DepMap_ID, GeneNameID, values) %>% distinct %>%
  filter(!is.na(values), !is.na(GeneNameID))

rnai_annot <- rnai_long %>% select(GeneNameID) %>% distinct %>% 
  left_join(dep_annot %>% select(GeneNameID, GeneNames), by = "GeneNameID") %>% distinct

# Save data ---------------------------------------------------------------
save(rnai_long, file = paste0(RNAi_dir, "data/rnai_long.rda"))
save(rnai_annot, file = paste0(RNAi_dir, "data/rnai_annot.rda"))

# Move data to where the rest of your DepMap data are stored (eg. DepMap_data_dir)
file.copy(
  from = paste0(RNAi_dir, "data/rnai_long.rda"), 
  to = paste0(DepMap_data_dir, "data/rnai_long.rda")
  )

file.copy(
  from = paste0(RNAi_dir, "data/rnai_annot.rda"), 
  to = paste0(DepMap_data_dir, "data/rnai_annot.rda")
  )