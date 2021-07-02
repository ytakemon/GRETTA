convert_GeneName_to_GeneNameID <- function(GeneName){
  res <- dep_annot %>%
    filter(GeneNames %in% GeneName) %>%
    pull(GeneNameID)
  return(res)
}

GetDepMapID <- function(x){
  y <- protein_annot %>% filter(GygiNames %in% x) %>% pull(DepMap_ID)
  return(y)
}

getAllExpr <- function(gene, uniprot_id = NULL){
  Target_gene <- gene
  target_uniprot <- uniprot_id

  target_ID <- dep_annot %>% filter(GeneNames == Target_gene) %>% pull(GeneNameID)

  extract_target_protein <- NULL
  extract_target_protein_supp <- NULL
  if(any(protein$Gene_Symbol %in% Target_gene)){
    extract_target_protein <- protein %>%
    filter(Gene_Symbol %in% Target_gene) %>%
    select(-c(Protein_Id,Description, Group_ID), ends_with("_TenPx**")) %>%
    pivot_longer(., -c(Gene_Symbol,Uniprot, Uniprot_Acc)) %>%
    mutate(DepMap_ID = map_chr(name, GetDepMapID))

    extract_target_protein_supp <- extract_target_protein %>%
      filter(str_detect(DepMap_ID, "^TenPx"))
    extract_target_protein <- extract_target_protein %>%
      filter(str_detect(DepMap_ID, "^ACH-") & DepMap_ID %in% dep$DepMap_ID)

  if(!is.null(target_uniprot)){
    extract_target_protein <- filter(extract_target_protein, Uniprot_Acc == target_uniprot)
    extract_target_protein_supp <- filter(extract_target_protein_supp, Uniprot_Acc == target_uniprot)
  }}

  # Extract mRNA data if any
  extract_target_rna <- CCLE_exp %>%
    select(DepMap_ID, !!as.name(target_ID)) %>%
    filter(DepMap_ID %in% dep$DepMap_ID)

  # Add to summary
  new_group <- Groups %>% left_join(., extract_target_rna, by = "DepMap_ID") %>%
    rename(!!as.name(paste0(Target_gene,"_RNA")) := !!as.name(target_ID))

  if(!is.null(extract_target_protein)){
    Groups <-  left_join(Groups, extract_target_protein %>% select(DepMap_ID, value)) %>%
    rename(!!as.name(paste0(Target_gene,"_RNA")) := value)
  }

  return(new_group)
}

customHSD <- function(control, mut, df, type){
  # df <- Groups
  if(type == "protein"){
    df <- df %>% filter(!is.na(protein))
    fit <- aov(protein ~ Group, df)
  } else if(type == "RNA"){
    df <- df %>% filter(!is.na(RNA))
    fit <- aov(RNA ~ Group, df)
  }

  # R.3.6 uses term "compare", R4.0 uses term "contrast"
  if(version$major == "4"){
    res <- TukeyHSD(fit) %>% tidy() %>% filter(contrast %in% paste(mut, control, sep = "-"))
  } else if(version$major == "3"){
    res <- TukeyHSD(fit) %>% tidy() %>% filter(comparison %in% paste(mut, control, sep = "-"))
  }

  p <- res$adj.p.value
  stars <- if(nrow(res) == 0){
    "NA"
  } else if(p < 0.001){
    "***"
  } else if(p < 0.01){
    "**"
  } else if(p < 0.05){
    "*"
  } else {
    "NS"
  }

  stars
}

geom_signif_customHSD <- function(data, type){
  #data <- protein_data
  #type <- "protein"

  all_groups <- fct_relevel(as.character(unique(data$Group)), "Control","Amplified")
  all_groups <- as.character(levels(all_groups))

  #check
  # need type specified
  if(!(type %in% c("protein","RNA"))){
    print("type not specified")
    return()
  }
  # need to remove others
  if(any(all_groups == "Others")){
    stop("please filter out Group name Other")
  } else if(all(all_groups != "Control")){ #must have control
    return()
  } else if(!any(all_groups %in% c(paste0(Target_gene,"_mut_1"),paste0(Target_gene,"_mut_2"),paste0(Target_gene,"_mut_3")))){
    print("no mutant group detected")
    return()
  } else {
    # the good stuff.
    if(length(all_groups) == 5){ # all groups are found excluding Other
      return(geom_signif(na.rm = T, comparisons = list(
        c("Control", "Amplified"),
        c("Control", paste0(Target_gene,"_mut_1")),
        c("Control", paste0(Target_gene,"_mut_2")),
        c("Control", paste0(Target_gene,"_mut_3"))),
        annotations = c(
          customHSD("Control", "Amplified", data, type),
          customHSD("Control",  paste0(Target_gene,"_mut_1"), data, type),
          customHSD("Control",  paste0(Target_gene,"_mut_2"), data, type),
          customHSD("Control",  paste0(Target_gene,"_mut_3"), data, type)),
        map_signif_level = TRUE,
        step_increase = 0.2))
    } else if(length(all_groups) == 4){ # One group is missing
      if(!any(all_groups == "Amplified")) {
        return(geom_signif(na.rm = T, comparisons = list(
          c("Control", paste0(Target_gene,"_mut_1")),
          c("Control", paste0(Target_gene,"_mut_2")),
          c("Control", paste0(Target_gene,"_mut_3"))),
          annotations = c(
            customHSD("Control",  paste0(Target_gene,"_mut_1"), data, type),
            customHSD("Control",  paste0(Target_gene,"_mut_2"), data, type),
            customHSD("Control",  paste0(Target_gene,"_mut_3"), data, type)),
          map_signif_level = TRUE,
          step_increase = 0.2))
      } else if(!any(all_groups == paste0(Target_gene,"_mut_1"))){
        return(geom_signif(na.rm = T, comparisons = list(
          c("Control", "Amplified"),
          c("Control", paste0(Target_gene,"_mut_2")),
          c("Control", paste0(Target_gene,"_mut_3"))),
          annotations = c(
            customHSD("Control", "Amplified", data, type),
            customHSD("Control",  paste0(Target_gene,"_mut_2"), data, type),
            customHSD("Control",  paste0(Target_gene,"_mut_3"), data, type)),
          map_signif_level = TRUE,
          step_increase = 0.2))
      } else if(!any(all_groups == paste0(Target_gene,"_mut_2"))){
        return(geom_signif(na.rm = T, comparisons = list(
          c("Control", "Amplified"),
          c("Control", paste0(Target_gene,"_mut_1")),
          c("Control", paste0(Target_gene,"_mut_3"))),
          annotations = c(
            customHSD("Control", "Amplified", data, type),
            customHSD("Control",  paste0(Target_gene,"_mut_1"), data, type),
            customHSD("Control",  paste0(Target_gene,"_mut_3"), data, type)),
          map_signif_level = TRUE,
          step_increase = 0.2))
      } else if(!any(all_groups == paste0(Target_gene,"_mut_3"))){
        return(geom_signif(na.rm = T, comparisons = list(
          c("Control", "Amplified"),
          c("Control", paste0(Target_gene,"_mut_1")),
          c("Control", paste0(Target_gene,"_mut_2"))),
          annotations = c(
            customHSD("Control", "Amplified", data, type),
            customHSD("Control",  paste0(Target_gene,"_mut_1"), data, type),
            customHSD("Control",  paste0(Target_gene,"_mut_2"), data, type)),
          map_signif_level = TRUE,
          step_increase = 0.2))
      }
    } else if(length(all_groups) == 3){
      if(any(all_groups == "Amplified")){
        if(any(all_groups == paste0(Target_gene,"_mut_1"))){
          return(geom_signif(na.rm = T, comparisons = list(
            c("Control", "Amplified"),
            c("Control", paste0(Target_gene,"_mut_1"))),
            annotations = c(
              customHSD("Control", "Amplified", data, type),
              customHSD("Control",  paste0(Target_gene,"_mut_1"), data, type)),
            map_signif_level = TRUE,
            step_increase = 0.2))
        } else if(any(all_groups == paste0(Target_gene,"_mut_2"))){
          return(geom_signif(na.rm = T, comparisons = list(
            c("Control", "Amplified"),
            c("Control", paste0(Target_gene,"_mut_2"))),
            annotations = c(
              customHSD("Control", "Amplified", data, type),
              customHSD("Control",  paste0(Target_gene,"_mut_2"), data, type)),
            map_signif_level = TRUE,
            step_increase = 0.2))
        } else if(any(all_groups == paste0(Target_gene,"_mut_3"))){
          return(geom_signif(na.rm = T, comparisons = list(
            c("Control", "Amplified"),
            c("Control", paste0(Target_gene,"_mut_3"))),
            annotations = c(
              customHSD("Control", "Amplified", data, type),
              customHSD("Control",  paste0(Target_gene,"_mut_3"), data, type)),
            map_signif_level = TRUE,
            step_increase = 0.2))
        }
      } else if(any(all_groups == paste0(Target_gene,"_mut_1"))){
        if(any(all_groups == paste0(Target_gene,"_mut_2"))){
          return(geom_signif(na.rm = T, comparisons = list(
            c("Control", paste0(Target_gene,"_mut_1")),
            c("Control", paste0(Target_gene,"_mut_2"))),
            annotations = c(
              customHSD("Control",  paste0(Target_gene,"_mut_1"), data, type),
              customHSD("Control",  paste0(Target_gene,"_mut_2"), data, type)),
            map_signif_level = TRUE,
            step_increase = 0.2))
        } else if(any(all_groups == paste0(Target_gene,"_mut_3"))){
          return(geom_signif(na.rm = T, comparisons = list(
            c("Control", paste0(Target_gene,"_mut_1")),
            c("Control", paste0(Target_gene,"_mut_3"))),
            annotations = c(
              customHSD("Control",  paste0(Target_gene,"_mut_1"), data, type),
              customHSD("Control",  paste0(Target_gene,"_mut_3"), data, type)),
            map_signif_level = TRUE,
            step_increase = 0.2))
        }
      } else if(any(all_groups == paste0(Target_gene,"_mut_2"))){
        return(geom_signif(na.rm = T, comparisons = list(
          c("Control", paste0(Target_gene,"_mut_2")),
          c("Control", paste0(Target_gene,"_mut_3"))),
          annotations = c(
            customHSD("Control",  paste0(Target_gene,"_mut_2"), data, type),
            customHSD("Control",  paste0(Target_gene,"_mut_3"), data, type)),
          map_signif_level = TRUE,
          step_increase = 0.2))
      }
    } else if(length(all_groups) == 2){
      if(any(all_groups == paste0(Target_gene,"_mut_1"))){
        return(geom_signif(na.rm = T, comparisons = list(
          c("Control", paste0(Target_gene,"_mut_1"))),
          annotations = c(
            customHSD("Control",  paste0(Target_gene,"_mut_1"), data, type)),
          map_signif_level = TRUE,
          step_increase = 0.2))
      } else if(any(all_groups == paste0(Target_gene,"_mut_2"))){
        return(geom_signif(na.rm = T, comparisons = list(
          c("Control", paste0(Target_gene,"_mut_2"))),
          annotations = c(
            customHSD("Control",  paste0(Target_gene,"_mut_2"), data, type)),
          map_signif_level = TRUE,
          step_increase = 0.2))
      } else if(any(all_groups == paste0(Target_gene,"_mut_3"))){
        return(geom_signif(na.rm = T, comparisons = list(
          c("Control", paste0(Target_gene,"_mut_3"))),
          annotations = c(
            customHSD("Control",  paste0(Target_gene,"_mut_3"), data, type)),
          map_signif_level = TRUE,
          step_increase = 0.2))
      }
    }
  }
}

volcPlot <- function(df, i){

  TopGenes_neg <- df %>%
    filter(Screen == i, log2FC > 2 , Pval < 0.05 , Mutant_median > 0.2) %>%
    arrange(-Mutant_median)

  TopGenes_pos <- df %>%
    filter(Screen == i, log2FC < -2 , Pval < 0.05 , Control_median > 0.2) %>%
    arrange(-Control_median)

  if(nrow(TopGenes_neg) > 10){
    TopGenes_neg <- TopGenes_neg %>% slice(1:10)
  }
  if(nrow(TopGenes_pos) > 10){
    TopGenes_pos <- TopGenes_pos %>% slice(1:10)
  }

  TopGenes <- bind_rows(TopGenes_neg, TopGenes_pos)


  df %>%
    filter(Screen == i) %>%
    mutate(colours = case_when(
        (Pval < 0.05) & (log2FC > 2) ~ "SL",
        (Pval < 0.05) & (log2FC < -2) ~ "Adv",
        TRUE ~ "others"),
      colours = fct_relevel(colours, "others", "SL")) %>%
    ggplot(., aes(x = log2FC, y = -log10(Pval), colour = colours))+
    geom_point(size = 1)+
    labs(title = paste(i,"screen"),
         x = "Log2 fold change",
         y = "-log10(p-value)")+
    theme(legend.position = "none",
          text = element_text(size = 15))+
    scale_colour_manual(values=c("grey","dark red", "dark blue")) +
    geom_label_repel(aes(label = ifelse(GeneNames %in% TopGenes$GeneNames, GeneNames,"")),
                                box.padding   = 0.35,
                                point.padding = 0.5,
                                segment.color = 'grey50',
                                color = "black")

}

compare_dep <- function(gene){
  geneID <- dep_annot %>% filter(GeneNames %in% gene) %>% pull(GeneNameID)
  temp_dep <- dep %>% select(DepMap_ID, all_of(geneID)) %>%
    left_join(., Groups %>% select(DepMap_ID, Group)) %>%
    rename(DepProb = geneID)
  return(temp_dep)
}

get_protein_expr <- function(sample_id, hugo_symbol, uniprot_id = NULL){
  # sample_id <- Target_lines %>% filter(Group == "Mutant") %>% pull(DepMap_ID)
  # sample_id <- "xxxx"
  # hugo_symbol <- "KMT2D"

  if(!any(protein_nodup$Gene_Symbol %in% hugo_symbol)){
    return(NA)
  }

  extract_target_protein <- NULL
  extract_target_protein <- protein_nodup %>%
    filter(Gene_Symbol %in% hugo_symbol) %>%
    select(-c(Protein_Id,Description, Group_ID), ends_with("_TenPx**")) %>%
    pivot_longer(., -c(Gene_Symbol,Uniprot, Uniprot_Acc)) %>%
    mutate(DepMap_ID = map_chr(name, GetDepMapID))

  if(length(names(table(extract_target_protein$Uniprot_Acc))) == 1) {
    target_uniprot <- names(table(extract_target_protein$Uniprot_Acc))


    extract_target_protein <- extract_target_protein %>%
      filter(Uniprot_Acc == target_uniprot,
        str_detect(DepMap_ID, "^ACH-"),
        DepMap_ID %in% sample_id) %>%
      pull(value)

    if(length(extract_target_protein) == 0){
      return(NA)
    } else {
      return(extract_target_protein)
    }
  } else {
    print(paste0("Contains > 1 uniprot acc. id, please pick on of the following available: ", unique(extract_target_protein$Uniprot_Acc)))
  }
}
