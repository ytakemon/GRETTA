#' @title Perform genetic interaction screen
#' 
#' @description Compares dependency probabilities of mutant and control groups to determine 
#' whether the mutant group are resistant or sensitive to specific gene perturbations.
#' 
#' @param control_IDs string, A vector containing two or more DepMap_IDs, Default: NULL
#' @param mutant_IDs string, A vector containing two or more DepMap_IDs, Default: NULL
#' @param core_num integer, Number of cores to run analysis, Default: NULL
#' @param output_dir string, Full path to where output file should be saved, Default: NULL
#' @param data_dir string Path to GINIR_data
#' @param test logical, For test_that to shorten computational time for testing
#'
#' @return A data frame containing results from the genetic screen. A copy is also saved to the 
#' directory defined in `output_dir`.
#' 
#' @details Description of output data frame
#' * `GeneName_ID` - Hugo symbol with NCBI gene ID
#' * `GeneNames` - Hugo symbol 
#' * `_median`, `_mean`, `_sd`, `_iqr` - Control and mutant group's median, mean, standard deviation (sd), 
#' and interquartile range (iqr) of dependency probabilities. Dependency probabilities range from zero to one, 
#' where one indicates a essential gene (ie. KO of gene was lethal) and zero indicates a non-essential gene 
#' (KO of gene was not lethal)
#' * `Pval` - P-value from Mann Whitney U test between control and mutant groups.
#' * `Adj_pval` - BH-adjusted P-value.
#' * `log2FC_by_median` - Log2 normalized median fold change of dependency probabilities (mutant / control).
#' * `log2FC_by_mean` - Log2 normalized mean fold change of dependency probabilities (mutant / control).
#' * `CliffDelta` - Cliff's delta non-parametric effect size between mutant and control dependency probabilities. 
#' Ranges between -1 to 1.
#' * `dip_pval` - Hartigan's dip test p-value. Tests whether distribution of mutant dependency probability is unimodel.
#' If dip test is rejected (p-value < 0.05), this indicates that there is a multimodel dependency probability distribution and
#' that there may be another factor contributing to this separation. 
#' * `Interaction_score` - Combined value generated from signed p-values: `-log10(Pval) \* sign(log2FC_by_median)`
#' @md
#' 
#' @examples 
#' \dontrun{
#' 
#' Screen_results <- GINI_screen(
#' control_IDs = c("ACH-001354", "ACH-000274", "ACH-001799"), 
#' mutant_IDs = c("ACH-000911", "ACH-001957", "ACH-000075"), 
#' core_num = 2, 
#' output_dir = "~/Desktop/GINI_test_dir/",
#' data_dir = "/path/to/DepMap_data/",
#' test = TRUE) # turn on for shorter test runs
#' 
#' }
#' 
#' @rdname GINI_screen
#' @export 
#' @importFrom parallel detectCores
#' @importFrom doMC registerDoMC
#' @importFrom dplyr mutate filter group_by summarize pull rename
#' @importFrom forcats fct_relevel
#' @importFrom foreach `%dopar%` foreach
#' @importFrom rcompanion cliffDelta
#' @importFrom diptest dip.test
#' @importFrom broom tidy
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tidyselect everything
#' @importFrom readr write_csv
#' @importFrom stats median sd IQR wilcox.test p.adjust

GINI_screen <- function(control_IDs = NULL, mutant_IDs = NULL, core_num = NULL, output_dir = NULL, data_dir = NULL, test = FALSE){
  
  # Check that essential inputs are given:
  if(is.null(control_IDs)){
    stop("No control IDs detected")
  }
  if(is.null(mutant_IDs)){
    stop("No mutant IDs detected")
  }
  if(is.null(core_num)){
    cores_detected <-  parallel::detectCores()
    print("No cores specified")
    print(paste0("Detected: ", cores_detected," cores"))
    print(paste0("Using: ", cores_detected/2," cores"))
    doMC::registerDoMC(cores_detected/2)
  }
  if(is.null(output_dir)){
    output_dir <- paste0(getwd(),"/GINIR_",Sys.Date())
    print(paste0("No output directory specified. Creating: ", output_dir))
    dir.create(output_dir)
  }
  if(!dir.exists(output_dir)){
    stop("Output directory does not exist. Please provide full path to directory.")
  }
  if(is.null(data_dir)){
    stop(paste0("No directory to data was specified. Please provide path to DepMap data."))
  }
  if(!dir.exists(data_dir)){
    stop(paste0("DepMap data directory does not exists. Please check again and provide the full path to the DepMap data directory."))
  }
  
  # Set cores:
  if(!is.null(core_num)){
    doMC::registerDoMC(core_num)
  }
  
  # Check to see enough samples were given:
  if(length(control_IDs) < 2){
    stop("Not enough controls! Provide at least two.")
  }
  if(length(mutant_IDs) < 2){
    stop("Not enough mutants! Provide at least two.")
  }
  
  # Load necessary data
  dep <- dep_annot <- NULL # see: https://support.bioconductor.org/p/24756/
  load(paste0(data_dir, "/dep.rda"), envir = environment())
  load(paste0(data_dir, "/dep_annot.rda"), envir = environment())
  
  # Check to see if enough samples were given after filtering:
  Control_group_avail <- control_IDs[control_IDs %in% dep$DepMap_ID]
  Mutant_groups_avail <- mutant_IDs[mutant_IDs %in% dep$DepMap_ID]
  if(length(Control_group_avail) < 2){
    stop(paste0("Not enough controls were screened! Only the following control samples were screen: ", 
                paste0(Control_group_avail, collapse = ", ")))
  }
  if(length(Control_group_avail) < 2){
    stop(paste0("Not enough mutants were screened! Only the following mutant samples were screen: ", 
                paste0(Mutant_groups_avail, collapse = ", ")))
  }
  
  # Filter dep probs to only those that are used:
  select_dep <- dep %>%
    tidyr::pivot_longer(cols = tidyselect::matches("\\d"), names_to = "GeneNameID", values_to = "DepProb") %>%
    dplyr::mutate(CellType = case_when(
      DepMap_ID %in% Mutant_groups_avail ~ "Mutant",
      DepMap_ID %in% Control_group_avail ~ "Control",
      TRUE ~ "Others")) %>%
    dplyr::filter(.data$CellType != "Others") %>%  
    dplyr::mutate(CellType = forcats::fct_relevel(.data$CellType, "Control", "Mutant"))
  
  
  # Need to define function. A fix for a strange bug:
  `%dopar%` <- foreach::`%dopar%`
  
  # For testing short loops:
  if(test == TRUE){
    run <- 10
  } else if (test == FALSE){
    run <- length(unique(select_dep$GeneNameID))
  }
  
  # Begin loop
  All_res <- each <- NULL
  All_res <- foreach::foreach(each = 1:run, .combine = bind_rows) %dopar% {
    
    # Give feedback
    if(each == 1){
      print(paste0("Processing ", each, " of ", length(unique(dep$GeneNameID))))
    } else if(each == length(unique(dep$GeneNameID))){
      print(paste0("Processing ", each, " of ", length(unique(dep$GeneNameID))))
    } else if(each%%1000 == 0){
      print(paste0("Processing ", each, " of ", length(unique(dep$GeneNameID))))
    }
    
    # Get each gene
    geneID <- unique(select_dep$GeneNameID)[each]
    df <- select_dep %>% dplyr::filter(.data$GeneNameID == geneID)
    
    if(all(df$DepProb == 0)){
      populate <- rep(0,11)
      
    } else {
      # # MWU doesn't handle na or zero's well so
      # # FOR NOW remove zeros.
      # df <- df %>% filter(!is.na(DepProb)) %>%
      #   filter(DepProb != 0)
      
      stats <- df %>%
        dplyr::group_by(.data$CellType) %>%
        dplyr::summarize(Median = stats::median(.data$DepProb, na.rm = TRUE),
                         Mean = mean(.data$DepProb, na.rm = TRUE),
                         SD = stats::sd(.data$DepProb, na.rm = TRUE),
                         IQR = stats::IQR(.data$DepProb, na.rm = TRUE),
                         .groups = "drop")
      
      if((any(is.na(stats)) != TRUE) & (nrow(stats) == 2)){
        
        fit_pval <- stats::wilcox.test(DepProb ~ CellType, df,
                                paired = F,
                                alternative = "two.sided",
                                conf.int = T,
                                na.action = "na.omit")$p.value
        
        # If group size is < 3 cliffDelta will have error:
        # missing value where TRUE/FALSE needed
        # Important note: because celltype has a specific factor order specified for select_dep, a delta > 0 indicates an effect score greater in the Control (first level) and a delta < 0 means an effect score greater in the Mutant( second level)
        
        CliffDelta <- rcompanion::cliffDelta(DepProb ~ CellType, df)
        
        # Add diptest for uni/multi modality
        # null hypothesis if p > 0.05 the data is unimodal
        # alternative hyp if p < 0.05 the data is multimodal
        dip_pval <- df %>%
          dplyr::filter(.data$CellType == "Mutant") %>%
          dplyr::pull(.data$DepProb) %>%
          diptest::dip.test() %>% 
          broom::tidy() %>%
          dplyr::pull(.data$p.value)
        
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
      dplyr::mutate(!!sym(geneID) := populate) %>%
      tidyr::pivot_longer(-.data$Result) %>%
      tidyr::pivot_wider(names_from = .data$Result, values_from = .data$value) %>%
      dplyr::rename(GeneNameID = .data$name)
  } # End of for loop
  
  # Add mutant group name
  output <- All_res %>%
    dplyr::mutate(log2FC_by_median = log2(.data$Mutant_median / .data$Control_median),
           log2FC_by_mean = log2(.data$Mutant_mean / .data$Control_mean),
           Adj_pval = p.adjust(.data$Pval, method = "BH", length(.data$Pval)),
           Interaction_score = -log10(.data$Pval) * sign(.data$log2FC_by_median)) %>%
    dplyr::left_join(dep_annot %>% select(.data$GeneNameID, .data$GeneNames),
              by = "GeneNameID") %>%
    dplyr::select(.data$GeneNameID, .data$GeneNames, .data$Control_median:.data$Pval, .data$Adj_pval,
           .data$log2FC_by_median, .data$log2FC_by_mean, tidyselect::everything(), .data$Interaction_score)
  
  # save and return output
  output %>% readr::write_csv(file = paste0(output_dir,"/GINI_screening_results.csv"))
  return(output)
}
