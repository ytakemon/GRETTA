#' @title Perform permutation test n genetic interaction screen
#' 
#' @description Randomly samples dependency probabilities from a given mutant and control group screened.
#' 
#' @param control_id string, A vector containing two or more DepMap_id, Default: NULL
#' @param mutant_id string, A vector containing two or more DepMap_id, Default: NULL
#' @param n_perm integer, Number of permutations to run, Default: 100
#' @param core_num integer, Number of cores to run analysis, Default: NULL
#' @param output_dir string, Full path to where output file should be saved, Default: NULL
#' @param data_dir string Path to GRETTA_data
#' @param filename string name of file without the '.csv' extension. 
#'
#' @return A data frame containing results from the genetic screen. A copy is also saved to the 
#' directory defined in `output_dir`.
#' 
#' @details Description of output data frame
#' * `perm_test` - permutation iteration
#' * `_median`, `_mean` - Control and mutant group's median, mean, of dependency probabilities. Dependency probabilities range from zero to one, 
#' where one indicates a essential gene (ie. KO of gene was lethal) and zero indicates a non-essential gene 
#' (KO of gene was not lethal)
#' * `Pval` - P-value from Mann Whitney U test between control and mutant groups.
#' * `n_perm` - Number of permutations assigned
#' @md
#' 
#' @examples 
#' gretta_data_dir <- './GRETTA_example/'
#' gretta_output_dir <- './GRETTA_example_output/'
#' 
#' if(!dir.exists(gretta_data_dir)){
#'   download_example_data(".")
#' }
#' 
#' if(!dir.exists(gretta_data_dir)){
#'   download_example_data(".")
#' }
#' 
#' Screen_results <- GI_screen_perms(
#' control_id = c('ACH-001354', 'ACH-000274', 'ACH-001799'), 
#' mutant_id = c('ACH-000911', 'ACH-001957', 'ACH-000075'), 
#' n_perm = 10,
#' core_num = 2, 
#' output_dir = gretta_output_dir,
#' data_dir = gretta_data_dir)
#' 
#' @rdname GI_screen_perms
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

GI_screen_perms <- function(control_id = NULL, mutant_id = NULL, n_perm = 100,
                      core_num = NULL, output_dir = NULL,
                      data_dir = NULL, filename = NULL) {
  
  # Check that essential inputs are given:
  if (is.null(control_id)) {
    stop("No control IDs detected")
  }
  if (is.null(mutant_id)) {
    stop("No mutant IDs detected")
  }
  if (is.null(core_num)) {
    cores_detected <- parallel::detectCores()
    message("No cores specified")
    message("Detected: ", cores_detected, " cores")
    message("Using: ", cores_detected/2, " cores")
    doMC::registerDoMC(cores_detected/2)
  }
  if (is.null(output_dir)) {
    output_dir <- paste0(getwd(), "/GRETTA_", Sys.Date())
    say <- paste0("No output directory specified. Creating: ",
                  output_dir)
    message(say)
    dir.create(output_dir)
  }
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist. Please provide full path to directory.")
  }
  if (is.null(data_dir)) {
    stop("No directory to data was specified. Please provide path to DepMap data.")
  }
  if (!dir.exists(data_dir)) {
    stop("DepMap data directory does not exists. Please check again and provide the full path to the DepMap data directory.")
  }
  if (!is.null(filename)) {
    output_dir_and_filename <- paste0(output_dir,
                                      "/", filename, ".csv")
  } else {
    output_dir_and_filename <- paste0(output_dir,
                                      "/GRETTA_GI_screen_perms.csv")
  }
  
  # Set cores:
  if (!is.null(core_num)) {
    doMC::registerDoMC(core_num)
  }
  
  # Check to see enough samples were given:
  if (length(control_id) < 2) {
    stop("Not enough controls! Provide at least two.")
  }
  if (length(mutant_id) < 2) {
    stop("Not enough mutants! Provide at least two.")
  }
  
  # Load necessary data
  dep <- dep_annot <- NULL  # see: https://support.bioconductor.org/p/24756/
  load(paste0(data_dir, "/dep.rda"), envir = environment())
  load(paste0(data_dir, "/dep_annot.rda"), envir = environment())
  
  # Check to see if enough samples were given
  # after filtering:
  Control_group_avail <- control_id[control_id %in%
                                      dep$DepMap_ID]
  Mutant_groups_avail <- mutant_id[mutant_id %in%
                                     dep$DepMap_ID]
  if (length(Control_group_avail) < 2) {
    say <- paste0("Not enough controls were screened! Only the following control samples were screen: ",
                  paste0(Control_group_avail, collapse = ", "))
    stop(say)
  }
  if (length(Control_group_avail) < 2) {
    say <- paste0("Not enough mutants were screened! Only the following mutant samples were screen: ",
                  paste0(Mutant_groups_avail, collapse = ", "))
    stop(say)
  }
  
  # Filter dep probs to only those that are
  # used:
  select_dep <- dep %>%
    tidyr::pivot_longer(cols = tidyselect::matches("\\d"),
                        names_to = "GeneNameID", values_to = "DepProb") %>%
    dplyr::mutate(CellType = case_when(.data$DepMap_ID %in%
                                         Mutant_groups_avail ~ "Mutant", .data$DepMap_ID %in%
                                         Control_group_avail ~ "Control", TRUE ~
                                         "Others")) %>%
    dplyr::filter(.data$CellType != "Others") %>%
    dplyr::mutate(CellType = forcats::fct_relevel(.data$CellType,
                                                  "Control", "Mutant"))
  
  # Need to define function. A fix for a
  # strange bug:
  `%dopar%` <- foreach::`%dopar%`
  
  # Begin loop
  All_res <- each <- NULL
  All_res <- foreach::foreach(each = seq_len(n_perm),
                              .combine = bind_rows) %dopar% {
                                
                                # Give feedback
                                if (each == 1) {
                                  message("Processing ", each, " of ", n_perm)
                                } else if (each == n_perm) {
                                  message("Processing ", each, " of ", n_perm)
                                } else if (each%%1000 == 0) {
                                  message("Processing ", each, " of ", n_perm)
                                }
                                
                                # Create randomly sampled DepProb dataframe
                                dummy_geneID <- "A1BG_1"
                                df <- select_dep %>% 
                                  mutate(DepProb_randomize = (sample(.data$DepProb, size = length(.data$DepProb), replace = FALSE))) %>%
                                  filter(.data$GeneNameID == dummy_geneID) %>%
                                  select(-.data$DepProb) %>%
                                  rename(DepProb = .data$DepProb_randomize)
                                
                                df_post_filter_check <- df %>%
                                  count(.data$CellType)
                                
                                if (any(df_post_filter_check$n < 2)) {
                                  populate <- rep(NA, 5)
                                  
                                } else if (all(df$DepProb == 0)) {
                                  populate <- rep(0, 5)
                                  
                                } else if (all(df$DepProb == 1)) {
                                  populate <- rep(1, 5)
                                  
                                } else {
                                  # # MWU doesn't handle na or zero's
                                  # well so # FOR NOW remove zeros.  df
                                  # <- df %>% filter(!is.na(DepProb))
                                  # %>% filter(DepProb != 0)
                                  
                                  stats <- df %>%
                                    dplyr::group_by(.data$CellType) %>%
                                    dplyr::summarize(
                                      Median = stats::median(.data$DepProb, na.rm = TRUE), 
                                      Mean = mean(.data$DepProb, na.rm = TRUE), .groups = "drop")
                                  
                                  if ((any(is.na(stats)) != TRUE) & (nrow(stats) == 2)) {
                                    fit_pval <- stats::wilcox.test(DepProb ~ CellType, df, paired = FALSE, alternative = "two.sided",
                                                                   conf.int = TRUE, na.action = "na.omit")$p.value
                                    
                                  } else if ((any(is.na(stats)) == TRUE) &
                                             (nrow(stats) == 2)) {
                                    populate <- rep(0, 5)
                                  }
                                  
                                  if ((any(is.na(stats)) != TRUE) & (nrow(stats) ==
                                                                     2)) {
                                    populate <- as.numeric(c(unlist(stats)[-c(1,
                                                                              2)], fit_pval))
                                  } else {
                                    populate <- rep(0, 5)
                                  }
                                }
                                
                                tibble(Result = c("Control_median", "Mutant_median",
                                                  "Control_mean", "Mutant_mean",
                                                  "Pval")) %>%
                                  dplyr::mutate(!!sym(dummy_geneID) := populate) %>%
                                  tidyr::pivot_longer(-.data$Result) %>%
                                  tidyr::pivot_wider(names_from = .data$Result,
                                                     values_from = .data$value) %>%
                                  dplyr::rename(GeneNameID = .data$name)
                              }  # End of for loop
  
  # Add mutant group name
  output <- All_res %>%
    mutate(
      perm_test = 1:n_perm,
      n_perm = n_perm) %>%
    select(-.data$GeneNameID) %>%
    select(.data$perm_test, everything())
  
  # save and return output
  output %>%
    readr::write_csv(file = output_dir_and_filename)
  
  message("Permutations finished. Outputs were also written to: ",
          output_dir_and_filename)
  return(output)
}