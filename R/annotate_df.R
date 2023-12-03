#' @title Combine and annotate co-essential or co-expression dataframe
#' 
#' @description Combines and annotates co-essential data frame with inflection points to determine which genes are 
#' likely to be candidates co-essential genes. 
#' 
#' @param input_df data frame, A data frame output from `coessential_map()`, `rna_express()`, and `protein_express()`, Default: NULL
#' @param input_inflec data frame, A data frame output from get_inflection_points(), Default: NULL
#' @param top_n vector, The number of top/bottom candidates to select (eg. if 10 there will be top top 
#' co-essential and 10 anti-essential genes), Default: NULL
#' @param use_inflection logical, TRUE to use inflection points as threshold in input_inflec or 
#' FALSE to use top_n as cutoff , Default: TRUE
#' 
#' @return The input data frame containing co-essential correlation coefficients will be annotated with an additional column 
#' `Candidate_gene` to indicate whether the gene is considered to a possible co-essential gene.
#' 
#' @details Description of output data frame
#' * `Candidate_gene` - Logical; Whether the coefficient had a p-value < 0.05 and was at or above (positive curve) 
#' or below (negative curve) of the inflection point. 
#' * See `?coessential_map` for details.
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
#' load(paste0(
#' gretta_data_dir,'/sample_22Q2_ARID1A_coessential_result.rda'), 
#' envir = environment())
#' load(paste0(
#' gretta_data_dir,'/sample_22Q2_ARID1A_coessential_inflection.rda'), 
#' envir = environment())
#' 
#' annotated_df <- annotate_df(coess_df, coess_inflection_df)
#' 
#' @rdname annotate_df
#' @export 
#' @importFrom dplyr mutate filter pull rename arrange case_when
#' @importFrom tibble tibble

annotate_df <- function(input_df = NULL, input_inflec = NULL, top_n = NULL, use_inflection = TRUE) {
  # Checkpoint
  if (is.null(input_df)) {
    stop("No coessential dataframe found!")
  }
  if (is.null(input_inflec)) {
    stop("No inflection points found!")
  }
  if (!is.data.frame(input_df)) {
    stop("Input is not a dataframe, please check the input")
  }
  
  All_res <- NULL
  for(g in seq_len(length(unique(input_df$GeneNameID_A)))){
    # g <- 1
    gene <- unique(input_df$GeneNameID_A)[g]
    
    select_input_df <- input_df %>%
      dplyr::filter(.data$GeneNameID_A %in% gene)
    
    if(nrow(input_inflec) == 1){
      select_input_inflec <- input_inflec 
    } else {
      select_input_inflec <- input_inflec %>%
        dplyr::filter(.data$GeneNameID_A %in% gene)
    }
    
    res <- select_input_df %>%
      dplyr::arrange(-.data$Rank) %>%
      dplyr::mutate(Candidate_gene = dplyr::case_when(
        (.data$Padj_BH < 0.05) & (.data$Rank >= select_input_inflec$Inflection_point_pos_byRank) ~ TRUE,
        (.data$Padj_BH < 0.05) & (.data$Rank <= select_input_inflec$Inflection_point_neg_byRank) ~ TRUE,
        TRUE ~ FALSE
      ))
    
    All_res <- dplyr::bind_rows(All_res, res)
  }
  
  # Check if manual threshold is needed
  candidate <- All_res %>% dplyr::filter(.data$Candidate_gene)
  if(length(unique(input_df$GeneNameID_A)) == nrow(candidate)){
    use_inflection == FALSE
    top <- top_n
    middle <- 18333 - (top_n * 2)
    message("Only ", nrow(candidate), " candidates. Selecting candidates based on input \n")
  } 
  
  # If using a manual threshold
  if(use_inflection == FALSE){
    if (is.null(top_n)) {
      stop("Manually  selecting candidates... Error: No input in top_n \n")
    } else {
      message("Selecting top and bottom ", top_n, " candidates.")
    }
    All_res <- All_res %>% dplyr::ungroup() %>% dplyr::mutate(
      Candidate_gene = rep(c(rep(TRUE, top), rep(FALSE, middle), rep(TRUE, top)), length(unique(input_df$GeneNameID_A)))
    )
  } else {
    message("Selecting candidates based on inflection points.")
  }
  
  return(All_res)
}