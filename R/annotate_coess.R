#' @title Combine and annotate co-essential dataframe
#' 
#' @description Combines and annotates co-essential data frame with inflection points to determine which genes are 
#' likely to be candidates co-essential genes. 
#' 
#' @param input_ess data frame, A data frame output from `coessential_map()`, Default: NULL
#' @param input_inflec data frame, A data frame output from get_inflection_points(), Default: NULL
#' @param top_n vector, The number of top/bottom candidates to select (eg. if 10 there will be top top 
#' co-essential and 10 anti-essential genes), Default: NULL
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
#' annotated_df <- annotate_coess(coess_df, coess_inflection_df)
#' 
#' @rdname annotate_coess
#' @export 
#' @importFrom dplyr mutate filter pull rename arrange case_when
#' @importFrom tibble tibble

annotate_coess <- function(input_ess = NULL, input_inflec = NULL, top_n = NULL) {
  # Checkpoint
  if (is.null(input_ess)) {
    stop("No coessential dataframe found!")
  }
  if (is.null(input_inflec)) {
    stop("No inflection points found!")
  }
  if (!is.data.frame(input_ess)) {
    stop("Input is not a dataframe, please check the input")
  }
  
  All_res <- NULL
  for(g in seq_len(length(unique(input_ess$GeneNameID_A)))){
    # g <- 1
    gene <- unique(input_ess$GeneNameID_A)[g]
    
    select_input_ess <- input_ess %>%
      dplyr::filter(.data$GeneNameID_A %in% gene)
    
    select_input_inflec <- input_inflec %>%
      dplyr::filter(.data$GeneNameID_A %in% gene)
    
    res <- select_input_ess %>%
      dplyr::arrange(-.data$Rank) %>%
      dplyr::mutate(Candidate_gene = dplyr::case_when(
        (.data$Padj_BH < 0.05) & (.data$Rank >= select_input_inflec$Inflection_point_pos_byRank) ~ TRUE,
        (.data$Padj_BH < 0.05) & (.data$Rank <= select_input_inflec$Inflection_point_neg_byRank) ~ TRUE,
        TRUE ~ FALSE
      ))
    
    All_res <- dplyr::bind_rows(All_res, res)
  }
  
  candidate <- All_res %>% dplyr::filter(.data$Candidate_gene)
  top <- top_n
  middle <- 18333 - (top_n *2)
  if(length(unique(input_ess$GeneNameID_A)) == nrow(candidate)){
    if (is.null(top_n)) {
      stop("No candidates... Please add a number for top candidates")
    } else {
      message("No candidates... Manually selecting top ", top_n, " candidates.")
      
    }
    All_res <- All_res %>% mutate(
      Candidate_gene = rep(c(rep(TRUE, top), rep(FALSE, middle), rep(TRUE, top)), length(unique(input_ess$GeneNameID_A)))
    )
  }
  
  return(All_res)
}