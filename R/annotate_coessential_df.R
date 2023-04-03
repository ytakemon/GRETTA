#' @title Combine and annotate co-essential dataframe
#' 
#' @description Combines and annotates co-essential data frame with inflection points to determine which genes are 
#' likely to be candidates co-essential genes. 
#' 
#' @param input_ess data frame, A data frame output from `coessential_map()`, Default: NULL
#' @param input_inflec data frame, A data frame output from get_inflection_points(), Default: NULL
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
#' \dontrun{
#' 
#' annotated_coessential_df <- annotate_coessential_df(
#' input_ess = co_ess_res,
#' input_inflec = inflection_points)
#' 
#' }
#' 
#' @rdname annotate_coessential_df
#' @export 
#' @importFrom dplyr mutate filter pull rename arrange case_when
#' @importFrom tibble tibble

annotate_coessential_df <- function(input_ess = NULL, input_inflec = NULL) {
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
  
  res <- input_ess %>%
    dplyr::arrange(-.data$Rank) %>%
    dplyr::mutate(
      Candidate_gene = dplyr::case_when(
        (.data$Padj_BH < 0.05) & (.data$Rank >= input_inflec$Inflection_point_pos_byRank) ~
          TRUE, 
        (.data$Padj_BH < 0.05) & (.data$Rank <= input_inflec$Inflection_point_neg_byRank) ~
          TRUE, 
        TRUE ~ FALSE
      )
    )
  
  return(res)
}