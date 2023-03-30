#' @title Calculate inflection points of co-essential data frame
#' 
#' @description Calculates the inflection point of the positive and negative curves to determine threshold for 
#' co-essential genes. See `?coessential_map` for details.
#' 
#' @param input_coessential_df data frame, A data frame output from `coessential_map()`, Default: NULL.
#'
#' @return A data frame containing rank at which the threshold should be drawn for positive and negative co-essential genes.
#' 
#' @details Description of output data frame
#' * `Inflection_point_pos_byRank` - Rank threshold for positive curve.
#' * `Inflection_point_neg_byRank` - Rank threshold for negative curve.
#' @md
#' 
#' @examples 
#' \dontrun{
#' 
#' inflection_points <- get_inflection_points(input_coessential_df = coess_df)
#' 
#' }
#' 
#' @rdname get_inflection_points
#' @export 
#' @importFrom dplyr mutate filter pull rename arrange case_when
#' @importFrom tibble tibble
#' @importFrom RootsExtremaInflections inflexi

get_inflection_points <- function(input_coessential_df = NULL) {
  # Checkpoint
  if (is.null(input_coessential_df)) {
    stop("No coessential dataframe found!")
  }
  if (!is.data.frame(input_coessential_df)) {
    stop("Input is not a dataframe, please check the input")
  }
  
  # Calculate inflection point
  print("This may take a few mins...")
  inflection_points <- NULL
  
  # Positive curve
  subset_df_pos <- input_coessential_df %>%
    dplyr::arrange(.data$estimate, .data$Padj_BH) %>%
    dplyr::filter(.data$estimate > 0)
  print("Calculating inflection point of positive curve.\n")
  x_pos <- subset_df_pos$Rank
  y_pos <- subset_df_pos$estimate
  fit_pos <- RootsExtremaInflections::inflexi(
    x_pos, y_pos, 1, length(x_pos),
    5, 5, plots = F, doparallel = F
  )
  fit_pos$an
  fit_pos$finfl
  inflection_point_pos <- fit_pos$finfl[2]
  
  # Negative curve
  subset_df_neg <- input_coessential_df %>%
    dplyr::arrange(.data$estimate, .data$Padj_BH) %>%
    dplyr::filter(.data$estimate < 0)
  print("Calculating inflection point of negative curve.\n")
  x_neg <- subset_df_neg$Rank
  y_neg <- subset_df_neg$estimate
  fit_neg <- RootsExtremaInflections::inflexi(
    x_neg, y_neg, 1, length(x_neg),
    5, 5, plots = F, doparallel = F
  )
  fit_neg$an
  fit_neg$finfl
  inflection_point_neg <- fit_neg$finfl[2]
  
  res <- tibble::tibble(
    Inflection_point_pos_byRank = inflection_point_pos, Inflection_point_neg_byRank = inflection_point_neg
  )
}
