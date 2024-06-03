#' @title Appends permutation adjusted p-value
#' 
#' @description Appends permutation adjusted p-values to new column.
#' 
#' @param GI_screen_res dataframe, Output from GI_screen(), Default: NULL
#' @param perms dataframe, Output from Pval_perm_adjust(), Default: NULL
#' 
#' @return A dataframe with permutation adjusted p-values appeneded as a new column named Pval_perm_adj.
#' 
#' @details P-values are adjusted based on permuted random sampled background data. 
#' @md
#' 
#' @examples 
#' # Please see tutorial
#' 
#' @rdname Append_perm_pval
#' @export 
#' @importFrom dplyr mutate filter group_by summarize pull rename
#' @importFrom purrr map_dbl
Append_perm_pval <- function(GI_screen_res = NULL, perms = NULL) {
  res <- GI_screen_res %>%
    dplyr::mutate(
      Pval_perm_adj = purrr::map_dbl(.data$Pval, .f = Pval_perm_adjust, perms$Pval),
      Interaction_score = -log10(.data$Pval_perm_adj) * sign(.data$log2FC_by_median)
    )
  
  return(res)
}