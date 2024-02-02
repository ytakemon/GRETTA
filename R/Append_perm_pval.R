#' @title Appends permutation adjusted p-value
#' 
#' @description Randomly samples dependency probabilities from a given mutant and control group screened.
#' 
#' @param p_value dataframe, Output from GI_screen(), Default: NULL
#' @param perm_pvalues dataframe, Output from GI_screen_perms(), Default: NULL
#' 
#' @return A vector of adjusted p-values
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
      Pval_perm_adj = purrr::map_dbl(Pval, .f = Pval_perm_adjust, perms$Pval))
  
  return(res)
}