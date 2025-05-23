#' @title Perform permutation test n genetic interaction screen
#' 
#' @description Randomly samples dependency probabilities from a given mutant and control group screened.
#' 
#' @param p_value numeric, A vector of p-value generated by GI_screen(), Default: NULL
#' @param perm_pvalues numeric, A vector of p-values generated by GI_screen_perms(), Default: NULL
#' 
#' @return A vector of adjusted p-values
#' 
#' @details P-values are adjusted based on permuted random sampled background data. 
#' @md
#' 
#' @examples 
#' # Please see tutorial
#' 
#' @rdname Pval_perm_adjust
#' @export 

Pval_perm_adjust <- function(p_value = NULL, perm_pvalues = NULL) {
  if(is.na(p_value)|p_value == Inf){
    return(NA_integer_)

  } else if(p_value == 1){
    return(1)
    
  } else {
    x <- table(perm_pvalues <= p_value)
    if(length(x) < 2){
      res <- 1/(x[[1]]+1)
    } else {
      res <- x[["TRUE"]]/(x[["FALSE"]]+ x[["TRUE"]]+1)
    }
    return(res)
  }
}