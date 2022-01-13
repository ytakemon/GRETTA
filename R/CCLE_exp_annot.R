#' @title DepMap 20Q1: Cancer cell line RNA-seq expression annotation
#'
#' @description A data set containing gene name column headers from the "CCLE_exp" data frame and
#' its various forms that exist in the DepMap data set platform.
#'
#' @format A data frame with 19145 rows and 4 variables:
#' \describe{
#'   \item{\code{names}}{character "CCLE_exp" column names}
#'   \item{\code{GeneNames}}{character Hugo symbols}
#'   \item{\code{GeneID}}{character NCBI gene IDs}
#'   \item{\code{GeneNameID}}{character Hugo symbol separated by NCBI gene ID} 
#'}
#' @source \url{https://figshare.com/articles/dataset/DepMap_20Q1_Public/11791698}
"CCLE_exp_annot"
