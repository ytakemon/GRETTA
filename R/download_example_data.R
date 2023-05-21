#' @title Download minimal example data
#' 
#' @description Downloads minimal data for examples shown in the tutorial. This can also be downloaded manually from:
#' https://www.bcgsc.ca/downloads/ytakemon/GRETTA/22Q2_sample
#' 
#' @param path character. Path for example data. Default: NULL.
#' @param test logical. Used for examples. Default: FALSE.
#'
#' @return A data frame containing rank at which the threshold should be drawn for positive and negative co-essential genes.
#' 
#' @details All files are stored in './GRETTA_example'
#' @md
#' 
#' @examples 
#' download_example_data(path = ".", test = TRUE)
#' 
#' @rdname download_example_data
#' @export 
#' @importFrom tools R_user_dir
#' @importFrom BiocFileCache BiocFileCache bfcrpath

download_example_data <- function(path = NULL, test = FALSE) {
    
    if(test == FALSE){
      # URL to public example data
      url1 <- "https://www.bcgsc.ca/downloads/ytakemon/GRETTA/22Q2_sample/CCLE_exp_annot.rda"
      url2 <- "https://www.bcgsc.ca/downloads/ytakemon/GRETTA/22Q2_sample/CCLE_exp.rda"
      url3 <- "https://www.bcgsc.ca/downloads/ytakemon/GRETTA/22Q2_sample/copy_num_annot.rda"
      url4 <- "https://www.bcgsc.ca/downloads/ytakemon/GRETTA/22Q2_sample/copy_num.rda"
      url5 <- "https://www.bcgsc.ca/downloads/ytakemon/GRETTA/22Q2_sample/dep_annot.rda"
      url6 <- "https://www.bcgsc.ca/downloads/ytakemon/GRETTA/22Q2_sample/dep.rda"
      url7 <- "https://www.bcgsc.ca/downloads/ytakemon/GRETTA/22Q2_sample/mut_calls.rda"
      url8 <- "https://www.bcgsc.ca/downloads/ytakemon/GRETTA/22Q2_sample/protein_annot.rda"
      url9 <- "https://www.bcgsc.ca/downloads/ytakemon/GRETTA/22Q2_sample/protein_nodup.rda"
      url10 <- "https://www.bcgsc.ca/downloads/ytakemon/GRETTA/22Q2_sample/sample_22Q2_ARID1A_coessential_inflection.rda"
      url11 <- "https://www.bcgsc.ca/downloads/ytakemon/GRETTA/22Q2_sample/sample_22Q2_ARID1A_coessential_result.rda"
      url12 <- "https://www.bcgsc.ca/downloads/ytakemon/GRETTA/22Q2_sample/sample_22Q2_ARID1A_KO_screen.rda.rda"
      url13 <- "https://www.bcgsc.ca/downloads/ytakemon/GRETTA/22Q2_sample/sample_annot.rda"
      url14 <- "https://www.bcgsc.ca/downloads/ytakemon/GRETTA/22Q2_sample/gene_effect.rda"
      
      cache <- tools::R_user_dir("GRETTA", which = "cache")
      bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
      pathsToLoad <- BiocFileCache::bfcrpath(bfc, c(url1, url2, url3, url4, url5, url6, url7,
                                                    url8, url9, url10, url11, url12, url13, url14))
      
      mut_calls <- copy_num <- copy_num_annot <- dep <- dep_annot <- sample_annot <- 
        CCLE_exp <- CCLE_exp_annot <- protein_nodup <- protein_annot <- gene_effect <- screen_results <- coess_annotated_df <- coess_df <-
        coess_inflection_df <- NULL
      for (i in seq(length(pathsToLoad))) {
        load(pathsToLoad[i])
      }
      
      dir.create(paste0(path,"/GRETTA_example"))
      dir.create(paste0(path,"/GRETTA_example_output"))
      
      save(mut_calls, file = paste0(path,"/GRETTA_example/mut_calls.rda"))
      save(copy_num_annot, file = paste0(path,"/GRETTA_example/copy_num_annot.rda"))
      save(copy_num, file = paste0(path,"/GRETTA_example/copy_num.rda"))
      save(dep, file = paste0(path,"/GRETTA_example/dep.rda"))
      save(dep_annot, file = paste0(path,"/GRETTA_example/dep_annot.rda"))
      save(sample_annot, file = paste0(path,"/GRETTA_example/sample_annot.rda"))
      save(CCLE_exp, file = paste0(path,"/GRETTA_example/CCLE_exp.rda"))
      save(CCLE_exp_annot, file = paste0(path,"/GRETTA_example/CCLE_exp_annot.rda"))
      save(protein_nodup, file = paste0(path,"/GRETTA_example/protein_nodup.rda"))
      save(protein_annot, file = paste0(path,"/GRETTA_example/protein_annot.rda"))
      save(gene_effect, file = paste0(path,"/GRETTA_example/gene_effect.rda"))
      save(screen_results, file = paste0(path,"/GRETTA_example/sample_22Q2_ARID1A_KO_screen.rda"))
      save(coess_df, file = paste0(path,"/GRETTA_example/sample_22Q2_ARID1A_coessential_result.rda"))
      save(coess_inflection_df, file = paste0(path,"/GRETTA_example/sample_22Q2_ARID1A_coessential_inflection.rda"))
      
      message("Data saved to: ", path,"/GRETTA_example/" )
    } else {
      message("Data would be saved to: ./GRETTA_example/" )
    }
     
  }