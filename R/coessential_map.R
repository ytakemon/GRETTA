#' @title Perform co-essentially mapping
#' 
#' @description Performs multiple correlation coefficient analyses and determines cut to identify most likely co-essential genes.
#' 
#' @param input_gene string, A vector containing one Hugo Symbol, Default: NULL
#' @param input_disease string, A vector one or more disease contexts, Will perform pan-cancer analyses 
#' (all cell lines) by default, Default: NULL
#' @param input_cell_lines string, A vector DepMap_IDs for which co-essentiality mapping will be performed on. 
#' Will perform pan-cancer analyses (all cell lines) by default, Default: NULL
#' @param core_num integer, Number of cores to run analysis, Default: NULL
#' @param output_dir string, Full path to where output file should be saved, Default: NULL
#' @param data_dir string Path to GINIR_data
#' @param output_filename string name of file without the '.csv' extension. 
#' @param test logical, TRUE/FALSE whether you want to run only a small subset (first 10 genes) to ensure function will run properly 
#' prior to running all 18,333 genes. Default: FALSE.
#'
#' @return A data frame containing Pearson correlation coefficients. A copy is also saved to the 
#' directory defined in `output_dir`.
#' 
#' @details Description of output data frame
#' * `GeneNameID_A` - Hugo symbol with NCBI gene ID of query gene.
#' * `GeneNameID_B` - Hugo symbol with NCBI gene ID of all genes targeted in the DepMap KO screen.
#' * `estimate` - Correlation coefficient output from `?cor.test`.
#' * `statistic` - Pearson's correlation statistic. Output from `?cor.test`.
#' * `p.value` - P-value from Pearson's correlation statistic. Output from `?cor.test`.
#' * `parameter` - Degrees of freedom. Output from `?cor.test`.
#' * `Rank` - Rank by correlation coefficient. 
#' * `Padj_BH` - Benjamini-Hochberg adjusted p-value.
#' @md
#' 
#' @examples 
#' \dontrun{
#' 
#' Screen_results <- GINI_screen(
#' input_gene = 'ARID1A', 
#' output_dir = '~/Desktop/GINI_test_dir/',
#' data_dir = '/path/to/DepMap_data/',
#' test = TRUE) # turn on for shorter test runs
#' 
#' }
#' 
#' @rdname coessential_map
#' @export 
#' @importFrom parallel detectCores
#' @importFrom doMC registerDoMC
#' @importFrom dplyr mutate filter pull rename arrange case_when
#' @importFrom foreach `%dopar%` foreach
#' @importFrom rcompanion cliffDelta
#' @importFrom diptest dip.test
#' @importFrom broom tidy
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect everything
#' @importFrom readr write_csv
#' @importFrom stats cor.test
#' @importFrom tibble as_tibble

coessential_map <- function(
    input_gene = NULL, input_disease = NULL, input_cell_lines = NULL, core_num = NULL,
    output_dir = NULL, data_dir = NULL, output_filename = NULL, test = FALSE
) {
  
  # Check that essential inputs are given:
  if (is.null(input_gene)) {
    stop("No control IDs detected")
  }
  if (is.null(output_dir)) {
    output_dir <- paste0(getwd(), "/GINIR_", Sys.Date())
    print(paste0("No output directory specified. Creating: ", output_dir))
    dir.create(output_dir)
  }
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist. Please provide full path to directory.")
  }
  if (is.null(data_dir)) {
    stop(
      paste0("No directory to data was specified. Please provide path to DepMap data.")
    )
  }
  if (!dir.exists(data_dir)) {
    stop(
      paste0(
        "DepMap data directory does not exists. Please check again and provide the full path to the DepMap data directory."
      )
    )
  }
  if (!is.null(output_filename)) {
    output_dir_and_filename <- paste0(output_dir, "/", output_filename, ".csv")
  } else {
    output_dir_and_filename <- paste0(output_dir, "/GINI_coessentiality_network_results.csv")
  }
  
  if(is.null(c(input_disease, input_cell_lines))){
    
    # Do this for a default pan-cancer map
    # Load necessary data
    cat(
      "Performing default pan-cancer essentiality mapping for:", input_gene,"\n",
      "For this analysis, core_num is ignored."
    )
    
    fit <- NULL  # see: https://support.bioconductor.org/p/24756/
    load(
      paste0(data_dir, "pancan_coess_precomputed.rda"),
      envir = environment()
    )
    
    fit_est_long <- fit$r %>%
      tibble::as_tibble(.data, rownames = "GeneNameID_A") %>%
      tidyr::pivot_longer(-.data$GeneNameID_A, names_to = "GeneNameID_B", values_to = "estimate") %>%
      dplyr::filter(str_detect(.data$GeneNameID_A, paste0(input_gene,"_")))
    
    fit_tstat_long <- fit$t %>%
      tibble::as_tibble(., rownames = "GeneNameID_A") %>%
      tidyr::pivot_longer(-.data$GeneNameID_A, names_to = "GeneNameID_B", values_to = "statistic") %>%
      dplyr::filter(str_detect(.data$GeneNameID_A, paste0(input_gene,"_")))
    
    fit_pval_long <- fit$p %>%
      tibble::as_tibble(., rownames = "GeneNameID_A") %>%
      tidyr::pivot_longer(-.data$GeneNameID_A, names_to = "GeneNameID_B", values_to = "p.value") %>%
      dplyr::filter(str_detect(.data$GeneNameID_A, paste0(input_gene,"_")))
    
    fit_param_long <- fit$n %>%
      tibble::as_tibble(., rownames = "GeneNameID_A") %>%
      tidyr::pivot_longer(-.data$GeneNameID_A, names_to = "GeneNameID_B", values_to = "parameter") %>%
      dplyr::filter(str_detect(.data$GeneNameID_A, paste0(input_gene,"_")))
    
    cor_df <- dplyr::left_join(fit_est_long, fit_tstat_long) %>%
      dplyr::left_join(fit_pval_long) %>%
      dplyr::left_join(fit_param_long)
    
    # save and return output
    output <- cor_df %>%
      dplyr::arrange(-.data$estimate, .data$p.value) %>%
      dplyr::mutate(
        Rank = order(-.data$estimate, decreasing = F),
        Padj_BH = p.adjust(.data$p.value, method = "BH", n = (length(.data$p.value)))
      ) %>%
      readr::write_csv(file = output_dir_and_filename)
    
    print(
      paste0(
        "Coessentiality mapping finished. Outputs were also written to: ", output_dir_and_filename
      )
    )
    return(output) 
    
  } else {
    
    # Do this if not a default pan-cancer map
    # Detect cores if not defined
    if (is.null(core_num)) {
      cores_detected <- parallel::detectCores()
      print("No cores specified")
      print(paste0("Detected: ", cores_detected, " cores"))
      print(paste0("Using: ", cores_detected/2, " cores"))
      doMC::registerDoMC(cores_detected/2)
    }
    
    # Load necessary data
    gene_effect <- sample_annot <- NULL  # see: https://support.bioconductor.org/p/24756/
    load(
      paste0(data_dir, "/gene_effect.rda"),
      envir = environment()
    )
    load(
      paste0(data_dir, "/sample_annot.rda"),
      envir = environment()
    )
    
    # Set cores:
    if (!is.null(core_num)) {
      doMC::registerDoMC(core_num)
    }
    
    if (!is.null(input_disease)) {
      selected_cell_lines <- sample_annot %>%
        dplyr::filter(.data$DepMap_ID %in% gene_effect$DepMap_ID, .data$disease %in% input_disease) %>%
        dplyr::pull(.data$DepMap_ID)
    } else if (!is.null(input_cell_lines)) {
      selected_cell_lines <- sample_annot %>%
        dplyr::filter(
          .data$DepMap_ID %in% gene_effect$DepMap_ID, .data$DepMap_ID %in%
            input_cell_lines
        ) %>%
        dplyr::pull(.data$DepMap_ID)
    } else {
      stop(paste("Following may not be available:",input_disease, "\n or \n", input_cell_lines))
    }
    
    # account for differences between DepMap versions
    if (ncol(gene_effect) == 3) {
      AllGenes <- unique(gene_effect$GeneNameID)
      gene_effect_long <- gene_effect %>%
        dplyr::filter(.data$DepMap_ID %in% selected_cell_lines)
      
    } else if (ncol(gene_effect) > 3) {
      AllGenes <- colnames(gene_effect)[-1]  # removes DepMap_ID column  
      gene_effect_long <- gene_effect %>%
        tidyr::pivot_longer(
          cols = matches("\\d"),
          names_to = "GeneNameID", values_to = "Effect_score"
        ) %>%
        dplyr::filter(.data$DepMap_ID %in% selected_cell_lines)
      
    }
    
    Gene_A_GeneNameID <- get_GeneNameID(input_gene, data_dir = data_dir)
    Gene_A_effect <- gene_effect_long %>%
      dplyr::filter(.data$GeneNameID == Gene_A_GeneNameID)
    
    # Need to define function. A fix for a strange bug:
    `%dopar%` <- foreach::`%dopar%`
    
    # Begin loop
    print(
      "This may take a few mins... Consider running with a higher core numbers to speed up the analysis."
    )
    if (test == TRUE) {
      run <- 10
    } else {
      run <- length(unique(AllGenes))
    }
    res <- each <- NULL
    res <- foreach::foreach(each = 1:run, .combine = bind_rows) %dopar%
      {
        if (each == 1) {
          print(
            paste0(
              "Processing ", each, " of ", length(AllGenes),
              "\n"
            )
          )
        } else if (each == length(AllGenes)) {
          print(
            paste0(
              "Processing ", each, " of ", length(AllGenes),
              "\n"
            )
          )
        } else if (each%%1000 == 0) {
          print(
            paste0(
              "Processing ", each, " of ", length(AllGenes),
              "\n"
            )
          )
        }
        
        Gene_B_effect <- gene_effect_long %>%
          dplyr::filter(.data$GeneNameID == AllGenes[each])
        
        res_pearson <- cor.test(
          Gene_A_effect$Effect_score, Gene_B_effect$Effect_score, alternative = "two.sided",
          method = "pearson", na.action = "na.omit") %>%
          broom::tidy() %>%
          dplyr::mutate(GeneNameID_A = Gene_A_GeneNameID, GeneNameID_B = AllGenes[each]) %>%
          dplyr::select(.data$GeneNameID_A, .data$GeneNameID_B, tidyr::everything())
        
        res_pearson
      }
    
    # save and return output
    output <- res %>%
      dplyr::arrange(-.data$estimate, .data$p.value) %>%
      dplyr::mutate(
        Rank = order(-.data$estimate, decreasing = F),
        Padj_BH = p.adjust(.data$p.value, method = "BH", n = (length(.data$p.value)))) %>%
      dplyr::select(.data$GeneNameID_A, .data$GeneNameID_B, .data$estimate, 
                    .data$statistic, .data$parameter, .data$Rank, .data$Padj_BH) %>%
      readr::write_csv(file = output_dir_and_filename)
    
    print(
      paste0(
        "Coessentiality mapping finished. Outputs were also written to: ", output_dir_and_filename
      )
    )
    return(output) 
  }
}

