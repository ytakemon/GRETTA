#' @title Perform RNA co-expression analysis
#' 
#' @description Performs multiple correlation coefficient analyses and determines cut to identify most likely co-expressed genes.
#' 
#' @param input_genes string, A vector containing one or more Hugo Symbol, Default: NULL
#' @param input_disease string, A vector one or more disease contexts, Will perform pan-cancer analyses 
#' (all cell lines) by default, Default: NULL
#' @param input_cell_lines string, A vector DepMap_IDs for which co-essentiality mapping will be performed on. 
#' Will perform pan-cancer analyses (all cell lines) by default, Default: NULL
#' @param core_num integer, Number of cores to run analysis, Default: NULL
#' @param output_dir string, Full path to where output file should be saved, Default: NULL
#' @param data_dir string Path to GRETTA_data
#' @param filename string name of file without the '.csv' extension. 
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
#' gretta_data_dir <- './GRETTA_example/'
#' gretta_output_dir <- './GRETTA_example_output/'
#' 
#' if(!dir.exists(gretta_data_dir)){
#'   download_example_data(".")
#' }
#' 
#' \dontrun{
#' coess_df <- rna_coespress(
#' input_gene = 'ARID1A',
#' input_disease = 'Pancreatic Cancer',
#' core_num = 2,
#' data_dir = gretta_data_dir, 
#' output_dir = gretta_output_dir,
#' test = TRUE)
#' }
#' 
#' @rdname rna_coexpress
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
#' @importFrom stringr str_detect

rna_coexpress <- function(input_genes = NULL, input_disease = NULL,
                          input_cell_lines = NULL, core_num = NULL, output_dir = NULL,
                          data_dir = NULL, filename = NULL, test = FALSE) {
  # Check that essential inputs are given:
  if (is.null(input_genes)) {
    stop("No genes detected")
  }
  if (is.null(output_dir)) {
    output_dir <- paste0(getwd(), "/GRETTA_", Sys.Date())
    message(
      "No output directory specified. Creating: ",
      output_dir
    )
    dir.create(output_dir)
  }
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist. Please provide full path to directory.")
  }
  if (is.null(data_dir)) {
    stop("No directory to data was specified. Please provide path to DepMap data.")
  }
  if (!dir.exists(data_dir)) {
    stop("DepMap data directory does not exists. Please check again and provide the full path to the DepMap data directory.")
  }
  if (!is.null(filename)) {
    output_dir_and_filename <- paste0(
      output_dir,
      "/", filename, ".csv"
    )
  } else {
    output_dir_and_filename <- paste0(
      output_dir,
      "/GRETTA_rna_coexpress_results.csv"
    )
  }
  
  # Detect cores if not defined
  if (is.null(core_num)) {
    cores_detected <- parallel::detectCores()
    message("No cores specified")
    message("Detected: ", cores_detected, " cores")
    message("Using: ", cores_detected / 2, " cores")
    doMC::registerDoMC(cores_detected / 2)
  }
  
  # Load necessary data
  CCLE_exp <- CCLE_exp_annot <- sample_annot <- NULL # see: https://support.bioconductor.org/p/24756/
  load(paste0(data_dir, "/CCLE_exp.rda"), envir = environment())
  load(paste0(data_dir, "/CCLE_exp_annot.rda"), envir = environment())
  load(paste0(data_dir, "/sample_annot.rda"), envir = environment())
  
  # Set cores:
  if (!is.null(core_num)) {
    doMC::registerDoMC(core_num)
  }
  
  # Check if inputs are recognized
  if (!all(input_cell_lines %in% sample_annot$DepMap_ID)) {
    stop(
      input_cell_lines[!input_cell_lines %in% sample_annot$DepMap_ID],
      ", not recognized as a valid sample"
    )
  }
  if (!all(input_genes %in% CCLE_exp_annot$GeneNames)) {
    stop(
      input_genes[!input_genes %in% CCLE_exp_annot$GeneNames],
      ", not recognized. Please check spelling or remove gene name from input"
    )
  }
  
  if (!is.null(input_disease)) {
    selected_cell_lines <- sample_annot %>%
      dplyr::filter(.data$DepMap_ID %in%
                      CCLE_exp$DepMap_ID, .data$disease %in%
                      input_disease) %>%
      dplyr::pull(.data$DepMap_ID)
  } else if(is.null(input_cell_lines) & is.null(input_disease)){
    selected_cell_lines <- sample_annot %>%
      dplyr::pull(.data$DepMap_ID)
  } else if(is.null(input_cell_lines) & is.null(input_disease)){
    selected_cell_lines <- sample_annot %>%
      dplyr::pull(.data$DepMap_ID)
  } else {
    stop(
      "Following may not be available:",
      input_disease, "\n or \n", input_cell_lines
    )
  }
  
  # Provide only expr of genes of interst
  AllGenes <- CCLE_exp_annot$GeneNameID[-1]
  
  All_res <- NULL
  for(g in seq_len(length(input_genes))){
    # g <- 1
    select_gene <- input_genes[g]
    
    Gene_A_expr <- CCLE_exp %>%
      dplyr::select(.data$DepMap_ID, get_GeneNameID(select_gene,
                                                    data_dir = data_dir
      )) %>%
      dplyr::filter(.data$DepMap_ID %in% selected_cell_lines) %>%
      dplyr::rename(RNA_expr = 2)
    
    # Need to define function. A fix for a
    # strange bug:
    `%dopar%` <- foreach::`%dopar%`
    
    # Begin loop
    message("This may take a few mins... Consider running with a higher core numbers to speed up the analysis.")
    if (test == TRUE) {
      run <- 10
    } else {
      run <- length(unique(AllGenes))
    }
    res <- each <- NULL
    res <- foreach::foreach(
      each = seq_len(run),
      .combine = dplyr::bind_rows
    ) %dopar% {
      if (each == 1) {
        message(
          "Processing ", each, " of ",
          length(AllGenes), "\n"
        )
      } else if (each == length(AllGenes)) {
        message(
          "Processing ", each, " of ",
          length(AllGenes), "\n"
        )
      } else if (each %% 1000 == 0) {
        message(
          "Processing ", each, " of ",
          length(AllGenes), "\n"
        )
      }
      
      Gene_B_expr <- CCLE_exp %>%
        dplyr::select(.data$DepMap_ID, AllGenes[each]) %>%
        dplyr::filter(.data$DepMap_ID %in% selected_cell_lines) %>%
        dplyr::rename(RNA_expr = 2)
      
      res_pearson <- cor.test(Gene_A_expr$RNA_expr,
                              Gene_B_expr$RNA_expr,
                              alternative = "two.sided",
                              method = "pearson", na.action = "na.omit"
      ) %>%
        broom::tidy() %>%
        dplyr::mutate(
          GeneNameID_A = get_GeneNameID(select_gene,
                                        data_dir = data_dir
          ),
          GeneNameID_B = AllGenes[each]
        ) %>%
        dplyr::select(
          .data$GeneNameID_A, .data$GeneNameID_B,
          tidyr::everything()
        )
      
      res_pearson
    }
    All_res <- dplyr::bind_rows(All_res, res)
  }
  
  # save and return output
  output <- All_res %>%
    dplyr::arrange(.data$GeneNameID_A) %>%
    dplyr::group_by(.data$GeneNameID_A) %>%
    dplyr::arrange(.data$GeneNameID_A, -.data$estimate, .data$p.value) %>%
    dplyr::mutate(
      Rank = order(-.data$estimate,
                   decreasing = FALSE
      ),
      Padj_BH = p.adjust(.data$p.value,
                         method = "BH", n = (length(.data$p.value))
      )
    ) %>%
    dplyr::select(
      .data$GeneNameID_A, .data$GeneNameID_B,
      .data$estimate, .data$statistic, .data$parameter,
      .data$Rank, .data$p.value, .data$Padj_BH
    ) %>%
    dplyr::ungroup() %>%
    readr::write_csv(file = output_dir_and_filename)
  
  message(
    "RNA co-expression mapping finished. Outputs were also written to: ",
    output_dir_and_filename
  )
  return(output)
}
