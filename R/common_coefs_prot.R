#' @title Perform Pearson coefficient mapping of all pairs of genes in a protein co-expression analysis
#' 
#' @description Performs Pearson correlation coefficient analyses between all gene pairs found in a protein co-expression analysis. 
#' Input parameters should match those used to run `annotate_df()`.
#' 
#' @param input_genes string, A vector containing one or more "Hugo Symbol", Default: NULL
#' @param input_disease string, A vector one or more disease contexts, Will perform pan-cancer analyses 
#' (all cell lines) by default, Default: NULL
#' @param input_cell_lines string, A vector DepMap_IDs for which co-essentiality mapping will be performed on. 
#' Will perform pan-cancer analyses (all cell lines) by default, Default: NULL
#' @param core_num integer, Number of cores to run analysis, Default: NULL
#' @param output_dir string, Full path to where output file should be saved, Default: NULL
#' @param data_dir string Path to GRETTA_data
#' @param filename string name of file without the '.csv' extension. 
#' @param test logical, TRUE/FALSE whether you want to run only a small subset (first 10 genes) to ensure function will run properly 
#' prior to running all genes. Default: FALSE.
#'
#' @return A data frame containing Pearson correlation coefficients. A copy is also saved to the 
#' directory defined in `output_dir`.
#' 
#' @details Description of output data frame
#' * `GeneNameID_A` - Hugo symbol of query gene.
#' * `GeneNameID_B` - Hugo symbol of all genes quantified.
#' * `estimate` - Correlation coefficient output from `?cor.test`.
#' * `statistic` - Pearson's correlation statistic. Output from `?cor.test`.
#' * `p.value` - P-value from Pearson's correlation statistic. Output from `?cor.test`.
#' * `parameter` - Degrees of freedom. Output from `?cor.test`.
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
#' coess_df <- common_coefs_prot(
#' input_gene = c("ARID1A", "SMARCB1"),
#' input_disease = 'Pancreatic Cancer',
#' core_num = 2,
#' data_dir = gretta_data_dir, 
#' output_dir = gretta_output_dir,
#' test = TRUE)
#' }
#' 
#' @rdname common_coefs_prot
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

common_coefs_prot <- function(input_genes = NULL, input_disease = NULL,
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
      "/GRETTA_protein_coexpress_results.csv"
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
  protein_annot <- protein_nodup <- sample_annot <- NULL # see: https://support.bioconductor.org/p/24756/
  load(paste0(data_dir, "/sample_annot.rda"), envir = environment())
  load(paste0(data_dir, "/protein_nodup.rda"), envir = environment())
  load(paste0(data_dir, "/protein_annot.rda"), envir = environment())
  
  # Set cores:
  if (!is.null(core_num)) {
    doMC::registerDoMC(core_num)
  }
  
  # Check if inputs are recognized
  if (!all(input_cell_lines %in% sample_annot$DepMap_ID)) {
    stop(
      paste0(input_cell_lines[!input_cell_lines %in% sample_annot$DepMap_ID], collapse = ", "),
      ", not recognized as a valid sample"
    )
  }
  if (!all(input_genes %in% protein_nodup$Gene_Symbol)) {
    stop(
      paste0(input_genes[!input_genes %in% protein_nodup$Gene_Symbol], collapse = ", "),
      ", not recognized or protein expression is not available. Please check spelling or remove gene name from input"
    )
  }
  
  # Define cell lines
  if (!is.null(input_disease)) {
    selected_cell_lines <- sample_annot %>%
      dplyr::filter(.data$DepMap_ID %in%
                      protein_annot$DepMap_ID, .data$disease %in%
                      input_disease) %>%
      dplyr::pull(.data$DepMap_ID)
  } else if (!is.null(input_cell_lines)) {
    selected_cell_lines <- sample_annot %>%
      dplyr::filter(.data$DepMap_ID %in%
                      protein_annot$DepMap_ID, .data$DepMap_ID %in%
                      input_cell_lines) %>%
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
  AllGenes <- input_genes
  All_res <- NULL
  for (g in seq_len(length(input_genes))) {
    # g <- 1
    select_gene <- input_genes[g]
    
    Gene_A_expr <- protein_nodup %>%
      dplyr::filter(.data$Gene_Symbol %in% select_gene) %>%
      dplyr::select(
        .data$Gene_Symbol, .data$Description,
        .data$Uniprot, .data$Uniprot_Acc, dplyr::contains("_TenPx")
      ) %>%
      tidyr::pivot_longer(
        -c(
          .data$Gene_Symbol, .data$Description,
          .data$Uniprot, .data$Uniprot_Acc
        ),
        names_to = "Gygi_ID",
        values_to = "protein_expr"
      ) %>%
      dplyr::left_join(protein_annot, by = c(Gygi_ID = "GygiNames")) %>%
      dplyr::filter(.data$DepMap_ID %in% selected_cell_lines, ) %>%
      dplyr::select(.data$DepMap_ID, .data$Gene_Symbol, .data$protein_expr, .data$Gygi_ID) %>%
      filter(!is.na(.data$protein_expr))
    
    if(nrow(Gene_A_expr) < 3){
      stop(
        "There are only 3 or less lines with ", select_gene, " expression.",
        "Please remove from input_gene"
      )
    }
    
    # Need to define function. A fix for a
    # strange bug:
    `%dopar%` <- foreach::`%dopar%`
    
    # Begin loop
    message("Analysing: ", select_gene, ". \n")
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
      
      Gene_B_expr <- protein_nodup %>%
        dplyr::filter(.data$Gene_Symbol %in% AllGenes[each]) %>%
        dplyr::select(
          .data$Gene_Symbol, .data$Description,
          .data$Uniprot, .data$Uniprot_Acc, dplyr::contains("_TenPx")
        ) %>%
        tidyr::pivot_longer(
          -c(
            .data$Gene_Symbol, .data$Description,
            .data$Uniprot, .data$Uniprot_Acc
          ),
          names_to = "Gygi_ID",
          values_to = "protein_expr"
        ) %>%
        dplyr::left_join(protein_annot, by = c(Gygi_ID = "GygiNames")) %>%
        dplyr::filter(.data$DepMap_ID %in% selected_cell_lines, ) %>%
        dplyr::select(.data$DepMap_ID, .data$Gene_Symbol, .data$protein_expr, .data$Gygi_ID)
      
      # Check if enough Ns
      check <- Gene_B_expr %>% dplyr::filter(!is.na(.data$protein_expr))
      if(nrow(check) < 3){
        res_pearson <- tibble(
          GeneNameID_A = select_gene,
          GeneNameID_B = AllGenes[each],
          statistic = NA,
          p.value = NA,
          parameter = NA,
          conf.low = NA,
          conf.high = NA, 
          method = NA, 
          alternative = NA
        )
        
        res_pearson
      } else {
        res_pearson <- cor.test(Gene_A_expr$protein_expr,
                                Gene_B_expr$protein_expr,
                                alternative = "two.sided",
                                method = "pearson", na.action = "na.omit"
        ) %>%
          broom::tidy() %>%
          dplyr::mutate(
            GeneNameID_A = select_gene,
            GeneNameID_B = AllGenes[each]
          ) %>%
          dplyr::select(
            .data$GeneNameID_A, .data$GeneNameID_B,
            tidyr::everything()
          )
        
        res_pearson
      }
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
    filter(GeneNameID_A != GeneNameID_B) %>%
    arrange(estimate,statistic,p.value, GeneNameID_A) %>%
    filter(row_number() %% 2 == 1) %>% ## Select odd rows to remove duplicates
    arrange(GeneNameID_A, GeneNameID_B) %>%
    readr::write_csv(file = output_dir_and_filename)
  
  message(
    "Protein co-expression mapping finished. Outputs were also written to: ",
    output_dir_and_filename
  )
  return(output)
}
