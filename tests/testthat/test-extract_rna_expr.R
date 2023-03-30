test_that("Extraction for 20Q1 data works", {
  dir <- "/projects/marralab/ytakemon_prj/DepMap/GRETA_data/20Q1/data"
  
  res <- extract_rna_expr(input_samples = c("ACH-001642"), input_genes = c("TP53"), data_dir = dir)
  expect_equal(res[ ,2][[1]] %>% signif(., 3), 0.516)
})
