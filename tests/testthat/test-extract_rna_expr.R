test_that("Extraction for 20Q1 data works", {
  res <- extract_rna_expr(Input_samples = c("ACH-001642"), Input_genes = c("TP53"))
  
  expect_equal(res[ ,2][[1]] %>% signif(., 3), 0.516)
})
