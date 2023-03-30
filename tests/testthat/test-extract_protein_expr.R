test_that("extracting protein expr info works", {
  dir <- "/projects/marralab/ytakemon_prj/DepMap/GRETA_data/20Q1/data"
  
  res <- extract_protein_expr(input_samples = c("ACH-000004"), input_genes = c("ATM"), data_dir = dir)
  expect_equal(signif(res$protein_expr,3), -0.334)
})

