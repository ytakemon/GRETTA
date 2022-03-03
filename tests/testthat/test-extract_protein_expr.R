test_that("extracting protein expr info works", {
  dir <- "/projects/marralab/ytakemon_prj/DepMap/GINIR_data/data"
  
  res <- extract_protein_expr(Input_samples = c("ACH-000004"), Input_genes = c("ATM"), data_dir = dir)
  expect_equal(signif(res$protein_expr,3), -0.334)
})

