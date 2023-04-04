test_that("extracting protein expr info works", {
  dir <- "/projects/marralab/ytakemon_prj/DepMap/GRETTA_data/20Q1/data"
  
  res <- extract_prot(input_samples = c("ACH-000004"), input_genes = c("ATM"), data_dir = dir)
  expect_equal(signif(res$protein_expr,3), -0.334)
})

