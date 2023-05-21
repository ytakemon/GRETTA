test_that("extracting protein expr info works", {
  dir <- paste0("../testdata/GRETTA_example/")
  
  res <- extract_prot(input_samples = c("ACH-000004"), input_genes = c("CIC"), data_dir = dir)
  expect_equal(signif(res$protein_expr,3), -0.385)
})

