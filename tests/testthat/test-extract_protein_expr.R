test_that("extracting protein expr info works", {
  res <- extract_protein_expr(Input_samples = c("ACH-000004"), Input_genes = c("ATM"))
  expect_equal(signif(res$protein_expr,3), -0.334)
})

