test_that("Extraction for 20Q1 data works", {
  dir <- paste0("../testdata/GRETTA_example/")

  res <- extract_rna(input_samples = c("ACH-001642"), input_genes = c("ARID1A"), data_dir = dir)
  expect_equal(res[ ,2][[1]] %>% signif(., 3), 2.88)
})
