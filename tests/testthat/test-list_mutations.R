test_that("listing mutations work", {
  dir <- paste0("../testdata/GRETTA_example/")
  df <- list_mutations(gene = "ARID1A", data_dir = dir)
  expect_equal(nrow(df), 458)
})
