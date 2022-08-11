test_that("listing mutations work", {
  dir <- "/projects/marralab/ytakemon_prj/DepMap/GRETA_data/20Q1/data"
  
  df <- list_available_mutations(Gene = "TP53", data_dir = dir)
  expect_equal(nrow(df), 1298)
})
