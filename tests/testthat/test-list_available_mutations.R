test_that("listing mutations work", {
  dir <- "/projects/marralab/ytakemon_prj/DepMap/GINIR_data/data"
  
  df <- list_available_mutations(Gene = "TP53", data_dir = dir)
  expect_equal(nrow(df), 1298)
})
