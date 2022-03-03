test_that("get_GeneNameID() converts hugo symbols or NCBI gene IDs to DepMap gene IDs", {
  dir <- "/projects/marralab/ytakemon_prj/DepMap/GINIR_data/data"
  
  expect_identical(get_GeneNameID("A1CF", data_dir = dir), "A1CF_29974")
  expect_identical(get_GeneNameID("29974", data_dir = dir), "A1CF_29974")
})
