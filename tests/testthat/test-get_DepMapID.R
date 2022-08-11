test_that("get_DepMapID() converts common cell line names to DepMapID", {
  dir <- "/projects/marralab/ytakemon_prj/DepMap/GRETA_data/20Q1/data"
  
  expect_identical(get_DepMapID("JURKAT", data_dir = dir), "ACH-000995")
})
