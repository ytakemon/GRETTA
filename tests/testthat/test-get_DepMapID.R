test_that("get_DepMapID() converts common cell line names to DepMapID", {
  dir <- paste0("../testdata/GRETTA_example/")
  
  expect_identical(get_DepMapID("JURKAT", data_dir = dir), "ACH-000995")
})
