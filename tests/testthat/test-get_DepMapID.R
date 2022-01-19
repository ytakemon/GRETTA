test_that("get_DepMapID() converts common cell line names to DepMapID", {
  expect_identical(get_DepMapID("JURKAT"), "ACH-000995")
})
