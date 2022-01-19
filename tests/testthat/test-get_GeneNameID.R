test_that("get_GeneNameID() converts hugo symbols or NCBI gene IDs to DepMap gene IDs", {
  expect_identical(get_GeneNameID("A1CF"), "A1CF_29974")
  expect_identical(get_GeneNameID("29974"), "A1CF_29974")
})
