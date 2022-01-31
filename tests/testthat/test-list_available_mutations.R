test_that("listing mutations work", {
  df <- list_available_mutations(Gene = "TP53")
  expect_equal(nrow(df), 539)
})
