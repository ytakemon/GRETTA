test_that("searching mutations work", {
  df <- search_mutations(Gene = "TP53")
  expect_equal(nrow(df), 539)
})
