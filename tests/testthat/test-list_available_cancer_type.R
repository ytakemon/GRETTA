test_that("Test that list_available_cancer_types works", {
  
  # To list all types
  res1 <- list_available_cancer_types()
  expect_equal(length(res1), 37)
  
  # To list all subtypes
  res2 <- list_available_cancer_subtypes("Lung Cancer")
  expect_equal(length(res2), 11)
})
