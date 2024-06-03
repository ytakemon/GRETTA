test_that("Test that list_available_cancer_types works", {
  dir <- paste0("../testdata/GRETTA_example/")
  
  # To list all types
  res1 <- list_cancer_types(data_dir = dir)
  expect_equal(length(res1), 33)
  
  # To list all subtypes
  res2 <- list_cancer_subtypes(input_disease = "Lung Cancer", data_dir = dir)
  expect_equal(length(res2), 11)
})
