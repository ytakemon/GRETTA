test_that("Number of example control and mutant groups are correct", {
  dir <- paste0("../testdata/GRETTA_example/")
  
  test <- select_cell_lines(input_gene = "ARID1A", data_dir = dir) %>% dplyr::count(Group)
  test_control <- test %>% dplyr::filter(Group == "Control")
  test_amp <- test %>% dplyr::filter(Group == "Others")
  
  expect_equal(test_control$n, 758)
  expect_equal(test_amp$n, 73)
})