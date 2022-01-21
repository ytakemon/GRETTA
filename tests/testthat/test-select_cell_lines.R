test_that("Number of example control and mutant groups are correct", {
  test <- select_cell_lines(Input_gene = "A1CF") %>% dplyr::count(Group)
  test_control <- test %>% dplyr::filter(Group == "Control")
  test_amp <- test %>% dplyr::filter(Group == "Amplified")
  
  expect_equal(test_control$n, 610)
  expect_equal(test_amp$n, 35)
})
