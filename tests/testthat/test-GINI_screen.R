test_that("GINIR_screen works", {
  dir <- "/projects/marralab/ytakemon_prj/DepMap/GINIR_data/data"
  
  Screen_results <- GINI_screen(
    control_IDs = c("ACH-001354", "ACH-000274", "ACH-001799"), 
    mutant_IDs = c("ACH-000911", "ACH-001957", "ACH-000075"), 
    core_num = 2, 
    output_dir = getwd(),
    data_dir = dir,
    test = TRUE)
  
  expect_equal(signif(Screen_results$Interaction_score[1], 3), -0.398)
})
