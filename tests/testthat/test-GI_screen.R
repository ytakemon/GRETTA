test_that("GI_screen works", {
  dir <- paste0("../testdata/GRETTA_example/")
  out <- paste0("../testdata/GRETTA_example_output/")
  dir.create(out)
  
  Screen_results <- GI_screen(
    control_id = c("ACH-001354", "ACH-000274", "ACH-001799"), 
    mutant_id = c("ACH-000911", "ACH-001957", "ACH-000075"), 
    gene_list = c('ARID1A', 'ARID1B', 'SMARCA2'),
    core_num = 2, 
    output_dir = out,
    data_dir = dir,
    test = TRUE)
  
  expect_equal(signif(Screen_results$Interaction_score[1], 3), 0.155)
})
