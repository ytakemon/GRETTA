# library(pacman)
# p_load(here, devtools, fs, tidyverse)
# here() # Should show directory of GINIscreenR git repo
# # Setup Pandoc path incase it doesn't exist, migh only apply to me (Yuka)
# if(Sys.getenv("RSTUDIO_PANDOC") == ""){
#   Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio-server/bin/pandoc")
# }
#
# (a <- factor(c("character", "hits", "your", "eyeballs")))
# #> [1] character hits      your      eyeballs
# #> Levels: character eyeballs hits your
# (b <- factor(c("but", "integer", "where it", "counts")))
# #> [1] but      integer  where it counts
# #> Levels: but counts integer where it
# c(a, b)
# #> [1] character hits      your      eyeballs  but       integer
# #> [7] where it  counts
# #> 8 Levels: character eyeballs hits your but counts ... where it
#
# factor(c(as.character(a), as.character(b)))
# #> [1] character hits      your      eyeballs  but       integer
# #> [7] where it  counts
# #> 8 Levels: but character counts eyeballs hits integer ... your
# fbind <- function(a, b) {
#   factor(c(as.character(a), as.character(b)))
# }
#
# use_r("fbind")
# #> • Edit 'R/fbind.R'
# #> • Call `use_test()` to create a matching test file
#
# load_all()
# #> ℹ Loading foofactors
#
# library(GINIscreenR)
#
# a <- factor(c("character", "hits", "your", "eyeballs"))
# b <- factor(c("but", "integer", "where it", "counts"))
#
# fbind(a, b)
# #> [1] character hits      your      eyeballs  but       integer
# #> [7] where it  counts
# #> 8 Levels: but character counts eyeballs hits integer ... your
#
#
# use_testthat()
# #> ✔ Adding 'testthat' to Suggests field in DESCRIPTION
# #> ✔ Setting Config/testthat/edition field in DESCRIPTION to '3'
# #> ✔ Creating 'tests/testthat/'
# #> ✔ Writing 'tests/testthat.R'
# #> • Call `use_test()` to initialize a basic test file and open it for editing.
#
#
# devtools::build_readme()
