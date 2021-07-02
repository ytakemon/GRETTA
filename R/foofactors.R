library(pacman)
p_load(here, devtools, fs, tidyverse)
here() # Should show directory of GINIscreenR git repo
# Setup Pandoc path incase it doesn't exist, migh only apply to me (Yuka)
if(Sys.getenv("RSTUDIO_PANDOC") == ""){
  Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio-server/bin/pandoc")
}
