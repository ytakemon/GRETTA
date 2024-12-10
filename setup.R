
# Install the GitHub package
devtools::install_github("Marralab/GRETTA", ref="tire_kicking")

#install depmap data
if (!file.exists("data/dataset.csv")) {
    download.file("https://example.com/dataset.csv", destfile = "data/dataset.csv", mode = "wb")
}

# Test the package
library(GRETTA)
cat("GitHub package installed and loaded successfully.\n")
