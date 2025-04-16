DepMap_dir <- "./DepMap/RNAi/raw_data/"
download.file(url = "https://figshare.com/ndownloader/articles/9170975/versions/1", destfile = paste0(DepMap_dir,"/RNAi.zip"))
unzip(zipfile = paste0(DepMap_dir,"/23Q4.zip"), unzip = "/usr/bin/unzip", exdir = paste0(DepMap_dir))