DepMap_dir <- "./DepMap/22Q3/Data/"
download.file(url = "https://figshare.com/ndownloader/articles/22765112/versions/4", destfile = paste0(DepMap_dir, "/23Q2.zip"))
unzip(zipfile = paste0(DepMap_dir,"/23Q2.zip"), unzip = "/usr/bin/unzip", exdir = paste0(DepMap_dir,"/"))