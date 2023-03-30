DepMap_dir <- "/projects/marralab/ytakemon_prj/DepMap/21Q4/Data/"
download.file(url = "https://figshare.com/ndownloader/articles/16924132/versions/1", destfile = paste0(DepMap_dir,"/21Q4.zip"))
unzip(zipfile = paste0(DepMap_dir,"/21Q4.zip"), unzip = "/usr/bin/unzip", exdir = paste0(DepMap_dir,"/"))
