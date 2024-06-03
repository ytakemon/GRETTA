DepMap_dir <- "./DepMap/23Q4/"
download.file(url = "https://plus.figshare.com/ndownloader/articles/24667905/versions/2", destfile = paste0(DepMap_dir,"/23Q4.zip"))
unzip(zipfile = paste0(DepMap_dir,"/23Q4.zip"), unzip = "/usr/bin/unzip", exdir = paste0(DepMap_dir))