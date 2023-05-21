DepMap_dir <- "./DepMap/22Q2/Data/"
download.file(url = "https://figshare.com/ndownloader/articles/19700056/versions/2", destfile = paste0(DepMap_dir,"/22Q2.zip"))
unzip(zipfile = paste0(DepMap_dir,"/22Q2.zip"), unzip = "/usr/bin/unzip", exdir = paste0(DepMap_dir,"/"))