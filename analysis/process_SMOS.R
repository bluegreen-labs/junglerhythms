library(raster)

# list files
files <- list.files("data-raw/SMOS-IC/","*.nc", full.names = TRUE)

# split out meta-data
dates <- as.Date(unlist(lapply(strsplit(basename(files), "_"),"[[",1)),"%Y%m%d")
path <- unlist(lapply(strsplit(basename(files), "_"),"[[",3))
file_info <- data.frame(dates, path, filename = files)

#descending <- stack(files[path == "DES"])
#names(descending) <- dates[path == "DES"]

# TODO fill in the dates which are missing with NA layers
# to make sure that even steps are taken during smoothing
# check previous SMOS routines to address this

ascending <- stack(files[path == "ASC"])
names(ascending) <- dates[path == "ASC"]

arr <- as.array(ascending)

# traverse rows and columsn (along z dimension)
test2 <- apply(arr, c(1,2), mean, na.rm = TRUE)

plot(raster(test2))
