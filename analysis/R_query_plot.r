# demo query to python from R using google_earth_engine_subset
#
# (please make sure to have the below libraries installed)
# This script clones the project in your home directory,
# you can remove these files afterward.

# load libraries
library(ggplot2)

# change this depending on system settings
python_path = "/usr/bin/python"

# clone the gee_subset project
#relies on git being installed
# and will work out of the box for most
# on OSX or Linux.
#
# basic gee_subset requirements apply
# mainly, having a working GEE python API install
setwd("~")
#system("git clone https://github.com/khufkens/google_earth_engine_subsets.git")
path = "~/gee_subset/gee_subset/"

# set product parameters, such as
# product name, band(s) to query, start and end date of the range
# and the lcoation
#product = "MODIS/006/MCD43A4"
#band = "Nadir_Reflectance_Band1 Nadir_Reflectance_Band2 Nadir_Reflectance_Band4 Nadir_Reflectance_Band3"

product = "COPERNICUS/S2"
band = "B1 B2 B3 B4 QA60"

# store output in the R temporary directory
directory = tempdir()

sites <- read.table("~/yangambi_sites.csv", sep = ",", header = TRUE)

years <- c(2015:2018)

output <- do.call("rbind", apply(sites,1,function(site){
  lapply(years, function(year){

    start_date = paste0(year, "-01-01")
    end_date = paste0(year, "-12-31")
    location = paste(site['lat'], site['lon'])

    # make the gee_subset.py python call
    # time the duration of the call for reporting
    system(sprintf("%s %s/gee_subset.py -p %s -b %s -s %s -e %s -l %s -d %s -sc 10 -pd 0.1",
                 python_path,
                 path,
                 product,
                 band,
                 start_date,
                 end_date,
                 location,
                 directory
                 ), wait = TRUE)

    # read in the data stored in the temporary directory
    df = read.table( paste0( directory, "/site_",
                           tail( unlist( strsplit( product, "[/]" ) ), n=1 ),
                           "_gee_subset.csv" ), sep = ",", header = TRUE )

    # add site name
    df$site <- site['site']
    return(df)
  })
}))

output <- do.call("rbind",output)

write.table(output, "Yangambi_Sentinel_MIX_plots.csv",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = ",")
