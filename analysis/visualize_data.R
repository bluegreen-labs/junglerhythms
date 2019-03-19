# plot things

# create summary graphs for all species with Plant List accepted names
# combine a linear and circular plots

source('/data/Dropbox/Research_Projects/code_repository/bitbucket/junglerhythms/R/linear_plot.R')
source('/data/Dropbox/Research_Projects/code_repository/bitbucket/junglerhythms/R/circular_plot.r')

# load libraries
library(tidyverse)
library(ggthemes)
library(gridExtra)

# read in the original
original_data <- readRDS("data/jungle_rhythms_weekly_annotations.rds")

# read in the meta-data
meta_data <- read.table("data/species_meta_data.csv",
                        header = TRUE, sep = ",")

meta_data <- meta_data[which(tolower(meta_data$genus) == "millettia"),]

# list full species name
species_list <- paste(meta_data$genus,
                      meta_data$species)

lapply(species_list[1:length(species_list)], function(species){

  message(sprintf("processing: %s", species))

  p1 <- try(linear_plot(data = original_data,
                        species_name = species))
  p2 <- try(circle_plot(data = original_data,
                        threshold = 0.1,
                        species_name = species))

  if(inherits(p1,"try-error") || inherits(p2,"try-error")){
    warning(paste("failed to process:", species))
    return(NULL)
  }

  filename <- sprintf("~/Desktop/tmp/%s.png",tolower(species))
  filename <- gsub(" ","_", filename)

  # write stuff to file
  png(filename,1000,500)
    grid.arrange(p1, p2, ncol = 2)
  dev.off()
})
