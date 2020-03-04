#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
#----- load required libraries -------------------------------------------------#
# library(TSA)
library(tidyverse)
library(plyr)
library(ggplot2)
library(gtools)
library(grid)
library(gridExtra)
# library(circular)
# library(fields)
# library(DescTools)
# library(truncnorm)
#----- source required functions -----------------------------------------------#
source("analysis/timelines.R")
source("analysis/consec_years.R")
source("analysis/fourier_function_species_level.R")
source("analysis/climate_coherence_function.R")
#-------------------------------------------------------------------------------#


#----------------------------------------------------------------------
#--------   Phenology data - species correction Meise   ---------------
#----------------------------------------------------------------------
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
# if(zero_events_remove){
#   df <- df[which(df$value != 0),]
# }

df$join_id <- paste0("R",df$image,"-",df$image_row)
metadata <- read.csv("data/phenology_archives_species_long_format_20190626.csv",
                     header = TRUE, sep = ",")
metadata$join_id <- paste(metadata$image,metadata$row, sep = "-")

# test merge the two tables based upon the unique join ID
data <- merge(df, metadata, by = c("join_id"), all.x = TRUE)
data$species_full <- paste(data$genus_Meise, data$species_Meise)

# remove column id.x and rename id.y to id (--> in id.y, empty ids are renamed to EK1, EK2, etc...)
data$id <- as.character(data$id.y)
data = data[,!(names(data) %in% c("id.x","id.y"))]
# data$id <- as.character(data$id)
# data$species_full <- as.character(data$species_full)

# remove rows with NA's in year -> individuals with 'no_data' in the archive
data <- data[!(is.na(data$year)),]
#----------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------------
#--- get selected species and clean time series -----------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
# overview <- read.csv("data/species_meta_data.csv",
#                      header = TRUE,
#                      sep = ",",
#                      stringsAsFactors = FALSE)
# species_list <- overview$species_full
species_list <-c("Afzelia bipindensis","Panda oleosa","Albizia adianthifolia",
                 "Chrysophyllum africanum","Entandrophragma cylindricum")

# for selected species and phenophase: get extended timelines at ID level with 2 year-gaps filled with zero
timelines_id <- two_year_gaps(data = data,
                            species_name = species_list,
                            pheno = "leaf_turnover")
# for fourier, no gaps in timelines allowed
# get longest consecutive timeline; at species level
timelines_sp_consec <- consecutive_timeline_sp(data = timelines_id,
                             species_name = species_list,
                             pheno = "leaf_turnover")

#----------------------------------------------------------------------------------------------------------------------
#--- fourier at the species level -------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
fourier_sp <- peak_detection(data = timelines_sp_consec,
                            species_name = species_list,
                            pheno = "leaf_turnover",
                            perio_plot = FALSE)


#----------------------------------------------------------------------------------------------------------------------
#--- COHERENCE TEST with climate data at the species level ------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
# read climate data
climate <- read.csv("data/yangambi_km5_monthly_kasongo.csv")
climate$date <- paste(climate$year,climate$month,"15",sep = "-")
climate$date <- as.Date(climate$date, "%Y-%m-%d")
climate$date_monthly <- format(as.Date(climate$date), "%Y-%m")

climate.avg <- read.csv("data/ClimData_monthly_avg.csv",header = TRUE,sep = ",")
climate.avg = climate.avg[,(names(climate.avg) %in% c("Month","insol_JR","tmax_JR"))]
colnames(climate.avg)[1] <- "month"
climate <- merge(climate, climate.avg, by = "month", all.x = TRUE)
climate <- climate %>%
  dplyr::arrange(date)

# coherence analysis
coherence <- climate_coherence(data = timelines_sp_consec,
                               climate = climate,
                               species_name = species_list,
                               pheno = "leaf_turnover")
