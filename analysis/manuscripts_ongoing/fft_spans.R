#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
#-------------------------------------------------------------------------------#
# libraries
library(tidyverse)
library(plyr)
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

#----
# span detection at full data resolution of 48 observations per year
#-----
species_name = "Scorodophloeus zenkeri"
pheno = "flowers"
data_subset <- data %>%
  filter(species_full %in% species_name) %>%
  filter(phenophase %in% pheno)

data_subset$date <- as.Date(paste(data_subset$year,
                                  round((data_subset$week*7.6)-7),sep="-"), "%Y-%j")

# average by date
data_subset <- data_subset %>%
  group_by(species_full, date) %>%
  dplyr::summarise(mean_value = mean(value),
                   scaled_value = ifelse(any(value > 0), 1, 0))
data_subset$year <- format.Date(data_subset$date, "%Y")
# sort dataframe according to date
data_subset <- data_subset %>%
  dplyr::arrange(date)
#----------------------------------------------------------------------
#----------------------------------------------------------------------

# reduce years by one consecutively and get bandwidth to ~0.05 by changing spans
data_reduced <- data_subset %>%
  filter(date > '1950-12-31')

first_year <- as.numeric(format.Date(data_reduced$date[1], "%Y"))
last_year <- as.numeric(format.Date(last(data_reduced$date), "%Y"))
data_ts <- ts(data = data_reduced$mean_value, start = c(first_year,1), frequency = 48) #, end = c(last_year,48)

spec_fun <- function(x) spectrum(x, spans = NULL,#c(3,3),
                                 plot=F,demean=T,detrend=T)
spec_fun(data_ts)$bandwidth
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
spans_lookup <- list(observations = c(96,144, 192, 240, 288, 336, 384, 432, 480, 528,576,624,672,720,768,816,864,912,960, 1008),
                     spans_smooth = list(NULL, NULL, NULL, NULL, # smoothed spectrum of data -> smoothed periodogram bandwith to approx. 0.05
                                         NULL, NULL, NULL, NULL,
                                         3,3,3,3,
                                         3,3,3,
                                         c(3,3), c(3,3), c(3,3), c(3,3), c(3,3)),
                     # spans_smooth=list(NULL, NULL, NULL, NULL, # smoothed spectrum of data -> smoothed periodogram bandwith to approx. 0.1
                     #                   3,3,3,3,
                     #                   c(3,3), c(3,3), c(3,3),
                     #                   c(3,5), c(3,5), c(3,5),
                     #                   c(5,5), c(5,5), c(5,5),
                     #                   c(5,7), c(5,7), c(5,7)),
                     spans_super_smooth = list(c(5,7), c(7,9), c(11,11), c(13,13), # super-smoothed spectrum of data -> smoothed periodogram bandwith to approx. 1
                                               c(15,17), c(19,19), c(19,21), c(23,23),
                                               c(25,27), c(27,29), c(29,31),
                                               c(33,33), c(35,37), c(37,39),
                                               c(39,41), c(45,45), c(45,47),
                                               c(49,51), c(49,51), c(49,51)))
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------



#----
# span detection at full data resolution of monthly observations per year
#-----

# convert date
data_subset$yr_month <- format(data_subset$date, "%Y-%m")

# average by year_month
data_monthly_res <- data_subset %>%
  group_by(yr_month) %>%
  dplyr::summarise(monthly_value = mean(scaled_value),
                   yr = min(year))
# sort dataframe according to date
data_monthly_res <- data_monthly_res %>%
  dplyr::arrange(yr_month)

#------------------
# reduce years by one consecutively and get bandwidth to ~0.05 by changing spans
data_monthly_res <- data_monthly_res %>%
  filter(yr > 1955)

# data as timeseries, species_level
first_year <- as.numeric(data_monthly_res$yr[1])
data_ts <- ts(data = data_monthly_res$monthly_value, start = first_year, frequency = 12)

spec_fun <- function(x) spectrum(x, spans = c(5,5),
                                 plot=F,demean=T,detrend=T)
spec_fun(data_ts)$bandwidth
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
spans_lookup_monthly <- list(observations = c(24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 180, 192, 204, 216, 228, 240),
                     # spans_smooth = list(NULL, NULL, NULL, NULL, # smoothed spectrum of data -> smoothed periodogram bandwith to approx. 0.05
                     #                     NULL, NULL, NULL,
                     #                     3,3,3,3,
                     #                     3,3,3,
                     #                     c(3,3), c(3,3), c(3,3), c(3,3), c(3,3)),
                     spans_smooth=list(NULL, NULL, NULL, NULL, # smoothed spectrum of data -> smoothed periodogram bandwith to approx. 0.1
                                       3,3,3,
                                       c(3,3), c(3,3), c(3,3), c(3,3),
                                       c(3,5), c(3,5), c(3,5),
                                       c(5,5), c(5,5), c(5,5),
                                       c(5,7), c(5,7)),
                     spans_super_smooth = list(c(5,7), c(9,9), c(11,11), c(13,13),
                                               c(15,17), c(19,21), c(21,21), c(23,23),
                                               c(25,27), c(29,29), c(29,31), c(33,35),
                                               c(37,39), c(37,39), c(39,41), c(45,45),
                                               c(45,45), c(49,51), c(49,51)))
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
