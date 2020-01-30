#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
#-------------------------------------------------------------------------------#
# libraries
library(TSA)
library(tidyverse)
library(plyr)
library(ggplot2)
library(reshape)
library(RColorBrewer)
library(gtools)
library(grid)
library(circular)
library(fields)
library(DescTools)
library(truncnorm)
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
data$species_full <- as.character(paste(data$genus_Meise, data$species_Meise))

# remove column id.x and rename id.y to id (--> in id.y, empty ids are renamed to EK1, EK2, etc...)
data$id <- as.character(data$id.y)
data = data[,!(names(data) %in% c("id.x","id.y"))]

# remove rows with NA's in year -> individuals with 'no_data' in the archive
data <- data[!(is.na(data$year)),]
#----------------------------------------------------------------------


# #----------------------------------------------------------------------
# #-------- get the species list ----------------------------------------
# #----------------------------------------------------------------------
overview <- read.csv("data/species_meta_data.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)

# merge number of events to data, so you can filter on minimum number of events
overview <- overview[,(names(overview) %in% c("species_full",
                                              "deciduousness",
                                              "basal_area_site"))]

# merge with data to filter on deciduousness
data <- merge(data, overview, by = "species_full", all.x = TRUE)

# get species_list of species
# - in basal area list
# - with a minimum of 5 events and 5 (5*48=240) years of data
species_list <- data %>%
  filter(basal_area_site > 0) %>%
  filter(phenophase == "leaf_dormancy") %>%
  group_by(species_full) %>%
  dplyr::summarise(count_week = sum(value), # total events
                   total_week = length(value)) # total weeks
species_list <- species_list %>%
  filter(count_week >=5) %>%
  filter(total_week >= 240)
species_list <- species_list$species_full

###########################
# testing species list
species_list <- c("Pericopsis elata","Erythrophleum suaveolens")
###########################
#----------------------------------------------------------------------
#-----subset for testing -------------------------------------
#----------------------------------------------------------------------
data_subset <- data %>%
  filter(phenophase %in% "leaf_dormancy") %>%
  filter(species_full %in% species_list)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#----- converting date and sorting
#----------------------------------------------------------------------
# convert date
data_subset$date <- as.Date(paste(data_subset$year,
                                  round((data_subset$week*7.6)-7),
                                  sep="-"), "%Y-%j")

which(is.na(data_subset$year))

# sort dataframe according to date
data_subset <- data_subset %>%
  dplyr::arrange(date)
# # add a month column
# data_subset$month <- format.Date(data_subset$date, "%m")

# species means
data_summ <- data_subset %>%
  group_by(species_full, year, date) %>%
  dplyr::summarise(mean_value = mean(value, na.rm=TRUE))


test <- data_summ %>%
  group_by(species_full) %>%
  dplyr::summarise(count_events = sum(mean_value), # total events
                   total_week = length(mean_value),
                   total_years = total_week/48,
                   start = min(date),
                   end = max(date))
sp_years <- data_summ %>%
  group_by(species_full,year) %>%
  dplyr::summarise(count_events = sum(mean_value), # total events
                   total_week = length(mean_value))
# ---> missing years ...




spec_fun <- function(x) spectrum(x,#spans = 3,
                                 plot=T,demean=T,detrend=T) #spectrum function for normal smoother periodogram

ftest <- data_summ %>%
  group_by(species_full) %>%
  do(freq_list = matrix(spec_fun(.$mean_value)$spec)) %>%
  unnest(freq_list)

tapply(ftest$freq_list, ftest$species_full, length)

#--------- needs time-series of equal length for direct comparison
#--------- now the number of frequencies is different

# ftest <- ftest %>%
#   tidyr::spread(key = species_full,value = freq_list)
#
# spread(ftest, species_full, freq_list)


#,
     #spec_list = spec_fun(.$mean_value)$spec)


Fourier_df<-data.frame()
for (i in 1:length(data_ls)){
  d<-data.frame(freq=spec_fun(data_ls[[i]])$freq/12)
  d$spec=spec_fun(data_ls[[i]])$spec
  d$spec_norm=spec_fun(data_ls[[i]])$spec*(1/mean(spec_fun(data_ls[[i]])$spec))
  d$ID=as.factor(IDs[i])
  Fourier_df<-rbind(Fourier_df,d)
}




#----------------------------------------------------------------------
#----------------------------------------------------------------------
# this simple example from
# https://anomaly.io/detect-seasonality-using-fourier-transform-r/
# works!
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# compute the Fourier Transform
# p = periodogram(raw$Visite)
# p = periodogram(data_subset$value)
#
# dd = data.frame(freq=p$freq, spec=p$spec)
# order = dd[order(-dd$spec),]
# top2 = head(order, 2)
#
# # display the 2 highest "power" frequencies
# top2
#
# # convert frequency to time periods
# time = 1/top2$freq
# time_year = time/48
# time_year
# time_months = time/4
#----------------------------------------------------------------------
#----------------------------------------------------------------------




####################################################
#### Intitial observation of the Periodogram #######
###################################################

###To display smoothed periodograms - Figure 1c (main text)

##Parameters ##

# Spans for the Daniel kernel to smooth the spectrum The following
# list allows the spec_fun function to find appropriate spans for the Daniel
# kernel to successively apply to the spectrum to give a smoothed periodogram
# with bandwidth similar to 0.1 (span size is linked to the length of the
# oringal timeseries data, we give options here for data from 24 to 360 months
# long).

##Functions ##

#Spectrum fucntion
spec_fun <- function(x) spectrum(x,#spans = 3,
                                 plot=F,demean=T,detrend=T) #spectrum function for normal smoother periodogram

##Run fourier analysis using function "spectrum" on each list object (individual
##time series) to give Fourier outputs for each individual

fourier <- spec_fun(test)

Fourier_df <- data.frame()
d <- data.frame(freq = fourier$freq)
d$spec = fourier$spec
d$spec_norm = fourier$spec*(1/mean(fourier$spec))
Fourier_df <- rbind(Fourier_df,d)



##Plot periodograms for visual inspection

#Periogorams for each individual on one plot - Figure 1c
Periodograms<-ggplot(data=Fourier_df,aes(x=freq,y=spec_norm)) +
  geom_line()+
  theme_minimal(base_size=14)+
  ylab("Power (standardised spectrum)")+
  xlab("Cycle frequency (cycles per month")+
  ggtitle("Smoothed periodograms of phenology cycles by individual")

Periodograms

# #Summarise dominant peaks
# Fourier_df_sample_summary<-ddply(Fourier_df, .(ID), summarize,
#                                  maxFreq=freq[which.max(spec)],
#                                  maxSpec_norm=max(spec_norm))
#Dominant phenology cycle for sample
# Median_cycle<-1/median(Fourier_df_sample_summary$maxFreq)
# Median_cycle<-1/median(Fourier_df_sample_summary$maxFreq)
# Median_cycle
maxFreq = Fourier_df$freq[which.max(Fourier_df$spec)]
cycle = 1/maxFreq
cycle

