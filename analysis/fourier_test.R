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



###################################
#combine species-specific lists of individual time-series

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
data = data[,!(names(data) %in% "id.x")]
data <- data %>%
  rename("id" = id.y)
data$id <- as.character(data$id)
data$species_full <- as.character(data$species_full)
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#-----subset for testing -------------------------------------
#----------------------------------------------------------------------
data_subset <- data %>%
  filter(species_full %in% "Anonidium mannii") %>%
  filter(phenophase %in% "flowers") %>%
  filter(id %in% "304")
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#----- converting date and sorting
#----------------------------------------------------------------------
# convert date
data_subset$date <- as.Date(paste(data_subset$year,
                                  round((data_subset$week*7.6)-7),sep="-"), "%Y-%j")
# sort dataframe according to date
data_subset <- data_subset %>%
  dplyr::arrange(date)
# add a month column
data_subset$month <- format.Date(data_subset$date, "%m")

# raw data plot
p_lin <- ggplot(data_subset,
                aes(x = date,
                    y = value)) +
  geom_line() +
  scale_x_date(date_breaks = "1 years",
               date_labels = "%Y",
               limits = as.Date(c('1938-01-01','1951-12-31'))) +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey50", size = 0.3))

p_lin

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

#----------------------------------------------------------------------
#-----subset for testing -------------------------------------
#----------------------------------------------------------------------
data_subset <- data %>%
  filter(species_full %in% "Anonidium mannii") %>%
  filter(phenophase %in% "flowers")

### Data format introduction ###

# To work this code, empirical or simulated data should be arranged initially as a named list of
# individual time series (called "data_ls"). These should each be of class “time series” and
# have no gaps (any small gaps, for example up to 3 data points long, could be
# interpolated using a simple linear estimator). This will then also be arranged
# into data frame format (data_df) for some graphics / analyses.



test <- ts(data = data_subset$value, start = 1, end = 624, frequency = 1)

#Summarise phenology scores across sample for each month and year
#(to give proportions for boxplots)
data_df_sample_summary <- ddply(data_subset,.(year,month),summarize,
                              propPhenophase = length(value[value > 0])/length(id))

data_df_sample_summary$month <- mapvalues(data_df_sample_summary$month,from=c("01","02","03","04","05","06","07","08","09","10","11","12"), to=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
data_df_sample_summary$month <- factor(data_df_sample_summary$month,ordered=T,levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

#Write circular boxplot graphic
Phenology_circular_boxplot<-ggplot(data_df_sample_summary,aes(x=factor(month),y=propPhenophase))+
  geom_boxplot(fill="orange",outlier.colour="lightgrey")+
  coord_polar(start=0)+
  theme_minimal(base_size=16)+
  theme(legend.position="none",axis.title.x=element_blank(),axis.ticks=element_blank())+
  ylab("")
Phenology_circular_boxplot
##Figure 1b - timeseries plot of raw data for each individual
raw_plots <- ggplot(data_subset,aes(x=date,y=value))+ #,colour=ID
  geom_line()+
  scale_y_continuous(labels=NULL)+
  ylab("Canopy score")+
  xlab("")+
  theme_light(base_size=14)+
  theme(legend.position="none",strip.background = element_blank(), strip.text = element_blank())
  # facet_wrap(~ID,ncol=1)
raw_plots


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

# for (i in 1:length(data_ls)){
#   d<-data.frame(freq=spec_fun(data_ls[[i]])$freq/12)
#   d$spec=spec_fun(data_ls[[i]])$spec
#   d$spec_norm=spec_fun(data_ls[[i]])$spec*(1/mean(spec_fun(data_ls[[i]])$spec))
#   d$ID=as.factor(IDs[i])
#   Fourier_df<-rbind(Fourier_df,d)
# }

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

