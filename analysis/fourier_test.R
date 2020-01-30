#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
#-------------------------------------------------------------------------------#
# libraries
library(TSA)
library(tidyverse)
library(plyr)
library(ggplot2)
# library(reshape)
# library(RColorBrewer)
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
data$species_full <- paste(data$genus_Meise, data$species_Meise)

# remove column id.x and rename id.y to id (--> in id.y, empty ids are renamed to EK1, EK2, etc...)
data$id <- as.character(data$id.y)
data = data[,!(names(data) %in% c("id.x","id.y"))]
# data$id <- as.character(data$id)
# data$species_full <- as.character(data$species_full)

# remove rows with NA's in year -> individuals with 'no_data' in the archive
data <- data[!(is.na(data$year)),]
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#-----subset for testing -------------------------------------
#----------------------------------------------------------------------
data_subset <- data %>%
  filter(species_full %in% "Pericopsis elata") %>% # Irvingia grandifolia
  filter(phenophase %in% "leaf_dormancy") %>% #flowers #leaf_dormancy
  filter(id %in% "127") #587
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

# data_subset <- data_subset %>%
#   filter(date > "1947-01-01")

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
#----- figures showing raw data from Bush (same as our circular - linear plots)
#----------------------------------------------------------------------
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
##Figure 1b - timeseries plot of raw data for each individual
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
# facet_wrap(~ID,ncol=1)
grid.arrange(Phenology_circular_boxplot,p_lin, widths = c(1,2))

####################################################
#### Initial observation of the Periodogram #######
###################################################

# Spectrum function from Bush
spec_fun <- function(x) spectrum(x,#spans = 3,
                               plot=F,demean=T,detrend=T) #spectrum function for normal smoother periodogram

# data as timeseries
test <- ts(data = data_subset$value, start = 1938, frequency = 48) # 1938 # time series of the data

# fourier transform
Fourier_df <- data.frame(freq = spec_fun(test)$freq/48,
                         spec = spec_fun(test)$spec,
                         spec_norm = spec_fun(test)$spec*(1/mean(spec_fun(test)$spec)))


##Plot periodograms for visual inspection

#Periogorams for each individual on one plot - Figure 1c
Periodogram <- ggplot(data=Fourier_df) +
  geom_line(aes(x=freq,y=spec_norm))+
  theme_minimal(base_size=14)+
  ylab("Power (standardised spectrum)")+
  xlab("Cycle frequency (cycles per month)")+
  ggtitle("Smoothed periodograms of phenology cycles by individual")

Periodogram

# Summarise dominant peak
maxFreq = Fourier_df$freq[which.max(Fourier_df$spec)]
cycle = 1/maxFreq
cycle # if months then yes correct 24 = 2 years, if observations then no to short


#####################################################################################
#### Periodogram analysis with confidence intervals to find frequency and phase #####
#####################################################################################


#Function to calculate the null continuum (super-smoothed spectrum of data) for null hypothesis test
spec_null_fun <- function(x) spectrum(x,spans=c(55,57), # a large span from Bush
                                    plot=F,demean=T,detrend=T) #spectrum function for null hypothesis spectrum


# Summarise dominant peak
freq_dom <- (spec_fun(test)$freq[which.max(spec_fun(test)$spec)])/48 #frequency of the dominant peak
cycle_dom <- 1/freq_dom
spec_dom <- max(spec_fun(test)$spec) #spectrum of dominant peak
spec_dom_norm =  spec_dom*(1/mean(spec_fun(test)$spec)) #normalised spectrum so that mean power is equal to 1 enabling comparisons of "Power" between individuals
# stats compared to null hypothesis
dfree <- spec_fun(test)$df #degrees of freedom of the smoothed periodogram
lower_ci <- (dfree*spec_dom)/(qchisq(c(0.975),dfree)) #lower CI of dominant peak of smoothed periodogram
spec_null_dom <- spec_null_fun(test)$spec[which(abs((spec_null_fun(test)$freq)/48-freq_dom)==min(abs((spec_null_fun(test)$freq)/48-freq_dom)))] #the corresponding dominant peak from the null continuum

sig <- lower_ci > spec_null_dom #If the null spectrum is lower than the lower confidence interval, dominant peak can be considered "signficantly different" to null hypothesis
# Categorise dominant cycles and screen results for positive false results where dominant peak is same as full length of the data - "No cyclicity"
cycle_category <- as.factor(ifelse(cycle_dom > length(test)/2,"No cyclicity",
                                   ifelse(cycle_dom > 49,"Supra-annual",
                                          ifelse(cycle_dom >= 47,"Annual","Sub-annual"))))
#Is the dominant cycle both significant and represent a regular cycle?
sig_cycle <- as.factor(ifelse(cycle_category!="No cyclicity" & sig==TRUE,TRUE,FALSE))


#####################################################################################
### Find phase of significant dominant cycles and assess synchrony between individuals
#####################################################################################

# Calculate the modal dominant cycle length of the sample
Dominant_cycle <- as.numeric(round(cycle_dom[sig_cycle==TRUE],digits=0))

#Simulate a "template" time series using the sample modal cycle peaking at the beginning of January 1986 to act as guide in co-Fourier analysis
w = 2*pi*(1/Dominant_cycle) # wavelength in radians
simulated_phase_ts <- ts(4*cos(w*1:720),start=c(1938,1),end=c(1949,48),freq=48) #1:360 -> 1:720 -> no interuption in cosinus for longer time series
ts_cross <- ts.intersect(simulated_phase_ts,test)

# Function for co-Fourier analysis of a simulated cosine curve against the empirical data
spec_cross_fun <- function(x)spec.pgram(x,#spans=spans_lookup$spans_smooth[[which.min(abs(spans_lookup$months-length(x)))]],
                                        demean=T,detrend=T,plot=F)

#the phase difference (in radians) at the dominant frequency
crosstest <- spec_cross_fun(ts_cross)
crosstest_df <- data.frame(freq = crosstest$freq,
                         spec1 = crosstest$spec[,1],
                         spec1_norm = crosstest$spec[,1]*(1/mean(crosstest$spec[,1])),
                         spec2 = crosstest$spec[,2],
                         spec2_norm = crosstest$spec[,2]*(1/mean(crosstest$spec[,2])),
                         phase = crosstest$phase
                         )
plot(ts_cross)
##Plot periodograms for visual inspection
Periodogram <- ggplot(data=crosstest_df) +
  geom_line(aes(x=freq,y=spec1_norm))+
  geom_line(aes(x=freq,y=spec2_norm), col = "blue")+
  theme_minimal(base_size=14)+
  ylab("Power (standardised spectrum)")+
  xlab("Cycle frequency")
Periodogram
p_phase <- ggplot(data=crosstest_df) +
  geom_line(aes(x=freq,y=phase))+
  theme_minimal(base_size=14)+
  ylab("Power (standardised spectrum)")+
  xlab("Cycle frequency")
p_phase

phase_radians <- spec_cross_fun(ts_cross)$phase[which.max(spec_cross_fun(ts_cross)$spec[,1])]
# different formats of timing
phase_degree <- deg(phase_radians)
phase_week <- phase_degree*48/360
phase_doy <- round((phase_week*7.6)-7)
timing_peak <- as.Date(phase_doy, origin = "1938-01-01")
timing_peak <- format.Date(timing_peak, "%m-%d")

# return data as data frame
return(data.frame(cycle_dom = cycle_dom,
                  sig = sig,
                  cycle_category = cycle_category,
                  phase_doy = phase_doy,
                  timing_peak = timing_peak
                  ))
