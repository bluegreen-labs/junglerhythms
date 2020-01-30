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
  filter(species_full %in% "Erythrophleum suaveolens") %>%
  filter(phenophase %in% "leaf_dormancy")
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


####################################################
#### Initial observation of the Periodogram #######
###################################################



# Spectrum function from Bush
spec_fun <- function(x) spectrum(x,#spans = 3,
                                 plot=F,demean=T,detrend=T) #spectrum function for normal smoother periodogram


# ids with to few years or events
a <- data.frame(total_years = (tapply(data_subset$value, data_subset$id, length))/48,
                total_events = tapply(data_subset$value, data_subset$id, sum))
a$id <- rownames(a)
a <- a[!(a$total_years < 2 | a$total_events < 4),]
ids_ok <- a$id

data_subset <- data_subset[(data_subset$id %in% ids_ok),]

# data as timeseries
data_ls <- tapply(data_subset$value, data_subset$id, FUN = function(x)ts(x, start = 1937, frequency = 48))
str(data_ls)

# fourier transform
Fourier_df <- data.frame()
IDs<-names(data_ls)
for (i in 1:length(data_ls)){
  d <- data.frame(freq = spec_fun(data_ls[[i]])$freq,
                  spec = spec_fun(data_ls[[i]])$spec,
                  spec_norm = spec_fun(data_ls[[i]])$spec*(1/mean(spec_fun(data_ls[[i]])$spec)),
                  ID = as.factor(IDs[i]))
  Fourier_df<-rbind(Fourier_df,d)
}

# raw data plot
p_lin <- ggplot(data_subset) +
  geom_line(aes(x = date,
                y = value,
                colour = id)) +
  scale_x_date(date_breaks = "1 years",
               date_labels = "%Y",
               limits = as.Date(c('1937-01-01','1956-12-31'))) +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey50", size = 0.3)) +
  facet_wrap(~id,
             ncol=1)
#Periogorams for each individual
periodogram <- ggplot(Fourier_df) +
  geom_line(aes(x = freq,
                y = spec_norm,
                group = ID,
                colour=ID)) +
  theme_minimal(base_size=14)+
  ylab("Power (standardised spectrum)")+
  xlab("Cycle frequency") +
  labs(title = unique(data_subset$species_full))

grid.arrange(periodogram,p_lin, heights = c(1,1))


#Summarise dominant peaks
Fourier_df_sample_summary <- ddply(Fourier_df, .(ID), summarize,
                                   maxFreq = freq[which.max(spec)]/48,
                                   maxSpec_norm=max(spec_norm))

#Dominant phenology cycle for sample
Median_cycle <- 1/median(Fourier_df_sample_summary$maxFreq)
Median_cycle



#####################################################################################
#### Periodogram analysis with confidence intervals to find frequency and phase #####
#####################################################################################


#Function to calculate the null continuum (super-smoothed spectrum of data) for null hypothesis test
spec_null_fun <- function(x) spectrum(x,spans=c(55,57), # a large span from Bush
                                      plot=F,demean=T,detrend=T) #spectrum function for null hypothesis spectrum



#Function to extract key variables from spectrum outputs
confidence_fun <- function(p){
  ts <- data_ls[[p]]
  d <- data.frame(ID=IDs[p])
  d$length <- length(ts)
  d$freq_dom <- (spec_fun(ts)$freq[which.max(spec_fun(ts)$spec)])/48 #frequency of the dominant peak
  d$cycle_dom <- 1/d$freq_dom
  d$spec_dom <- max(spec_fun(ts)$spec) #spectrum of dominant peak
  d$spec_dom_norm =  d$spec_dom*(1/mean(spec_fun(ts)$spec)) #normalised spectrum so that mean power is equal to 1 enabling comparisons of "Power" between individuals
  dfree <- spec_fun(ts)$df #degrees of freedom of the smoothed periodogram
  lower_ci <- (dfree*d$spec_dom)/(qchisq(c(0.975),dfree)) #lower CI of dominant peak of smoothed periodogram
  spec_null_dom <- spec_null_fun(ts)$spec[which(abs((spec_null_fun(ts)$freq)/48-d$freq_dom)==min(abs((spec_null_fun(ts)$freq)/48-d$freq_dom)))] #the corresponding dominant peak from the null continuum
  d$sig <- lower_ci>spec_null_dom #If the null spectrum is lower than the lower confidence interval, dominant peak can be considered "signficantly different" to null hypothesis
  data.frame(d)
}

## Run loop to apply spectrum functions to each individual time series in data_ls
confidence_results <- ldply(1:length(data_ls),.fun=confidence_fun)

#Categorise dominant cycles and screen results for positive false results where dominant peak is same as full length of the data - "No cyclicity"
confidence_results$cycle_category<-as.factor(ifelse(1/confidence_results$freq_dom > confidence_results$length/2,"No cyclicity",
                                                    ifelse(1/confidence_results$freq_dom > 49,"Supra-annual",
                                                           ifelse(1/confidence_results$freq_dom >= 47,"Annual","Sub-annual"))))

#Is the dominant cycle both significant and represent a regular cycle?
confidence_results$sig_cycle <- as.factor(ifelse(confidence_results$cycle_category!="No cyclicity" & confidence_results$sig==TRUE,TRUE,FALSE))




#####################################################################################
### Find phase of significant dominant cycles and assess synchrony between individuals
#####################################################################################
first_year_events <- data_subset %>%
  # filter(value > 0) %>%
  group_by(species_full, id) %>%
  dplyr::summarise(first_year = min(year))

# Function for co-Fourier analysis of a simulated cosine curve against the empirical data
spec_cross_fun <- function(x)spec.pgram(x,#spans=spans_lookup$spans_smooth[[which.min(abs(spans_lookup$months-length(x)))]],
                                        demean=T,detrend=T,plot=F)
#Function to extract data from co-Fourier anaysis regarding phase difference in radians
#Simulate a "template" time series using the sample modal cycle peaking at the beginning of January 1937 to act as guide in co-Fourier analysis
phase_fun <- function(p){
  w[p] = 2*pi*(1/confidence_results$cycle_dom[p]) #wavelength in radians
  simulated_phase_ts[p] <- ts(4*cos(w[p]*1:360),start=c(first_year_events$first_year[p]),freq=48) #
  ts_cross <- ts.intersect(simulated_phase_ts[p],data_ls[[p]])
  d<-data.frame(ID=IDs[p])
  # d$cycle_used <- w[p]
  d$phase_radians<-spec_cross_fun(ts_cross)$phase[which.max(spec_cross_fun(ts_cross)$spec[,1])] #the phase difference (in radians) at the dominant frequency
  d
}

w = 2*pi*(1/144.000000)
testsimulated_phase_ts <- ts(cos(w*1:360),start=c(1938),freq=48)
plot(testsimulated_phase_ts)#, xlim = c(1938,1970))

## Run to apply the co-fourier spectrum function to each individual time series
# in combination with simulated guide time series for sample modal cycle to identify relative phase difference for each
#the phase difference (in radians) at the dominant frequency
phase_results <- ldply(1:length(data_ls),.fun=phase_fun)
# different formats of timing
phase_results$phase_degree <- deg(phase_results$phase_radians)
phase_results$phase_degree <-ifelse(phase_results$phase_degree < 0, phase_results$phase_degree +360, phase_results$phase_degree)
phase_results$phase_week <- phase_results$phase_degree*48/360
phase_results$phase_doy <- round((phase_results$phase_week*7.6)-7)
phase_results$timing_peak <- as.Date(phase_results$phase_doy, origin = "1938-01-01")
phase_results$timing_peak <- format.Date(phase_results$timing_peak, "%m-%d")

##Combine confidence and phase results
results<-merge(confidence_results,phase_results,"ID")


