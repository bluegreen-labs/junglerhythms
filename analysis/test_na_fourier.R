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

species_name = "Albizia adianthifolia" #Erythrophleum suaveolens , Irvingia grandifolia
pheno = "leaf_turnover"
data_subset <- data %>%
  filter(species_full %in% species_name) %>%
  filter(phenophase %in% pheno)

# convert date
data_subset$date <- as.Date(paste(data_subset$year,
                                  round((data_subset$week*7.6)-7),sep="-"), "%Y-%j")
# add a month column
data_subset$month <- format.Date(data_subset$date, "%m")

# sort dataframe according to date
data_subset <- data_subset %>%
  dplyr::arrange(id,date)

unique(data_subset$id)

data_ind <- data_subset %>%
  filter(id %in% "353")
a <- unique(data_ind$year) # 1943 missing



# grow dataset to full range
years <- sort(rep(min(data_ind$year):max(data_ind$year), 48))
b <- as.numeric(unique(years))
days <- rep(round((1:48 * 7.6) - 7),length(unique(years)))
df_extended <- data.frame(date = as.Date(paste(years, days, sep = "-"), "%Y-%j"))

data_ind_full <- merge(data_ind, df_extended, by = "date", all.y = TRUE)
data_ind_full$species_full <- unique(na.omit(data_ind_full$species_full)) #fill species in empty years
data_ind_full$id <- unique(na.omit(data_ind_full$id))
data_ind_full$year <- format.Date(data_ind_full$date, "%Y")

raw_plots <- ggplot(data_ind_full,
                    aes(x = date,
                        y = value))+
  geom_line()+
  scale_y_continuous(labels=NULL)+
  ylab("Observations") +
  xlab("") +
  ggtitle(paste(species_name, "-", pheno)) +
  theme_light(base_size = 14) +
  theme(legend.position = "none") +
  facet_wrap(~id,
             ncol = 1,
             strip.position = "left")
raw_plots

# group consecutive missing years
missing_years <- setdiff(b,a)
consec_years <- cumsum(c(1, abs(missing_years[-length(missing_years)] - missing_years[-1]) > 1))
missing_years <- as.data.frame(cbind(missing_years,consec_years))
colnames(missing_years)[1] <- 'year'
test <- missing_years %>%
  group_by(consec_years) %>%
  dplyr::summarise(lgh_consec_years = length(consec_years))
missing_years <- merge(missing_years, test, by = "consec_years", all.x = TRUE)
data_ind_full <- merge(data_ind_full, missing_years, by = "year", all.x = TRUE)

# if only 2 consecutive years of missing data -> fill with zeros. If longer, keep NA
data_ind_full$value <- ifelse(is.na(data_ind_full$value) & data_ind_full$lgh_consec_years < 3 , 0, data_ind_full$value)

raw_plots <- ggplot(data_ind_full,
                    aes(x = date,
                        y = value))+
  geom_line()+
  scale_y_continuous(labels=NULL)+
  ylab("Observations") +
  xlab("") +
  ggtitle(paste(species_name, "-", pheno)) +
  theme_light(base_size = 14) +
  theme(legend.position = "none") +
  facet_wrap(~id,
             ncol = 1,
             strip.position = "left")
raw_plots


### repeat to only keep longest timeline for fourier analysis
# remove rows with NA's in values -> individuals with 'no_data' in the archive
timelines <- data_ind_full[!(is.na(data_ind_full$value)),]
years_timeline <- as.numeric(unique(timelines$year))
consec_timelines <- cumsum(c(1, abs(years_timeline[-length(years_timeline)] - years_timeline[-1]) > 1))
years_timeline <- as.data.frame(cbind(years_timeline,consec_timelines))
colnames(years_timeline)[1] <- 'year'
test <- years_timeline %>%
  group_by(consec_timelines) %>%
  dplyr::summarise(lgh_consec_timelines = length(consec_timelines))
years_timeline <- merge(years_timeline, test, by = "consec_timelines", all.x = TRUE)

timelines <- merge(timelines, years_timeline, by = "year", all.x = TRUE)

timelines$value <- ifelse(timelines$lgh_consec_timelines == max(timelines$lgh_consec_timelines) , timelines$value, NA)
timelines <- timelines[!(is.na(timelines$value)),]



raw_plots <- ggplot(timelines,
                    aes(x = date,
                        y = value))+
  geom_line()+
  scale_y_continuous(labels=NULL)+
  ylab("Observations") +
  xlab("") +
  ggtitle(paste(species_name, "-", pheno)) +
  theme_light(base_size = 14) +
  theme(legend.position = "none") +
  facet_wrap(~id,
             ncol = 1,
             strip.position = "left")
raw_plots



first_year <- as.numeric(format.Date(timelines$date[1], "%Y"))
data_ts_full <- ts(data = timelines$value, start = c(first_year,1), frequency = 48) #, end = c(last_year,48)


spec_fun <- function(x) spectrum(x, #spans = c(3,3), # adding spans changes significance of dominant frequencies found!!!
                                 plot=F, demean=T, detrend=T)





# fourier transform
fourier_df_ind_full <- data.frame(freq = spec_fun(data_ts_full)$freq/48,
                             spec = spec_fun(data_ts_full)$spec,
                             spec_norm = spec_fun(data_ts_full)$spec*(1/mean(spec_fun(data_ts_full)$spec)),
                             freq_null = spec_null_fun(data_ts_full)$freq/48,
                             spec_null = spec_null_fun(data_ts_full)$spec,
                             spec_norm_null = spec_null_fun(data_ts_full)$spec*(1/mean(spec_null_fun(data_ts_full)$spec)))


Periodograms_full <- ggplot(fourier_df_ind_full) +
  geom_line(aes(x = freq,
                y = spec_norm)) +
  geom_line(aes(x = freq_null,
                y = spec_norm_null),
            color = "blue") +
  scale_y_continuous(labels=NULL) +
  ylab("Power (standardised spectrum)") +
  xlab("Cycle frequency") +
  ggtitle("") +
  theme_light(base_size=14)+
  theme(legend.position = "none")
Periodograms_full


freq_dom = (spec_fun(data_ts_full)$freq[which.max(spec_fun(data_ts_full)$spec)])/48 # frequency of the dominant peak
cycle_dom = 1/freq_dom
cycle_dom
spec_dom = max(spec_fun(data_ts_full)$spec) # spectrum of dominant peak
spec_dom_norm =  spec_dom*(1/mean(spec_fun(data_ts_full)$spec)) # normalised spectrum so that mean power is equal to 1 enabling comparisons of "Power" between individuals
dfree = spec_fun(data_ts_full)$df # degrees of freedom of the smoothed periodogram
lower_ci = (dfree*spec_dom)/(qchisq(c(0.975),dfree)) # lower CI of dominant peak of smoothed periodogram
spec_null_dom <- spec_null_fun(data_ts_full)$spec[which(abs((spec_null_fun(data_ts_full)$freq)/48-freq_dom)==min(abs((spec_null_fun(data_ts_full)$freq)/48-freq_dom)))] # the corresponding dominant peak from the null continuum
sig = lower_ci > spec_null_dom # If the null spectrum is lower than the lower confidence interval, dominant peak can be considered "signficantly different" to null hypothesis
sig
cycle_category <- as.factor(ifelse(1/freq_dom > length(data_ts_full)/2,"No cyclicity",
                                   ifelse(1/freq_dom > 49,"Supra-annual",
                                          ifelse(1/freq_dom >= 47,"Annual","Sub-annual"))))
# Is the dominant cycle both significant and represent a regular cycle?
sig_cycle <- as.factor(ifelse(cycle_category!="No cyclicity" & sig==TRUE,TRUE,FALSE))



