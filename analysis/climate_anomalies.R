#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
# graphics.off()
#----- load required packages --------------------------------------------------#
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(tseries)

#----------------------------------------------------------------------
#-------- input that you can change   ---------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------
drought_mm = 150
#----------------------------------------------------------------------



# read climate data
climate <- read.csv("data/yangambi_km5_precip_monthly_kasongo.csv")
climate <- climate %>%
  filter(year > 1936,
         year < 1959,
         measurement == "precipitation") %>%
  # mutate(col = ifelse(value < drought_mm,"dry","wet"))
  mutate(col = ifelse(month %in% c(1,2,6,7,12), "dry", "wet"))
climate$date <- paste(climate$year,climate$month,"15",sep = "-")
climate$date <- as.Date(climate$date, "%Y-%m-%d")
climate$yr_month <- format(as.Date(climate$date), "%Y-%m")

# # check with yearly summaries of climate data
# climate_check <- read.csv("data/yangambi_km5_precip_yearly_kasongo.csv")
# climate_check <- climate_check %>%
#   filter(year > 1936,
#          year < 1960)
# test <- as.data.frame(tapply(climate$value,list(climate$year),sum))
# colnames(test)[1] <- "precip_sum"
# test$year <- rownames(test)
# test <- merge(climate_check,test,by = "year")
# test$diff <- test$precip_sum - test$value
# summary(test$diff)

# remove years where data is inconsistent with the yearly summary
# climate$value <- ifelse(climate$year %in% c("1940","1943","1949"), NA, climate$value)

#----------------------------------------------------------------------
# difference between average year and actual year --> anomalies study
#----------------------------------------------------------------------
# average year
mmonthly <- as.data.frame(tapply(climate$value,list(climate$month),mean, na.rm = TRUE))
colnames(mmonthly)[1] <- "mmonthly"
mmonthly$month <- rownames(mmonthly)
climate <- merge(climate,mmonthly,by = "month")
# difference
climate$diff <- climate$value - climate$mmonthly
# rename value to precip
climate <- climate %>%
  rename("precip" = value)
#----------------------------------------------------------------------



p_climate <- ggplot() +
  geom_col(data = climate, aes(x = date,
                               y = diff,
                               colour = col,
                               fill = col)) +
  scale_colour_manual(values = c("lightcoral","lightblue"),
                      aesthetics = c("fill","colour")) +
  labs(y = "precip. (mm/month)",
       x = "Year") +
  scale_x_date(date_breaks = "1 year",
               date_labels = "%Y",
               limits = as.Date(c('1935-01-01','1960-01-01'))) +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="none"
  )







#----------------------------------------------------------------------
#--------   Phenology data - species correction Meise   ---------------
#----------------------------------------------------------------------
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
# df <- df[which(df$value != 0),]
df$join_id <- paste0("R",df$image,"-",df$image_row)

metadata <- read.csv("data/phenology_archives_species_long_format_20190619.csv",
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
#----------------------------------------------------------------------

# #----------------------------------------------------------------------
# #-------- read in census data  ----------------------------------------
# #-------- get species-specific basal area at plot and site level ------
# #----------------------------------------------------------------------
# # read in census data and convert basal area
# census <- read.csv("data/YGB_ForestPlotsNET_corrected_indet.csv",
#                    header = TRUE,
#                    sep = ",",
#                    stringsAsFactors = FALSE)
# census <- census %>%
#   rename("species_full" = Species)
# census$species_full <- ifelse(census$species_full == "Unknown", NA, census$species_full)
#
# # remove individuals without C1DBH4, these are new recruits for census2 + understory if understory.remove = TRUE
# # and calculate basal_area for each individual
# census <- census[!(is.na(census$C1DBH4)),]
# census$basal_area = pi*(census$C1DBH4/2000)^2
# # only keep mixed plots
# # summarize at species level
# census <- census %>%
#   filter(grepl("MIX", Plot)) %>%
#   group_by(species_full) %>%
#   dplyr::summarise(basal_area_site = sum(basal_area)/5,
#                    abundance = length(basal_area)/5,
#                    dbh_max_mm_census = max(C1DBH4))
# basal_area_total <- sum(census$basal_area_site)
# census$basal_area_percentage <- census$basal_area_site/basal_area_total *100
# #----------------------------------------------------------------------
# # read in trait data from Steven Janssens
# traits <- read.csv("data/Dataset_traits_African_trees.csv",
#                    header = TRUE,
#                    sep = ",",
#                    stringsAsFactors = FALSE)
# #----------------------------------------------------------------------
# census <- merge(census, traits, by = "species_full", all.x = TRUE)


overview <- read.csv("data/species_meta_data.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)
species_list <- overview$species_full
#------------------------------------------------------------------------
# filtering species for analysis
#------------------------------------------------------------------------
# species_list <- census %>%
#   filter(deciduousness == "deciduous")
# species_list <- species_list$species_full
#
# # subset on species
# data_subset <- data %>%
#   filter(species_full %in% "Scorodophloeus zenkeri",#species_list,  Erythrophleum suaveolens, Irvingia grandifolia, Diospyros iturensis
#          phenophase == "leaf_dormancy") # leaf_dormancy
#
#
# # convert date to year-month --> to be able to match with monthly climate data
# data_subset$date <- as.Date(
#   paste(data_subset$year,
#         round(data_subset$week*7.6),sep="-"), "%Y-%j")
# data_subset$yr_month <- format(as.Date(data_subset$date), "%Y-%m")
#
# # phenological data to monthly resolution
# data_subset <- data_subset %>%
#   group_by(yr_month) %>%
#   dplyr::summarise(mean_value = mean(value),
#                    sum_value = sum(value),
#                    nr_indiv = length(value))
#
# # remove if data = NA
# data_subset <- data_subset[!(is.na(data_subset$yr_month)),]
#
# #------------------------------------------------------------------------
# # merge with climate data
# #------------------------------------------------------------------------
# data_subset <- merge(data_subset, climate, by = "yr_month", all.x = TRUE)

#------------------------------------------------------------------------
# cross-correlation function
#------------------------------------------------------------------------
## check if time-series is stationary: https://towardsdatascience.com/cross-correlation-of-currency-pairs-in-r-ccf-d27eec2d4b91
# adf.test(data_subset$precip) # p value lower then 0.05 = stationary
# kpss.test(data_subset$precip) # p value higher then 0.05 = stationary
# adf.test(data_subset$diff)
# kpss.test(data_subset$diff)
## Check for autocorrelation in a timeseries
# autocorr <- acf(data_subset$precip)
# autocorr <- acf(data_subset$diff)

# adf.test(data_subset$mean_value)
# kpss.test(data_subset$mean_value)
#
#
#
# corr <- ccf(data_subset$precip, data_subset$mean_value, calc.ci = TRUE, pl=TRUE) # precip as first value: to check if something in precip has an effect on phenology at time t
# ci <- 0.95
# ci_value <- qnorm((1 + ci)/2)/sqrt(corr$n.used)
# #abline(h = ci_value, col = "red")
# corr <- ccf(data_subset$diff, data_subset$mean_value, pl=TRUE)



#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------- DORMANCY x actual precipitation --------------------------------------------
#----------------------------------------------------------------------------------------

output_dormancy <- data.frame(matrix(0,nrow=length(species_list),ncol=19))
colnames(output_dormancy)[1] <- "species_full"
colnames(output_dormancy)[2:14] <- c("t-6","t-5","t-4","t-3","t-2","t-1","t","t+1","t+2","t+3","t+4","t+5","t+6")
colnames(output_dormancy)[15] <- "ci"

ci <- 0.95

for(i in 1:length(species_list)){
  data_sp <- data %>%
    filter(species_full %in% species_list[i], #"Staudtia kamerunensis",
           phenophase == "leaf_dormancy")
  # convert date to year-month --> to be able to match with monthly climate data
  data_sp$date <- as.Date(
    paste(data_sp$year,
          round(data_sp$week*7.6),sep="-"), "%Y-%j")
  data_sp$yr_month <- format(as.Date(data_sp$date), "%Y-%m")

  # phenological data to monthly resolution
  data_sp <- data_sp %>%
    group_by(yr_month) %>%
    dplyr::summarise(mean_value = mean(value),
                     sum_value = sum(value),
                     nr_indiv = length(value))

  # remove if data = NA
  data_sp <- data_sp[!(is.na(data_sp$yr_month)),]


  # merge with climate data
  data_sp <- merge(data_sp, climate, by = "yr_month", all.x = TRUE)

  # cross-correlation + confidence interval
  if(max(data_sp$mean_value)>0){
    corr.precip <- ccf(data_sp$precip, data_sp$mean_value,lag = 6, pl = FALSE)
    # plot(corr.precip, main = paste(species_list[i], " - dormancy x precip"))
    ci_value <- qnorm((1 + ci)/2)/sqrt(corr.precip$n.used)
  } else {
    corr.precip$acf <- "NA"
    ci_value <- "NA"
  }

  output_dormancy[i,1] <- species_list[i]
  output_dormancy[i,2:14] <- as.numeric(corr.precip$acf)
  output_dormancy[i,15] <- as.numeric(ci_value)
}

# maximum correlation before and after time 0 of phenophase
colnames(output_dormancy)[16] <- "max.corr.before"
colnames(output_dormancy)[17] <- "max.corr.after"
absmax <- function(x) { x[which.max( abs(x) )][1]}
output_dormancy$max.corr.before <- apply(output_dormancy[,2:7], 1, absmax) # if lag-time changes, this should be adapted --> this is 6 lags max
output_dormancy$max.corr.after <- apply(output_dormancy[,9:14], 1, absmax)
# only keep significant correlations
output_dormancy$max.corr.before <- ifelse(abs(output_dormancy$max.corr.before) > output_dormancy$ci, output_dormancy$max.corr.before, "NA")
output_dormancy$max.corr.after <- ifelse(abs(output_dormancy$max.corr.after) > output_dormancy$ci, output_dormancy$max.corr.after, "NA")

# get the timing of the correlations
colnames(output_dormancy)[18] <- "max.corr.before.timing"
colnames(output_dormancy)[19] <- "max.corr.after.timing"
max.corr.before.timing <- apply(output_dormancy[,2:7], 1, function(x) { x[which.max( abs(x) )]}) # if lag-time changes, this should be adapted --> this is 6 lags max
max.corr.before.timing <- names(unlist(max.corr.before.timing))
max.corr.after.timing <- apply(output_dormancy[,9:14], 1, function(x) { x[which.max( abs(x) )]})
max.corr.after.timing <- names(unlist(max.corr.after.timing))

# remove NA's, because this messes up the colnames of the timing
output_dormancy <- output_dormancy[complete.cases(output_dormancy), ]
output_dormancy$max.corr.before.timing <- max.corr.before.timing
output_dormancy$max.corr.after.timing <- max.corr.after.timing
# only keep timing of significant correlations
output_dormancy$max.corr.before.timing <- ifelse(output_dormancy$max.corr.before %in% "NA", "NA", output_dormancy$max.corr.before.timing)
output_dormancy$max.corr.after.timing <- ifelse(output_dormancy$max.corr.after %in% "NA", "NA", output_dormancy$max.corr.after.timing)



dormancy_precip <- output_dormancy[,(names(output_dormancy) %in% c("species_full",
                                                                   "max.corr.before",
                                                                   "max.corr.after",
                                                                    "max.corr.before.timing",
                                                                    "max.corr.after.timing"))]
colnames(dormancy_precip) <- c("species_full",
                               "ccf_dormancy_precip_lead",
                               "ccf_dormancy_precip_lag",
                               "ccf_dormancy_precip_lead_timing",
                               "ccf_dormancy_precip_lag_timing")

# change the order of the columns
dormancy_precip <- dormancy_precip[c("species_full",
                                     "ccf_dormancy_precip_lead",
                                     "ccf_dormancy_precip_lead_timing",
                                     "ccf_dormancy_precip_lag",
                                     "ccf_dormancy_precip_lag_timing")]


#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------- TURNOVER x actual precipitation --------------------------------------------
#----------------------------------------------------------------------------------------

output_turnover <- data.frame(matrix(0,nrow=length(species_list),ncol=19))
colnames(output_turnover)[1] <- "species_full"
colnames(output_turnover)[2:14] <- c("t-6","t-5","t-4","t-3","t-2","t-1","t","t+1","t+2","t+3","t+4","t+5","t+6")
colnames(output_turnover)[15] <- "ci"

ci <- 0.95

for(i in 1:length(species_list)){
  data_sp <- data %>%
    filter(species_full %in% species_list[i], #"Staudtia kamerunensis",
           phenophase == "leaf_turnover")
  # convert date to year-month --> to be able to match with monthly climate data
  data_sp$date <- as.Date(
    paste(data_sp$year,
          round(data_sp$week*7.6),sep="-"), "%Y-%j")
  data_sp$yr_month <- format(as.Date(data_sp$date), "%Y-%m")

  # phenological data to monthly resolution
  data_sp <- data_sp %>%
    group_by(yr_month) %>%
    dplyr::summarise(mean_value = mean(value),
                     sum_value = sum(value),
                     nr_indiv = length(value))

  # remove if data = NA
  data_sp <- data_sp[!(is.na(data_sp$yr_month)),]


  # merge with climate data
  data_sp <- merge(data_sp, climate, by = "yr_month", all.x = TRUE)

  # cross-correlation + confidence interval
  if(max(data_sp$mean_value)>0){
    corr.precip <- ccf(data_sp$precip, data_sp$mean_value,lag = 6, pl = FALSE)
    # plot(corr.precip, main = paste(species_list[i], " - turnover x precip"))
    ci_value <- qnorm((1 + ci)/2)/sqrt(corr.precip$n.used)
  } else {
    corr.precip$acf <- "NA"
    ci_value <- "NA"
  }

  output_turnover[i,1] <- species_list[i]
  output_turnover[i,2:14] <- as.numeric(corr.precip$acf)
  output_turnover[i,15] <- as.numeric(ci_value)
}

# maximum correlation before and after time 0 of phenophase
colnames(output_turnover)[16] <- "max.corr.before"
colnames(output_turnover)[17] <- "max.corr.after"
absmax <- function(x) { x[which.max( abs(x) )][1]}
output_turnover$max.corr.before <- apply(output_turnover[,2:7], 1, absmax) # if lag-time changes, this should be adapted --> this is 6 lags max
output_turnover$max.corr.after <- apply(output_turnover[,9:14], 1, absmax)
# only keep significant correlations
output_turnover$max.corr.before <- ifelse(abs(output_turnover$max.corr.before) > output_turnover$ci, output_turnover$max.corr.before, "NA")
output_turnover$max.corr.after <- ifelse(abs(output_turnover$max.corr.after) > output_turnover$ci, output_turnover$max.corr.after, "NA")

# get the timing of the correlations
colnames(output_turnover)[18] <- "max.corr.before.timing"
colnames(output_turnover)[19] <- "max.corr.after.timing"
max.corr.before.timing <- apply(output_turnover[,2:7], 1, function(x) { x[which.max( abs(x) )]}) # if lag-time changes, this should be adapted --> this is 6 lags max
max.corr.before.timing <- names(unlist(max.corr.before.timing))
max.corr.after.timing <- apply(output_turnover[,9:14], 1, function(x) { x[which.max( abs(x) )]})
max.corr.after.timing <- names(unlist(max.corr.after.timing))
# remove NA's, because this messes up the colnames of the timing
output_turnover <- output_turnover[complete.cases(output_turnover), ]
output_turnover$max.corr.before.timing <- max.corr.before.timing
output_turnover$max.corr.after.timing <- max.corr.after.timing
# only keep timing of significant correlations
output_turnover$max.corr.before.timing <- ifelse(output_turnover$max.corr.before %in% "NA", "NA", output_turnover$max.corr.before.timing)
output_turnover$max.corr.after.timing <- ifelse(output_turnover$max.corr.after %in% "NA", "NA", output_turnover$max.corr.after.timing)



turnover_precip <- output_turnover[,(names(output_turnover) %in% c("species_full",
                                                                   "max.corr.before",
                                                                   "max.corr.after",
                                                                   "max.corr.before.timing",
                                                                   "max.corr.after.timing"))]
colnames(turnover_precip) <- c("species_full",
                               "ccf_turnover_precip_lead",
                               "ccf_turnover_precip_lag",
                               "ccf_turnover_precip_lead_timing",
                               "ccf_turnover_precip_lag_timing")

# change the order of the columns
turnover_precip <- turnover_precip[c("species_full",
                                     "ccf_turnover_precip_lead",
                                     "ccf_turnover_precip_lead_timing",
                                     "ccf_turnover_precip_lag",
                                     "ccf_turnover_precip_lag_timing")]



#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------- DORMANCY x anomalies in precipitation --------------------------------------
#----------------------------------------------------------------------------------------

output_dormancy_anom <- data.frame(matrix(0,nrow=length(species_list),ncol=41))
colnames(output_dormancy_anom)[1] <- "species_full"
colnames(output_dormancy_anom)[2:38] <- c("t-18","t-17","t-16","t-15","t-14","t-13",
                                          "t-12","t-11","t-10","t-9","t-8","t-7",
                                          "t-6","t-5","t-4","t-3","t-2","t-1","t",
                                          "t+1","t+2","t+3","t+4","t+5","t+6",
                                          "t+7","t+8","t+9","t+10","t+11","t+12",
                                          "t+13","t+14","t+15","t+16","t+17","t+18")
colnames(output_dormancy_anom)[39] <- "ci"

ci <- 0.95

for(i in 1:length(species_list)){
  data_sp <- data %>%
    filter(species_full %in% species_list[i],
           phenophase == "leaf_dormancy")
  # convert date to year-month --> to be able to match with monthly climate data
  data_sp$date <- as.Date(
    paste(data_sp$year,
          round(data_sp$week*7.6),sep="-"), "%Y-%j")
  data_sp$yr_month <- format(as.Date(data_sp$date), "%Y-%m")

  # phenological data to monthly resolution
  data_sp <- data_sp %>%
    group_by(yr_month) %>%
    dplyr::summarise(mean_value = mean(value),
                     sum_value = sum(value),
                     nr_indiv = length(value))

  # remove if data = NA
  data_sp <- data_sp[!(is.na(data_sp$yr_month)),]


  # merge with climate data
  data_sp <- merge(data_sp, climate, by = "yr_month", all.x = TRUE)

  # cross-correlation + confidence interval
  if(max(data_sp$mean_value)>0){
    corr.anom <- ccf(data_sp$diff, data_sp$mean_value,lag = 18, pl = FALSE)
    # plot(corr.anom, main = paste(species_list[i], " - dormancy x anom"))
    ci_value <- qnorm((1 + ci)/2)/sqrt(corr.anom$n.used)
  } else {
    corr.anom$acf <- "NA"
    ci_value <- "NA"
  }

  output_dormancy_anom[i,1] <- species_list[i]
  output_dormancy_anom[i,2:38] <- as.numeric(corr.anom$acf)
  output_dormancy_anom[i,39] <- as.numeric(ci_value)
}

# maximum correlation before time 0 of phenophase
colnames(output_dormancy_anom)[40] <- "max.corr.before"
absmax <- function(x) { x[which.max( abs(x) )][1]}
output_dormancy_anom$max.corr.before <- apply(output_dormancy_anom[,2:19], 1, absmax) # if lag-time changes, this should be adapted --> this is 6 lags max
# only keep significant correlations
output_dormancy_anom$max.corr.before <- ifelse(abs(output_dormancy_anom$max.corr.before) > output_dormancy_anom$ci, output_dormancy_anom$max.corr.before, "NA")

# get the timing of the correlations
colnames(output_dormancy_anom)[41] <- "max.corr.before.timing"
max.corr.before.timing <- apply(output_dormancy_anom[,2:19], 1, function(x) { x[which.max( abs(x) )]}) # if lag-time changes, this should be adapted --> this is 6 lags max
max.corr.before.timing <- names(unlist(max.corr.before.timing))
# remove NA's, because this messes up the colnames of the timing
output_dormancy_anom <- output_dormancy_anom[complete.cases(output_dormancy_anom), ]
output_dormancy_anom$max.corr.before.timing <- max.corr.before.timing
# only keep timing of significant correlations
output_dormancy_anom$max.corr.before.timing <- ifelse(output_dormancy_anom$max.corr.before %in% "NA", "NA", output_dormancy_anom$max.corr.before.timing)

dormancy_anom <- output_dormancy_anom[,(names(output_dormancy_anom) %in% c("species_full",
                                                                           "max.corr.before",
                                                                           "max.corr.before.timing"))]
colnames(dormancy_anom) <- c("species_full",
                             "ccf_dormancy_anom_lead",
                             "ccf_dormancy_anom_lead_timing")

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------- TURNOVER x anomalies in precipitation --------------------------------------
#----------------------------------------------------------------------------------------
output_turnover_anom <- data.frame(matrix(0,nrow=length(species_list),ncol=41))
colnames(output_turnover_anom)[1] <- "species_full"
colnames(output_turnover_anom)[2:38] <- c("t-18","t-17","t-16","t-15","t-14","t-13",
                                          "t-12","t-11","t-10","t-9","t-8","t-7",
                                          "t-6","t-5","t-4","t-3","t-2","t-1","t",
                                          "t+1","t+2","t+3","t+4","t+5","t+6",
                                          "t+7","t+8","t+9","t+10","t+11","t+12",
                                          "t+13","t+14","t+15","t+16","t+17","t+18")
colnames(output_turnover_anom)[39] <- "ci"

ci <- 0.95

for(i in 1:length(species_list)){
  data_sp <- data %>%
    filter(species_full %in% species_list[i],
           phenophase == "leaf_turnover")
  # convert date to year-month --> to be able to match with monthly climate data
  data_sp$date <- as.Date(
    paste(data_sp$year,
          round(data_sp$week*7.6),sep="-"), "%Y-%j")
  data_sp$yr_month <- format(as.Date(data_sp$date), "%Y-%m")

  # phenological data to monthly resolution
  data_sp <- data_sp %>%
    group_by(yr_month) %>%
    dplyr::summarise(mean_value = mean(value),
                     sum_value = sum(value),
                     nr_indiv = length(value))

  # remove if data = NA
  data_sp <- data_sp[!(is.na(data_sp$yr_month)),]


  # merge with climate data
  data_sp <- merge(data_sp, climate, by = "yr_month", all.x = TRUE)

  # cross-correlation + confidence interval
  if(max(data_sp$mean_value)>0){
    corr.anom <- ccf(data_sp$diff, data_sp$mean_value,lag = 18, pl = FALSE)
    # plot(corr.anom, main = paste(species_list[i], " - turnover x anom"))
    ci_value <- qnorm((1 + ci)/2)/sqrt(corr.anom$n.used)
  } else {
    corr.anom$acf <- "NA"
    ci_value <- "NA"
  }

  output_turnover_anom[i,1] <- species_list[i]
  output_turnover_anom[i,2:38] <- as.numeric(corr.anom$acf)
  output_turnover_anom[i,39] <- as.numeric(ci_value)
}

# maximum correlation before time 0 of phenophase
colnames(output_turnover_anom)[40] <- "max.corr.before"
absmax <- function(x) { x[which.max( abs(x) )][1]}
output_turnover_anom$max.corr.before <- apply(output_turnover_anom[,2:19], 1, absmax) # if lag-time changes, this should be adapted --> this is 6 lags max
# only keep significant correlations
output_turnover_anom$max.corr.before <- ifelse(abs(output_turnover_anom$max.corr.before) > output_turnover_anom$ci, output_turnover_anom$max.corr.before, "NA")

# get the timing of the correlations
colnames(output_turnover_anom)[41] <- "max.corr.before.timing"
max.corr.before.timing <- apply(output_turnover_anom[,2:19], 1, function(x) { x[which.max( abs(x) )]}) # if lag-time changes, this should be adapted --> this is 6 lags max
max.corr.before.timing <- names(unlist(max.corr.before.timing))
# remove NA's, because this messes up the colnames of the timing
output_turnover_anom <- output_turnover_anom[complete.cases(output_turnover_anom), ]
output_turnover_anom$max.corr.before.timing <- max.corr.before.timing
# only keep timing of significant correlations
output_turnover_anom$max.corr.before.timing <- ifelse(output_turnover_anom$max.corr.before %in% "NA", "NA", output_turnover_anom$max.corr.before.timing)

turnover_anom <- output_turnover_anom[,(names(output_turnover_anom) %in% c("species_full",
                                                                           "max.corr.before",
                                                                           "max.corr.before.timing"))]
colnames(turnover_anom) <- c("species_full",
                             "ccf_turnover_anom_lead",
                             "ccf_turnover_anom_lead_timing")

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------- Merge  and write table, will be loaded in summary_table.R ------------------
#----------------------------------------------------------------------------------------

output <- merge(overview[1], dormancy_precip, by = "species_full", all.x = TRUE)
output <- merge(output, dormancy_anom, by = "species_full", all.x = TRUE)
output <- merge(output, turnover_precip, by = "species_full", all.x = TRUE)
output <- merge(output, turnover_anom, by = "species_full", all.x = TRUE)

write.table(output, "data/timeseries_correlations.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")
