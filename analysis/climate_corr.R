#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
# graphics.off()
#----- load required packages --------------------------------------------------#
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(tseries)
library(grid)

#----------------------------------------------------------------------
#-------- input that you can change   ---------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------
drought_mm = 150
#----------------------------------------------------------------------



# read climate data
climate <- read.csv("data/yangambi_km5_monthly_kasongo.csv")
climate <- climate %>%
  filter(year > 1936,
         year < 1959) %>%
  # mutate(col = ifelse(precip < drought_mm,"dry","wet"))
  mutate(col = ifelse(month %in% c(1,2,6,7,12), "dry", "wet"))
climate$date <- paste(climate$year,climate$month,"15",sep = "-")
climate$date <- as.Date(climate$date, "%Y-%m-%d")
climate$yr_month <- format(as.Date(climate$date), "%Y-%m")

climate.avg <- read.csv("data/ClimData_monthly_avg.csv",header = TRUE,sep = ",")
climate.avg = climate.avg[,(names(climate.avg) %in% c("Month","insol_JR","tmax_JR"))]
colnames(climate.avg)[1] <- "month"
climate <- merge(climate, climate.avg, by = "month", all.x = TRUE)

# # # check with yearly summaries of climate data
# climate_check <- read.csv("data/yangambi_km5_precip_yearly_kasongo.csv")
# climate_check <- climate_check %>%
#   filter(year > 1936,
#          year < 1960)
# test <- as.data.frame(tapply(climate$precip,list(climate$year),sum))
# colnames(test)[1] <- "precip_sum"
# test$year <- rownames(test)
# test <- merge(climate_check,test,by = "year")
# test$diff <- test$precip_sum - test$precip
# summary(test$diff)
#
# # # corr temp and sun
# climate_check <- read.csv("data/yangambi_km5_monthly_kasongo.csv")
# climate_check <- climate_check %>%
#   filter(year > 1950,
#          year < 2012)
# cor.test(climate_check$precip, climate_check$sun, method = "pearson")
# cor.test(climate_check$tmax, climate_check$sun, method = "pearson")
# cor.test(climate_check$precip, climate_check$sun, method = "pearson")

# remove years where data is inconsistent with the yearly summary
# climate$precip <- ifelse(climate$year %in% c("1940","1943","1949"), NA, climate$precip)


# stationarity and autocorr precip
adf.test(climate$precip) # p value lower then 0.05 = stationary
kpss.test(climate$precip) # p value higher then 0.05 = stationary
autocorr <- acf(climate$precip)


# p_climate <- ggplot() +
#   geom_col(data = climate, aes(x = date,
#                                y = precip,
#                                colour = col,
#                                fill = col)) +
#   scale_colour_manual(values = c("lightcoral","lightblue"),
#                       aesthetics = c("fill","colour")) +
#   labs(y = "precip. (mm/month)",
#        x = "Year") +
#   scale_x_date(date_breaks = "1 year",
#                date_labels = "%Y",
#                limits = as.Date(c('1935-01-01','1960-01-01'))) +
#   theme_minimal() +
#   theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
#         panel.grid.minor.x =  element_blank(),
#         panel.grid.minor.y = element_blank(),
#         panel.background = element_blank(),
#         plot.background = element_rect(fill = 'white', colour = 'white'),
#         strip.text = element_text(hjust = 0),
#         axis.line.x = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1),
#         legend.position="none"
#   )




#----------------------------------------------------------------------
#--------   Phenology data - species correction Meise   ---------------
#----------------------------------------------------------------------
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
# df <- df[which(df$value != 0),]
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
#----------------------------------------------------------------------

# #----------------------------------------------------------------------
# #-------- get the species list ----------------------------------------
# #----------------------------------------------------------------------
overview <- read.csv("data/species_meta_data.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)
species_list <- overview$species_full

# merge number of events to data, so you can filter on minimum number of events
overview <- overview[,(names(overview) %in% c("species_full",
                                             "total_nr_events_leaf_dormancy",
                                             "total_nr_events_leaf_turnover"))]
data <- merge(data, overview, by = c("species_full"), all.x = TRUE)




#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------- DORMANCY x actual precipitation --------------------------------------------
#----------------------------------------------------------------------------------------

output_dormancy_precip <- data.frame(matrix(0,nrow=length(species_list),ncol=17))
colnames(output_dormancy_precip)[1] <- "species_full"
colnames(output_dormancy_precip)[2:14] <- c("t-6","t-5","t-4","t-3","t-2","t-1","t","t+1","t+2","t+3","t+4","t+5","t+6")
colnames(output_dormancy_precip)[15] <- "ci"

ci <- 0.95

for(i in 1:length(species_list)){
  data_sp <- data %>%
    filter(species_full %in% species_list[i],
           phenophase == "leaf_dormancy",
           total_nr_events_leaf_dormancy >= 5)
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
    corr.precip$acf <- NA
    ci_value <- NA
  }
  output_dormancy_precip[i,1] <- species_list[i]
  output_dormancy_precip[i,2:14] <- as.numeric(corr.precip$acf)
  output_dormancy_precip[i,15] <- as.numeric(ci_value)
}

colnames(output_dormancy_precip)[16] <- "corr"
colnames(output_dormancy_precip)[17] <- "corr.timing"
output_dormancy_precip$corr <- ifelse(abs(output_dormancy_precip[,8]) > output_dormancy_precip$ci, output_dormancy_precip[,8], # 0
                               ifelse(abs(output_dormancy_precip[,7]) > output_dormancy_precip$ci, output_dormancy_precip[,7], # t-1
                                      ifelse(abs(output_dormancy_precip[,6]) > output_dormancy_precip$ci, output_dormancy_precip[,6], # t-2
                                             ifelse(abs(output_dormancy_precip[,5]) > output_dormancy_precip$ci, output_dormancy_precip[,5], # t-3
                                                    NA))))
output_dormancy_precip$corr.timing <- ifelse(abs(output_dormancy_precip[,8]) > output_dormancy_precip$ci, "0", # 0
                                      ifelse(abs(output_dormancy_precip[,7]) > output_dormancy_precip$ci, "-1", # t-1
                                             ifelse(abs(output_dormancy_precip[,6]) > output_dormancy_precip$ci, "-2", # t-2
                                                    ifelse(abs(output_dormancy_precip[,5]) > output_dormancy_precip$ci, "-3", # t-3
                                                           NA))))



dormancy_precip <- output_dormancy_precip[,(names(output_dormancy_precip) %in% c("species_full",
                                                                   "corr",
                                                                   "corr.timing"))]
colnames(dormancy_precip) <- c("species_full",
                               "corr_dormancy_precip",
                               "corr_dormancy_precip_timing")

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------- DORMANCY x sun hours --------------------------------------------
#----------------------------------------------------------------------------------------

output_dormancy_sun <- data.frame(matrix(0,nrow=length(species_list),ncol=17))
colnames(output_dormancy_sun)[1] <- "species_full"
colnames(output_dormancy_sun)[2:14] <- c("t-6","t-5","t-4","t-3","t-2","t-1","t","t+1","t+2","t+3","t+4","t+5","t+6")
colnames(output_dormancy_sun)[15] <- "ci"

ci <- 0.95

for(i in 1:length(species_list)){
  data_sp <- data %>%
    filter(species_full %in% species_list[i],
           phenophase == "leaf_dormancy",
           total_nr_events_leaf_dormancy >= 5)
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
  # data_sp <- data_sp %>%
  #   filter(year >= 1951)

  # cross-correlation + confidence interval
  if(max(data_sp$mean_value)>0){
    corr.insol_JR <- ccf(data_sp$insol_JR, data_sp$mean_value,lag = 6, pl = FALSE)
    # plot(corr.insol_JR, main = paste(species_list[i], " - dormancy x insol_JR"))
    ci_value <- qnorm((1 + ci)/2)/sqrt(corr.insol_JR$n.used)
  } else {
    corr.insol_JR$acf <- NA
    ci_value <- NA
  }

  output_dormancy_sun[i,1] <- species_list[i]
  output_dormancy_sun[i,2:14] <- as.numeric(corr.insol_JR$acf)
  output_dormancy_sun[i,15] <- as.numeric(ci_value)
}

# direct correlation at time t
colnames(output_dormancy_sun)[16] <- "corr"
colnames(output_dormancy_sun)[17] <- "corr.timing"
output_dormancy_sun$corr <- ifelse(abs(output_dormancy_sun[,8]) > output_dormancy_sun$ci, output_dormancy_sun[,8], # 0
                               ifelse(abs(output_dormancy_sun[,7]) > output_dormancy_sun$ci, output_dormancy_sun[,7], # t-1
                                      ifelse(abs(output_dormancy_sun[,6]) > output_dormancy_sun$ci, output_dormancy_sun[,6], # t-2
                                             ifelse(abs(output_dormancy_sun[,5]) > output_dormancy_sun$ci, output_dormancy_sun[,5], # t-3
                                                    ifelse(abs(output_dormancy_sun[,4]) > output_dormancy_sun$ci, output_dormancy_sun[,4], # t-4
                                                           ifelse(abs(output_dormancy_sun[,3]) > output_dormancy_sun$ci, output_dormancy_sun[,3], # t-5
                                                                  NA))))))
output_dormancy_sun$corr.timing <- ifelse(abs(output_dormancy_sun[,8]) > output_dormancy_sun$ci, "0", # 0
                                      ifelse(abs(output_dormancy_sun[,7]) > output_dormancy_sun$ci, "-1", # t-1
                                             ifelse(abs(output_dormancy_sun[,6]) > output_dormancy_sun$ci, "-2", # t-2
                                                    ifelse(abs(output_dormancy_sun[,5]) > output_dormancy_sun$ci, "-3", # t-3
                                                           ifelse(abs(output_dormancy_sun[,4]) > output_dormancy_sun$ci, "-4", # t-4
                                                                  ifelse(abs(output_dormancy_sun[,3]) > output_dormancy_sun$ci, "-5", # t-5
                                                                         NA))))))


dormancy_insol_JR <- output_dormancy_sun[,(names(output_dormancy_sun) %in% c("species_full",
                                                                     "corr",
                                                                     "corr.timing"))]
colnames(dormancy_insol_JR) <- c("species_full",
                                 "corr_dormancy_insol_JR",
                                 "corr_dormancy_insol_JR_timing")

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------- DORMANCY x tmax hours --------------------------------------------
#----------------------------------------------------------------------------------------

output_dormancy_tmax <- data.frame(matrix(0,nrow=length(species_list),ncol=17))
colnames(output_dormancy_tmax)[1] <- "species_full"
colnames(output_dormancy_tmax)[2:14] <- c("t-6","t-5","t-4","t-3","t-2","t-1","t","t+1","t+2","t+3","t+4","t+5","t+6")
colnames(output_dormancy_tmax)[15] <- "ci"

ci <- 0.95

for(i in 1:length(species_list)){
  data_sp <- data %>%
    filter(species_full %in% species_list[i],
           phenophase == "leaf_dormancy",
           total_nr_events_leaf_dormancy >= 5)
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
  # data_sp <- data_sp %>%
  #   filter(year >= 1951)

  # cross-correlation + confidence interval
  if(max(data_sp$mean_value)>0){
    corr.tmax_JR <- ccf(data_sp$tmax_JR, data_sp$mean_value,lag = 6, pl = FALSE)
    # plot(corr.tmax_JR, main = paste(species_list[i], " - dormancy x tmax_JR"))
    ci_value <- qnorm((1 + ci)/2)/sqrt(corr.tmax_JR$n.used)
  } else {
    corr.tmax_JR$acf <- NA
    ci_value <- NA
  }
  output_dormancy_tmax[i,1] <- species_list[i]
  output_dormancy_tmax[i,2:14] <- as.numeric(corr.tmax_JR$acf)
  output_dormancy_tmax[i,15] <- as.numeric(ci_value)
}

colnames(output_dormancy_tmax)[16] <- "corr"
colnames(output_dormancy_tmax)[17] <- "corr.timing"
output_dormancy_tmax$corr <- ifelse(abs(output_dormancy_tmax[,8]) > output_dormancy_tmax$ci, output_dormancy_tmax[,8], # 0
                               ifelse(abs(output_dormancy_tmax[,7]) > output_dormancy_tmax$ci, output_dormancy_tmax[,7], # t-1
                                      ifelse(abs(output_dormancy_tmax[,6]) > output_dormancy_tmax$ci, output_dormancy_tmax[,6], # t-2
                                             ifelse(abs(output_dormancy_tmax[,5]) > output_dormancy_tmax$ci, output_dormancy_tmax[,5], # t-3
                                                    ifelse(abs(output_dormancy_tmax[,4]) > output_dormancy_tmax$ci, output_dormancy_tmax[,4], # t-4
                                                           ifelse(abs(output_dormancy_tmax[,3]) > output_dormancy_tmax$ci, output_dormancy_tmax[,3], # t-5
                                                                  NA))))))
output_dormancy_tmax$corr.timing <- ifelse(abs(output_dormancy_tmax[,8]) > output_dormancy_tmax$ci, "0", # 0
                                      ifelse(abs(output_dormancy_tmax[,7]) > output_dormancy_tmax$ci, "-1", # t-1
                                             ifelse(abs(output_dormancy_tmax[,6]) > output_dormancy_tmax$ci, "-2", # t-2
                                                    ifelse(abs(output_dormancy_tmax[,5]) > output_dormancy_tmax$ci, "-3", # t-3
                                                           ifelse(abs(output_dormancy_tmax[,4]) > output_dormancy_tmax$ci, "-4", # t-4
                                                                  ifelse(abs(output_dormancy_tmax[,3]) > output_dormancy_tmax$ci, "-5", # t-5
                                                                         NA))))))


dormancy_tmax_JR <- output_dormancy_tmax[,(names(output_dormancy_tmax) %in% c("species_full",
                                                                    "corr",
                                                                    "corr.timing"))]
colnames(dormancy_tmax_JR) <- c("species_full",
                                "corr_dormancy_tmax_JR",
                                "corr_dormancy_tmax_JR_timing")



#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------- TURNOVER x actual precipitation --------------------------------------------
#----------------------------------------------------------------------------------------

output_turnover_precip <- data.frame(matrix(0,nrow=length(species_list),ncol=17))
colnames(output_turnover_precip)[1] <- "species_full"
colnames(output_turnover_precip)[2:14] <- c("t-6","t-5","t-4","t-3","t-2","t-1","t","t+1","t+2","t+3","t+4","t+5","t+6")
colnames(output_turnover_precip)[15] <- "ci"

ci <- 0.95

for(i in 1:length(species_list)){
  data_sp <- data %>%
    filter(species_full %in% species_list[i],
           phenophase == "leaf_turnover",
           total_nr_events_leaf_turnover >= 5)
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
    corr.precip$acf <- NA
    ci_value <- NA
  }
  output_turnover_precip[i,1] <- species_list[i]
  output_turnover_precip[i,2:14] <- as.numeric(corr.precip$acf)
  output_turnover_precip[i,15] <- as.numeric(ci_value)
}

colnames(output_turnover_precip)[16] <- "corr"
colnames(output_turnover_precip)[17] <- "corr.timing"
output_turnover_precip$corr <- ifelse(abs(output_turnover_precip[,8]) > output_turnover_precip$ci, output_turnover_precip[,8], # 0
                               ifelse(abs(output_turnover_precip[,7]) > output_turnover_precip$ci, output_turnover_precip[,7], # t-1
                                      ifelse(abs(output_turnover_precip[,6]) > output_turnover_precip$ci, output_turnover_precip[,6], # t-2
                                             ifelse(abs(output_turnover_precip[,5]) > output_turnover_precip$ci, output_turnover_precip[,5], # t-3
                                                    NA))))
output_turnover_precip$corr.timing <- ifelse(abs(output_turnover_precip[,8]) > output_turnover_precip$ci, "0", # 0
                                      ifelse(abs(output_turnover_precip[,7]) > output_turnover_precip$ci, "-1", # t-1
                                             ifelse(abs(output_turnover_precip[,6]) > output_turnover_precip$ci, "-2", # t-2
                                                    ifelse(abs(output_turnover_precip[,5]) > output_turnover_precip$ci, "-3", # t-3
                                                           NA))))


turnover_precip <- output_turnover_precip[,(names(output_turnover_precip) %in% c("species_full",
                                                                   "corr",
                                                                   "corr.timing"))]
colnames(turnover_precip) <- c("species_full",
                               "corr_turnover_precip",
                               "corr_turnover_precip_timing")

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------- turnover x sun hours --------------------------------------------
#----------------------------------------------------------------------------------------

output_turnover_sun <- data.frame(matrix(0,nrow=length(species_list),ncol=17))
colnames(output_turnover_sun)[1] <- "species_full"
colnames(output_turnover_sun)[2:14] <- c("t-6","t-5","t-4","t-3","t-2","t-1","t","t+1","t+2","t+3","t+4","t+5","t+6")
colnames(output_turnover_sun)[15] <- "ci"

ci <- 0.95

for(i in 1:length(species_list)){
  data_sp <- data %>%
    filter(species_full %in% species_list[i],
           phenophase == "leaf_turnover",
           total_nr_events_leaf_turnover >= 5)
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
  # data_sp <- data_sp %>%
  #   filter(year >= 1951)

  # cross-correlation + confidence interval
  if(max(data_sp$mean_value)>0){
    corr.insol_JR <- ccf(data_sp$insol_JR, data_sp$mean_value,lag = 6, pl = FALSE)
    # plot(corr.insol_JR, main = paste(species_list[i], " - turnover x insol_JR"))
    ci_value <- qnorm((1 + ci)/2)/sqrt(corr.insol_JR$n.used)
  } else {
    corr.insol_JR$acf <- NA
    ci_value <- NA
  }
  output_turnover_sun[i,1] <- species_list[i]
  output_turnover_sun[i,2:14] <- as.numeric(corr.insol_JR$acf)
  output_turnover_sun[i,15] <- as.numeric(ci_value)
}

colnames(output_turnover_sun)[16] <- "corr"
colnames(output_turnover_sun)[17] <- "corr.timing"
output_turnover_sun$corr <- ifelse(abs(output_turnover_sun[,8]) > output_turnover_sun$ci, output_turnover_sun[,8], # 0
                               ifelse(abs(output_turnover_sun[,7]) > output_turnover_sun$ci, output_turnover_sun[,7], # t-1
                                      ifelse(abs(output_turnover_sun[,6]) > output_turnover_sun$ci, output_turnover_sun[,6], # t-2
                                             ifelse(abs(output_turnover_sun[,5]) > output_turnover_sun$ci, output_turnover_sun[,5], # t-3
                                                    ifelse(abs(output_turnover_sun[,4]) > output_turnover_sun$ci, output_turnover_sun[,4], # t-4
                                                           ifelse(abs(output_turnover_sun[,3]) > output_turnover_sun$ci, output_turnover_sun[,3], # t-5
                                                                  NA))))))
output_turnover_sun$corr.timing <- ifelse(abs(output_turnover_sun[,8]) > output_turnover_sun$ci, "0", # 0
                                      ifelse(abs(output_turnover_sun[,7]) > output_turnover_sun$ci, "-1", # t-1
                                             ifelse(abs(output_turnover_sun[,6]) > output_turnover_sun$ci, "-2", # t-2
                                                    ifelse(abs(output_turnover_sun[,5]) > output_turnover_sun$ci, "-3", # t-3
                                                           ifelse(abs(output_turnover_sun[,4]) > output_turnover_sun$ci, "-4", # t-4
                                                                  ifelse(abs(output_turnover_sun[,3]) > output_turnover_sun$ci, "-5", # t-5
                                                                         NA))))))


turnover_insol_JR <- output_turnover_sun[,(names(output_turnover_sun) %in% c("species_full",
                                                                     "corr",
                                                                     "corr.timing"))]
colnames(turnover_insol_JR) <- c("species_full",
                                 "corr_turnover_insol_JR",
                                 "corr_turnover_insol_JR_timing")
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------- turnover x tmax hours --------------------------------------------
#----------------------------------------------------------------------------------------

output_turnover_tmax <- data.frame(matrix(0,nrow=length(species_list),ncol=17))
colnames(output_turnover_tmax)[1] <- "species_full"
colnames(output_turnover_tmax)[2:14] <- c("t-6","t-5","t-4","t-3","t-2","t-1","t","t+1","t+2","t+3","t+4","t+5","t+6")
colnames(output_turnover_tmax)[15] <- "ci"

ci <- 0.95

for(i in 1:length(species_list)){
  data_sp <- data %>%
    filter(species_full %in% species_list[i],
           phenophase == "leaf_turnover",
           total_nr_events_leaf_turnover >= 5)
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
  # data_sp <- data_sp %>%
  #   filter(year >= 1951)

  # cross-correlation + confidence interval
  if(max(data_sp$mean_value)>0){
    corr.tmax_JR <- ccf(data_sp$tmax_JR, data_sp$mean_value,lag = 6, pl = FALSE)
    # plot(corr.tmax_JR, main = paste(species_list[i], " - turnover x tmax_JR"))
    ci_value <- qnorm((1 + ci)/2)/sqrt(corr.tmax_JR$n.used)
  } else {
    corr.tmax_JR$acf <- NA
    ci_value <- NA
  }
  output_turnover_tmax[i,1] <- species_list[i]
  output_turnover_tmax[i,2:14] <- as.numeric(corr.tmax_JR$acf)
  output_turnover_tmax[i,15] <- as.numeric(ci_value)
}

colnames(output_turnover_tmax)[16] <- "corr"
colnames(output_turnover_tmax)[17] <- "corr.timing"
output_turnover_tmax$corr <- ifelse(abs(output_turnover_tmax[,8]) > output_turnover_tmax$ci, output_turnover_tmax[,8], # 0
                               ifelse(abs(output_turnover_tmax[,7]) > output_turnover_tmax$ci, output_turnover_tmax[,7], # t-1
                                      ifelse(abs(output_turnover_tmax[,6]) > output_turnover_tmax$ci, output_turnover_tmax[,6], # t-2
                                             ifelse(abs(output_turnover_tmax[,5]) > output_turnover_tmax$ci, output_turnover_tmax[,5], # t-3
                                                    ifelse(abs(output_turnover_tmax[,4]) > output_turnover_tmax$ci, output_turnover_tmax[,4], # t-4
                                                           ifelse(abs(output_turnover_tmax[,3]) > output_turnover_tmax$ci, output_turnover_tmax[,3], # t-5
                                                                  NA))))))
output_turnover_tmax$corr.timing <- ifelse(abs(output_turnover_tmax[,8]) > output_turnover_tmax$ci, "0", # 0
                                      ifelse(abs(output_turnover_tmax[,7]) > output_turnover_tmax$ci, "-1", # t-1
                                             ifelse(abs(output_turnover_tmax[,6]) > output_turnover_tmax$ci, "-2", # t-2
                                                    ifelse(abs(output_turnover_tmax[,5]) > output_turnover_tmax$ci, "-3", # t-3
                                                           ifelse(abs(output_turnover_tmax[,4]) > output_turnover_tmax$ci, "-4", # t-4
                                                                  ifelse(abs(output_turnover_tmax[,3]) > output_turnover_tmax$ci, "-5", # t-5
                                                                         NA))))))


turnover_tmax_JR <- output_turnover_tmax[,(names(output_turnover_tmax) %in% c("species_full",
                                                                    "corr",
                                                                    "corr.timing"))]
colnames(turnover_tmax_JR) <- c("species_full",
                                "corr_turnover_tmax_JR",
                                "corr_turnover_tmax_JR_timing")

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------- Merge  and write table, will be loaded in summary_table.R ------------------
#----------------------------------------------------------------------------------------
output <- merge(overview[1], dormancy_precip, by = "species_full", all.x = TRUE)
output <- merge(output, turnover_precip, by = "species_full", all.x = TRUE)
output <- merge(output, dormancy_insol_JR, by = "species_full", all.x = TRUE)
output <- merge(output, turnover_insol_JR, by = "species_full", all.x = TRUE)
output <- merge(output, dormancy_tmax_JR, by = "species_full", all.x = TRUE)
output <- merge(output, turnover_tmax_JR, by = "species_full", all.x = TRUE)

write.table(output, "data/timeseries_correlations.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# get an idea of how many species show correlations, at what average timing etc...
# in groups of deciduous and evergreen
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------


overview <- read.csv("data/species_meta_data.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)
# clear
overview$deciduousness <- ifelse(overview$species_full %in% "Pericopsis elata", "deciduous*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Trichilia welwitschii", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Copaifera mildbraedii", "deciduous*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Tridesmostemon omphalocarpoides", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Omphalocarpum lecomteanum", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Fernandoa adolfi-friderici", "deciduous*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Tabernaemontana crassa", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Trichilia tessmannii", "deciduous*",overview$deciduousness)
# not so sure, limited data
overview$deciduousness <- ifelse(overview$species_full %in% "Trichilia gilletii", "evergreen* (?)",overview$deciduousness)
# not so sure, unclear phenological data
overview$deciduousness <- ifelse(overview$species_full %in% "Radlkofera calodendron", "evergreen* (?)",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Gilletiodendron mildbraedii", "evergreen* (?)",overview$deciduousness)

## two stars, in literature found as evergreen or (sometimes) deciduous
## selected a class based on the actual data
overview$deciduousness <- ifelse(overview$species_full %in% "Celtis mildbraedii", "evergreen**",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Combretum lokele", "deciduous**",overview$deciduousness)

overview$deciduousness <- ifelse(overview$species_full %in% "Homalium africanum", "evergreen**",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Quassia silvestris", "evergreen**",overview$deciduousness)

# not so sure, unclear phenological data
overview$deciduousness <- ifelse(overview$species_full %in% "Homalium longistylum", "evergreen** (?)",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Irvingia gabonensis", "deciduous** (?)",overview$deciduousness)



# only keep the columns you want to work with
overview = overview[,(names(overview) %in% c("species_full",
                                             "deciduousness",
                                             "site_years_with_leaf_dormancy",
                                             "site_years_with_leaf_turnover",
                                             "total_nr_events_leaf_dormancy",
                                             "total_nr_events_leaf_turnover"
))]

df <- merge(overview, output, by = "species_full", all.x = TRUE)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# corr summary figure
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# dormancy
df_dorm <- df %>%
  # filter(grepl("deciduous",deciduousness)) %>%
  filter(site_years_with_leaf_dormancy > 0) #%>%
  # filter(total_nr_events_leaf_dormancy >=5)

df_dorm$precip_phase <- ifelse(df_dorm$total_nr_events_leaf_dormancy <5, "few-events",
                               ifelse(df_dorm$corr_dormancy_precip_timing == "0" & df_dorm$corr_dormancy_precip < 0, "no-lag-neg",
                                   ifelse(df_dorm$corr_dormancy_precip_timing == "0" & df_dorm$corr_dormancy_precip > 0, "no-lag-pos",
                                          ifelse(df_dorm$corr_dormancy_precip_timing %in% c("-1","-2","-3") & df_dorm$corr_dormancy_precip < 0, "lag-neg",
                                                 ifelse(df_dorm$corr_dormancy_precip_timing %in% c("-1","-2","-3") & df_dorm$corr_dormancy_precip > 0, "lag-pos",
                                                        NA)))))
df_dorm$precip_phase <- ifelse(is.na(df_dorm$precip_phase) & df_dorm$total_nr_events_leaf_dormancy >= 5, "h-no-corr", df_dorm$precip_phase)


df_dorm$insol_phase <- ifelse(df_dorm$total_nr_events_leaf_dormancy <5, "few-events",
                              ifelse(df_dorm$corr_dormancy_insol_JR_timing == "0" & df_dorm$corr_dormancy_insol_JR < 0, "no-lag-neg",
                                   ifelse(df_dorm$corr_dormancy_insol_JR_timing == "0" & df_dorm$corr_dormancy_insol_JR > 0, "no-lag-pos",
                                          ifelse(df_dorm$corr_dormancy_insol_JR_timing %in% c("-1","-2","-3") & df_dorm$corr_dormancy_insol_JR < 0, "lag-neg",
                                                 ifelse(df_dorm$corr_dormancy_insol_JR_timing %in% c("-1","-2","-3") & df_dorm$corr_dormancy_insol_JR > 0, "lag-pos",
                                                        NA)))))
df_dorm$insol_phase <- ifelse(is.na(df_dorm$insol_phase) & df_dorm$total_nr_events_leaf_dormancy >= 5, "h-no-corr", df_dorm$insol_phase)

df_dorm$tmax_phase <- ifelse(df_dorm$total_nr_events_leaf_dormancy <5, "few-events",
                             ifelse(df_dorm$corr_dormancy_tmax_JR_timing == "0" & df_dorm$corr_dormancy_tmax_JR < 0, "no-lag-neg",
                                  ifelse(df_dorm$corr_dormancy_tmax_JR_timing == "0" & df_dorm$corr_dormancy_tmax_JR > 0, "no-lag-pos",
                                         ifelse(df_dorm$corr_dormancy_tmax_JR_timing %in% c("-1","-2","-3") & df_dorm$corr_dormancy_tmax_JR < 0, "lag-neg",
                                                ifelse(df_dorm$corr_dormancy_tmax_JR_timing %in% c("-1","-2","-3") & df_dorm$corr_dormancy_tmax_JR > 0, "lag-pos",
                                                       NA)))))
df_dorm$tmax_phase <- ifelse(is.na(df_dorm$tmax_phase) & df_dorm$total_nr_events_leaf_dormancy >= 5, "h-no-corr", df_dorm$tmax_phase)


# turnover
df_turn <- df %>%
  # filter(grepl("deciduous",deciduousness)) %>%
  filter(site_years_with_leaf_turnover > 0) #%>%
  # filter(total_nr_events_leaf_turnover >=5)

df_turn$precip_phase <- ifelse(df_turn$total_nr_events_leaf_turnover <5, "few-events",
                               ifelse(df_turn$corr_turnover_precip_timing == "0" & df_turn$corr_turnover_precip < 0, "no-lag-neg",
                               ifelse(df_turn$corr_turnover_precip_timing == "0" & df_turn$corr_turnover_precip > 0, "no-lag-pos",
                                      ifelse(df_turn$corr_turnover_precip_timing %in% c("-1","-2","-3") & df_turn$corr_turnover_precip < 0, "lag-neg",
                                             ifelse(df_turn$corr_turnover_precip_timing %in% c("-1","-2","-3") & df_turn$corr_turnover_precip > 0, "lag-pos",
                                                    NA)))))
df_turn$precip_phase <- ifelse(is.na(df_turn$precip_phase) & df_turn$total_nr_events_leaf_turnover >= 5, "h-no-corr", df_turn$precip_phase)


df_turn$insol_phase <- ifelse(df_turn$total_nr_events_leaf_turnover <5, "few-events",
                              ifelse(df_turn$corr_turnover_insol_JR_timing == "0" & df_turn$corr_turnover_insol_JR < 0, "no-lag-neg",
                              ifelse(df_turn$corr_turnover_insol_JR_timing == "0" & df_turn$corr_turnover_insol_JR > 0, "no-lag-pos",
                                     ifelse(df_turn$corr_turnover_insol_JR_timing %in% c("-1","-2","-3") & df_turn$corr_turnover_insol_JR < 0, "lag-neg",
                                            ifelse(df_turn$corr_turnover_insol_JR_timing %in% c("-1","-2","-3") & df_turn$corr_turnover_insol_JR > 0, "lag-pos",
                                                   NA)))))
df_turn$insol_phase <- ifelse(is.na(df_turn$insol_phase) & df_turn$total_nr_events_leaf_turnover >= 5, "h-no-corr", df_turn$insol_phase)

df_turn$tmax_phase <- ifelse(df_turn$total_nr_events_leaf_turnover <5, "few-events",
                             ifelse(df_turn$corr_turnover_tmax_JR_timing == "0" & df_turn$corr_turnover_tmax_JR < 0, "no-lag-neg",
                             ifelse(df_turn$corr_turnover_tmax_JR_timing == "0" & df_turn$corr_turnover_tmax_JR > 0, "no-lag-pos",
                                    ifelse(df_turn$corr_turnover_tmax_JR_timing %in% c("-1","-2","-3") & df_turn$corr_turnover_tmax_JR < 0, "lag-neg",
                                           ifelse(df_turn$corr_turnover_tmax_JR_timing %in% c("-1","-2","-3") & df_turn$corr_turnover_tmax_JR > 0, "lag-pos",
                                                  NA)))))
df_turn$tmax_phase <- ifelse(is.na(df_turn$tmax_phase) & df_turn$total_nr_events_leaf_turnover >= 5, "h-no-corr", df_turn$tmax_phase)


#-----------------------------
# evergreen - dormancy
#-----------------------------
df_ever_dorm <- df_dorm %>%
  filter(grepl("evergreen",deciduousness))

ed_precip <- as.data.frame(tapply(df_ever_dorm$corr_dormancy_precip, list(df_ever_dorm$precip_phase), length)) #ed = ever dorm
counts_ever_dorm <- length(df_ever_dorm$species_full)
colnames(ed_precip) <- "value"
ed_precip$value <- ed_precip$value / counts_ever_dorm *100
ed_precip$relation <- rownames(ed_precip)
ed_precip$variable <- "precipitation"
# # no corr with precip -> make empty dataframe
# ed_precip <- data.frame(
#   value = c(0,0,0,0), #,11/2/94,11/2/94),
#   relation = c("no-lag-neg","no-lag-pos","lag-neg","lag-pos"), #,"enough-events-pos","enough-events-neg"),
#   variable = "precipitation")

ed_insol <- as.data.frame(tapply(df_ever_dorm$corr_dormancy_insol_JR, list(df_ever_dorm$insol_phase), length)) #ed = ever dorm
counts_ever_dorm <- length(df_ever_dorm$species_full)
colnames(ed_insol) <- "value"
ed_insol$value <- ed_insol$value / counts_ever_dorm *100
ed_insol$relation <- rownames(ed_insol)
ed_insol$variable <- "sunhours"

ed_tmax <- as.data.frame(tapply(df_ever_dorm$corr_dormancy_tmax_JR, list(df_ever_dorm$tmax_phase), length)) #ed = ever dorm
counts_ever_dorm <- length(df_ever_dorm$species_full)
colnames(ed_tmax) <- "value"
ed_tmax$value <- ed_tmax$value / counts_ever_dorm *100
ed_tmax$relation <- rownames(ed_tmax)
ed_tmax$variable <- "tmax"

ed_summary <- rbind(ed_precip, ed_insol, ed_tmax)

ed_summary$value <- ifelse(ed_summary$relation %in% c("few-events","h-no-corr"), ed_summary$value / 2, ed_summary$value)
hack <- ed_summary %>%
  filter(relation %in% "few-events")
hack$relation <- "few-events-neg"
hack2 <- ed_summary %>%
  filter(relation %in% "h-no-corr")
hack2$relation <- "h-no-corr-neg"
ed_summary <- rbind(ed_summary, hack)
ed_summary <- rbind(ed_summary, hack2)

p_ed <- ggplot(ed_summary,
       aes(x = variable,
           y = ifelse(relation %in% c("no-lag-neg","lag-neg","few-events-neg","h-no-corr-neg"), -value, value), #,"enough-events-neg"
           fill = relation)) +
  geom_col() +
  scale_y_continuous(limits = c(-70,70),
                     breaks = c(-50,-25,0,25,50),
                     labels = c("", "","","",""),
                     sec.axis = dup_axis(name = "test",
                                         breaks = c(-25,25),
                                         labels = c("negative correlations", "positive correlations"))) +
  scale_x_discrete(limits = rev(levels(as.factor(ed_summary$variable)))) +
  coord_flip() +
  scale_fill_manual(values = c("grey90", "grey90","#dfc27d", "#dfc27d", "#80cdc1","#018571","#018571")) + # no label for lag-pos, so 1 #80cdc1 removed
  geom_hline(yintercept = 0, color =c("white")) +
  annotate("text", x = 3.2, y = 65, label = paste("n = ",counts_ever_dorm), size = 3) +
  labs(y = "",
       x = "Dormancy")  +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black",vjust = 5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 3),
        legend.position = "none",
        plot.margin = unit(c(0.5,0,0,0.5),"cm")
  )



#-----------------------------
# evergreen - turnover
#-----------------------------
df_ever_turn <- df_turn %>%
  filter(grepl("evergreen",deciduousness))

et_precip <- as.data.frame(tapply(df_ever_turn$corr_dormancy_precip, list(df_ever_turn$precip_phase), length)) #et = ever dorm
counts_ever_turn <- length(df_ever_turn$species_full)
colnames(et_precip) <- "value"
et_precip$value <- et_precip$value / counts_ever_turn *100
et_precip$relation <- rownames(et_precip)
et_precip$variable <- "precipitation"

et_insol <- as.data.frame(tapply(df_ever_turn$corr_dormancy_insol_JR, list(df_ever_turn$insol_phase), length)) #et = ever dorm
counts_ever_turn <- length(df_ever_turn$species_full)
colnames(et_insol) <- "value"
et_insol$value <- et_insol$value / counts_ever_turn *100
et_insol$relation <- rownames(et_insol)
et_insol$variable <- "sunhours"

et_tmax <- as.data.frame(tapply(df_ever_turn$corr_dormancy_tmax_JR, list(df_ever_turn$tmax_phase), length)) #et = ever dorm
counts_ever_turn <- length(df_ever_turn$species_full)
colnames(et_tmax) <- "value"
et_tmax$value <- et_tmax$value / counts_ever_turn *100
et_tmax$relation <- rownames(et_tmax)
et_tmax$variable <- "tmax"

et_summary <- rbind(et_precip, et_insol, et_tmax)

et_summary$value <- ifelse(et_summary$relation %in% c("few-events","h-no-corr"), et_summary$value / 2, et_summary$value)
hack <- et_summary %>%
  filter(relation %in% "few-events")
hack$relation <- "few-events-neg"
hack2 <- et_summary %>%
  filter(relation %in% "h-no-corr")
hack2$relation <- "h-no-corr-neg"
et_summary <- rbind(et_summary, hack)
et_summary <- rbind(et_summary, hack2)


p_et <- ggplot(et_summary,
       aes(x = variable,
           y = ifelse(relation %in% c("no-lag-neg","lag-neg","few-events-neg","h-no-corr-neg"), -value, value),
           fill = relation)) +
  geom_col() +
  scale_y_continuous(limits = c(-70,70),
                     breaks = c(-50,-25,0,25,50),
                     labels = c("", "","","","")) +
  scale_x_discrete(limits = rev(levels(as.factor(et_summary$variable)))) +
  coord_flip() +
  scale_fill_manual(values = c("grey90", "grey90","#dfc27d", "#dfc27d", "#80cdc1", "#80cdc1", "#018571","#018571")) +
  geom_hline(yintercept = 0, color =c("white")) +
  annotate("text", x = 3.2, y = 65, label = paste("n = ",counts_ever_turn), size = 3) +
  labs(y = "",
       x = "Turnover") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        # axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 3),
        legend.position = "none",
        plot.margin = unit(c(0,0,0.2,0.5),"cm")
  )


#-----------------------------
# deciduous - dormancy
#-----------------------------
df_dec_dorm <- df_dorm %>%
  filter(grepl("deciduous",deciduousness))

dd_precip <- as.data.frame(tapply(df_dec_dorm$corr_dormancy_precip, list(df_dec_dorm$precip_phase), length)) #dd = dec dorm
counts_dec_dorm <- length(df_dec_dorm$species_full)
colnames(dd_precip) <- "value"
dd_precip$value <- dd_precip$value / counts_dec_dorm *100
dd_precip$relation <- rownames(dd_precip)
dd_precip$variable <- "precipitation"

dd_insol <- as.data.frame(tapply(df_dec_dorm$corr_dormancy_insol_JR, list(df_dec_dorm$insol_phase), length)) #dd = dec dorm
counts_dec_dorm <- length(df_dec_dorm$species_full)
colnames(dd_insol) <- "value"
dd_insol$value <- dd_insol$value / counts_dec_dorm *100
dd_insol$relation <- rownames(dd_insol)
dd_insol$variable <- "sunhours"

dd_tmax <- as.data.frame(tapply(df_dec_dorm$corr_dormancy_tmax_JR, list(df_dec_dorm$tmax_phase), length)) #dd = dec dorm
counts_dec_dorm <- length(df_dec_dorm$species_full)
colnames(dd_tmax) <- "value"
dd_tmax$value <- dd_tmax$value / counts_dec_dorm *100
dd_tmax$relation <- rownames(dd_tmax)
dd_tmax$variable <- "tmax"

dd_summary <- rbind(dd_precip, dd_insol, dd_tmax)

dd_summary$value <- ifelse(dd_summary$relation %in% c("few-events","h-no-corr"), dd_summary$value / 2, dd_summary$value)
hack <- dd_summary %>%
  filter(relation %in% "few-events")
hack$relation <- "few-events-neg"
hack2 <- dd_summary %>%
  filter(relation %in% "h-no-corr")
hack2$relation <- "h-no-corr-neg"
dd_summary <- rbind(dd_summary, hack)
dd_summary <- rbind(dd_summary, hack2)


p_dd <- ggplot(dd_summary,
       aes(x = variable,
           y = ifelse(relation %in% c("no-lag-neg","lag-neg","few-events-neg","h-no-corr-neg"), -value, value),
           fill = relation)) +
  geom_col() +
  scale_y_continuous(limits = c(-70,70),
                     breaks = c(-50,-25,0,25,50),
                     labels = c("", "","","","")) +
  scale_x_discrete(limits = rev(levels(as.factor(dd_summary$variable)))) +
  coord_flip() +
  scale_fill_manual(values = c("grey90", "grey90","#dfc27d", "#dfc27d", "#80cdc1", "#80cdc1", "#018571","#018571")) +
  geom_hline(yintercept = 0, color =c("white")) +
  annotate("text", x = 3.2, y = 65, label = paste("n = ",counts_dec_dorm), size = 3) +
  labs(y = "",
       x = "Dormancy") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        # axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 3),
        legend.position = "none",
        plot.margin = unit(c(0,0,0.2,0.5),"cm")
  )

#-----------------------------
# deciduous - turnover
#-----------------------------
df_dec_turn <- df_turn %>%
  filter(grepl("deciduous",deciduousness))

dt_precip <- as.data.frame(tapply(df_dec_turn$corr_dormancy_precip, list(df_dec_turn$precip_phase), length)) #dt = dec turn
counts_dec_turn <- length(df_dec_turn$species_full)
colnames(dt_precip) <- "value"
dt_precip$value <- dt_precip$value / counts_dec_turn *100
dt_precip$relation <- rownames(dt_precip)
dt_precip$variable <- "precipitation"

dt_insol <- as.data.frame(tapply(df_dec_turn$corr_dormancy_insol_JR, list(df_dec_turn$insol_phase), length)) #dt = dec turn
counts_dec_turn <- length(df_dec_turn$species_full)
colnames(dt_insol) <- "value"
dt_insol$value <- dt_insol$value / counts_dec_turn *100
dt_insol$relation <- rownames(dt_insol)
dt_insol$variable <- "sunhours"

dt_tmax <- as.data.frame(tapply(df_dec_turn$corr_dormancy_tmax_JR, list(df_dec_turn$tmax_phase), length)) #dt = dec turn
counts_dec_turn <- length(df_dec_turn$species_full)
colnames(dt_tmax) <- "value"
dt_tmax$value <- dt_tmax$value / counts_dec_turn *100
dt_tmax$relation <- rownames(dt_tmax)
dt_tmax$variable <- "tmax"

dt_summary <- rbind(dt_precip, dt_insol, dt_tmax)

dt_summary$value <- ifelse(dt_summary$relation %in% c("few-events","h-no-corr"), dt_summary$value / 2, dt_summary$value)
hack <- dt_summary %>%
  filter(relation %in% "few-events")
hack$relation <- "few-events-neg"
hack2 <- dt_summary %>%
  filter(relation %in% "h-no-corr")
hack2$relation <- "h-no-corr-neg"
dt_summary <- rbind(dt_summary, hack)
dt_summary <- rbind(dt_summary, hack2)


p_dt <- ggplot(dt_summary,
       aes(x = variable,
           y = ifelse(relation %in% c("no-lag-neg","lag-neg","few-events-neg","h-no-corr-neg"), -value, value),
           fill = relation)) +
  geom_col() +
  scale_y_continuous(limits = c(-70,70),
                     breaks = c(-50,-25,0,25,50),
                     labels = c(50,25,0,25,50)) +
  scale_x_discrete(limits = rev(levels(as.factor(dt_summary$variable)))) +
  coord_flip() +
  scale_fill_manual(values = c("grey90", "grey90","#dfc27d", "#dfc27d", "#80cdc1", "#80cdc1", "#018571","#018571"),
                    breaks = c("no-lag-neg","lag-neg","h-no-corr-neg","few-events-neg"),
                    labels = c(" in-phase   "," lag   "," no correlation   "," too few events   ")) +
  geom_hline(yintercept = 0, color =c("white")) +
  annotate("text", x = 3.2, y = 65, label = paste("n = ",counts_dec_turn), size = 3) +
  labs(y = "% species with neg. or pos. correlations",
       x = "Turnover") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_text(),
        axis.title.x = element_text(),
        axis.title.y = element_text(vjust = 3),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.margin = unit(c(0,0,0.2,0.5),"cm")
  )


#-----------------------------
# summary plot
#-----------------------------
# p_all <- grid.arrange(p_ed, p_et, p_dd, p_dt, heights = c(1,1,1,1))

# p_ed <- ggplot_gtable(ggplot_build(p_ed))
# p_et <- ggplot_gtable(ggplot_build(p_et))
# p_dd <- ggplot_gtable(ggplot_build(p_dd))
# p_dt <- ggplot_gtable(ggplot_build(p_dt))
#
# p_et$heights <-p_ed$heights
# p_dd$heights <-p_ed$heights
# p_dt$heights <-p_ed$heights


p_all <- grid.arrange(arrangeGrob(p_ed, p_et, heights = c(1,0.65),
                                 left = textGrob("Evergreen", gp=gpar(fontsize=12), rot = 90, hjust = 0.7)), #, hjust = 0.05, vjust = 2
                     arrangeGrob(p_dd, p_dt, heights = c(0.5,1),
                                 left = textGrob("Deciduous", gp=gpar(fontsize=12), rot = 90, hjust = 0.02)),
                     ncol = 1,
                     heights = c(0.85,1))

pdf("~/Desktop/figure3_corr.pdf",6,4) # 5,10)
plot(p_all)
dev.off()



p_all <- grid.arrange(arrangeGrob(p_ed, p_et, heights = c(1,1),
                                  left = grobTree(rectGrob(gp=gpar(fill="grey", border = "white"),hjust = 0.05),
                                                  textGrob("Evergreen", gp=gpar(fontsize=12), rot = 90))), #, hjust = 0.05, vjust = 2
                      arrangeGrob(p_dd, p_dt, heights = c(1,1),
                                  left = textGrob("Deciduous", gp=gpar(fontsize=12), rot = 90)),
                      ncol = 1,
                      heights = c(1,1))
