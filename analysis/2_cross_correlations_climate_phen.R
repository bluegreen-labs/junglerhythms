#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
#----- load required libraries -------------------------------------------------#
library(tidyverse)
library(plyr)
library(ggplot2)
library(gtools)
library(grid)
library(gridExtra)
library(tseries)
library(showtext)
font_add_google(
  "Lato",
  regular.wt = 300,
  bold.wt = 700)
#----- source required functions -----------------------------------------------#
source("R/timeline_gap_fill.R")
source("R/consec_years.R")
source("R/climate_ccf_function.R")
#-------------------------------------------------------------------------------#


#----------------------------------------------------------------------
#--------   Phenology data --------------------------------------------
#----------------------------------------------------------------------
data <- readRDS("data/jungle_rhythms_data_cleaned.rds")
#----------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------------
#--- Climate data -----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
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


#----------------------------------------------------------------------------------------------------------------------
#--- get selected species and clean time series -----------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
overview <- read.csv("data/species_meta_data_phase2.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)
# get extra info on species from 'overview' file
# in summary figure, only species with events are included
overview <- overview %>%
  dplyr::select(species_full,
                deciduousness,
                total_nr_events_leaf_dormancy, # species with too few events (< 5) not included for cross correlation analysis
                total_nr_events_leaf_turnover)

species_list <- overview$species_full

#--- leaf turnover ----------------------------------------------------------------------------------------------------
# for selected species and phenophase: get extended timelines at ID level with 2 year-gaps filled with zero
timelines_id_turn <- missing_year_gaps(data = data,
                                       species_name = species_list,
                                       pheno = "leaf_turnover",
                                       gapfill_missingyears = 0)
# # for fourier, no gaps in timelines allowed
# # get longest consecutive timeline; at species level
# timelines_sp_consec_turn <- consecutive_timeline_sp(data = timelines_id_turn,
#                                                species_name = species_list,
#                                                pheno = "leaf_turnover")
# timelines_sp_consec_turn <- merge(overview, timelines_sp_consec_turn, by = "species_full", all.x = TRUE)

timelines_sp_turn <- timelines_id_turn %>%
  group_by(species_full, date) %>%
  dplyr::summarise(mean_value = mean(value, na.rm = TRUE))
timelines_sp_turn$mean_value[is.nan(timelines_sp_turn$mean_value)] <- NA
timelines_sp_turn$scaled_value <- ifelse(timelines_sp_turn$mean_value > 0, 1,
                                         ifelse(timelines_sp_turn$mean_value == 0, 0, NA))
timelines_sp_turn <- merge(overview, timelines_sp_turn, by = "species_full", all.x = TRUE)
timelines_sp_turn$phenophase <- "leaf_turnover"
#--- leaf dormancy ----------------------------------------------------------------------------------------------------
timelines_id_dorm <- missing_year_gaps(data = data,
                                       species_name = species_list,
                                       pheno = "leaf_dormancy",
                                       gapfill_missingyears = 0)
# timelines_sp_consec_dorm <- consecutive_timeline_sp(data = timelines_id_dorm,
#                                                     species_name = species_list,
#                                                     pheno = "leaf_dormancy")
# timelines_sp_consec_dorm <- merge(overview, timelines_sp_consec_dorm, by = "species_full", all.x = TRUE)
timelines_sp_dorm <- timelines_id_dorm %>%
  group_by(species_full, date) %>%
  dplyr::summarise(mean_value = mean(value, na.rm = TRUE))
timelines_sp_dorm$mean_value[is.nan(timelines_sp_dorm$mean_value)] <- NA
timelines_sp_dorm$scaled_value <- ifelse(timelines_sp_dorm$mean_value > 0, 1,
                                         ifelse(timelines_sp_dorm$mean_value == 0, 0, NA))
timelines_sp_dorm <- merge(overview, timelines_sp_dorm, by = "species_full", all.x = TRUE)
timelines_sp_dorm$phenophase <- "leaf_dormancy"
#----------------------------------------------------------------------------------------------------------------------

# event_count <- function(
#   data = data,
#   species_name = "Afzelia bipindensis"){
#   el_output <- data.frame()
#
#   for (j in 1:length(species_name)){
#     data_subset <- data %>%
#       filter(species_full %in% species_name[j])
#     # sort dataframe according to date
#     data_subset <- data_subset %>%
#       dplyr::arrange(date)
#     # get first differences
#     diff_values <- diff(data_subset$scaled_value)
#     # get matching info
#     start <- data_subset[which(diff_values == 1) + 1,]
#     a <- length(start$scaled_value)
#     el_sp <- data.frame(species_full = species_name[j],
#                         nr_events = length(start$scaled_value))
#     el_output <- rbind(el_output,el_sp)
#   }
#   return(el_output)
# }
#
# dorm_events <- event_count(data = timelines_sp_consec_dorm,
#                          species_name = unique(timelines_sp_consec_dorm$species_full))
# dorm_events <- dorm_events %>%
#   dplyr::rename(nr_dorm_events_consec = nr_events)
#
# turn_events <- event_count(data = timelines_sp_consec_turn,
#                               species_name = unique(timelines_sp_consec_turn$species_full))
# turn_events <- turn_events %>%
#   dplyr::rename(nr_turn_events_consec = nr_events)

#----------------------------------------------------------------------------------------------------------------------
#--- cross correlations climate - phenology timeseries ----------------------------------------------------------------
#--- phenology timeseries with too few events (< 5) are not analysed within the function (set to NA) ------------------
#----------------------------------------------------------------------------------------------------------------------
# leaf turnover
crosscorr <- climate_ccf(data = timelines_sp_turn,
                         climate = climate,
                         climate.variable = "precip",
                         species_name = species_list,
                         pheno = "leaf_turnover")
crosscorr2 <- climate_ccf(data = timelines_sp_turn,
                          climate = climate,
                          climate.variable = "sun",
                          species_name = species_list,
                          pheno = "leaf_turnover")
crosscorr3 <- climate_ccf(data = timelines_sp_turn,
                          climate = climate,
                          climate.variable = "temp",
                          species_name = species_list,
                          pheno = "leaf_turnover")
# leaf dormancy
crosscorr4 <- climate_ccf(data = timelines_sp_dorm,
                          climate = climate,
                          climate.variable = "precip",
                          species_name = species_list,
                          pheno = "leaf_dormancy")
crosscorr5 <- climate_ccf(data = timelines_sp_dorm,
                          climate = climate,
                          climate.variable = "sun",
                          species_name = species_list,
                          pheno = "leaf_dormancy")
crosscorr6 <- climate_ccf(data = timelines_sp_dorm,
                          climate = climate,
                          climate.variable = "temp",
                          species_name = species_list,
                          pheno = "leaf_dormancy")
# merge
# output <- merge(overview[1], dormancy_precip, by = "species_full", all.x = TRUE)
output <- merge(crosscorr, crosscorr2, by = "species_full", all.x = TRUE)
output <- merge(output, crosscorr3, by = "species_full", all.x = TRUE)
output <- merge(output, crosscorr4, by = "species_full", all.x = TRUE)
output <- merge(output, crosscorr5, by = "species_full", all.x = TRUE)
output <- merge(output, crosscorr6, by = "species_full", all.x = TRUE)

# is.na didn;t work in the big ifelse statement
output$corr_leaf_turnover_precip <- ifelse(is.na(output$corr_leaf_turnover_precip), 0, output$corr_leaf_turnover_precip)
output$corr_leaf_turnover_sun <- ifelse(is.na(output$corr_leaf_turnover_sun), 0, output$corr_leaf_turnover_sun)
output$corr_leaf_turnover_temp <- ifelse(is.na(output$corr_leaf_turnover_temp), 0, output$corr_leaf_turnover_temp)


output$turn_groups <- ifelse(output$corr_leaf_turnover_precip < 0 &  output$corr_leaf_turnover_sun > 0 &    output$corr_leaf_turnover_temp > 0, "group1",
                             ifelse(output$corr_leaf_turnover_precip == 0 & output$corr_leaf_turnover_sun > 0 &    output$corr_leaf_turnover_temp > 0, "group1",
                                    ifelse(output$corr_leaf_turnover_precip == 0 & output$corr_leaf_turnover_sun > 0 &    output$corr_leaf_turnover_temp == 0, "group1",
                                           ifelse(output$corr_leaf_turnover_precip == 0 & output$corr_leaf_turnover_sun == 0 &   output$corr_leaf_turnover_temp > 0, "group1",
                                                  ifelse(output$corr_leaf_turnover_precip < 0 &  output$corr_leaf_turnover_sun == 0 &   output$corr_leaf_turnover_temp > 0, "group1",
                                                         ifelse(output$corr_leaf_turnover_precip < 0 &  output$corr_leaf_turnover_sun > 0 &    output$corr_leaf_turnover_temp == 0, "group1",
                                                                ifelse(output$corr_leaf_turnover_precip < 0 &  output$corr_leaf_turnover_sun == 0 &   output$corr_leaf_turnover_temp == 0, "group1",

                                                                       ifelse(output$corr_leaf_turnover_precip > 0 &  output$corr_leaf_turnover_sun < 0 &    output$corr_leaf_turnover_temp < 0, "group2",
                                                                              ifelse(output$corr_leaf_turnover_precip == 0 & output$corr_leaf_turnover_sun < 0 &    output$corr_leaf_turnover_temp < 0, "group2",
                                                                                     ifelse(output$corr_leaf_turnover_precip == 0 & output$corr_leaf_turnover_sun < 0 &    output$corr_leaf_turnover_temp == 0, "group2",
                                                                                            ifelse(output$corr_leaf_turnover_precip == 0 & output$corr_leaf_turnover_sun == 0 &   output$corr_leaf_turnover_temp < 0, "group2",
                                                                                                   ifelse(output$corr_leaf_turnover_precip > 0 &  output$corr_leaf_turnover_sun == 0 &   output$corr_leaf_turnover_temp < 0, "group2",
                                                                                                          ifelse(output$corr_leaf_turnover_precip > 0 &  output$corr_leaf_turnover_sun < 0 &    output$corr_leaf_turnover_temp == 0, "group2",
                                                                                                                 ifelse(output$corr_leaf_turnover_precip > 0 &  output$corr_leaf_turnover_sun == 0 &   output$corr_leaf_turnover_temp == 0, "group2", NA))))))))))))))
output$corr_leaf_dormancy_precip <- ifelse(is.na(output$corr_leaf_dormancy_precip), 0, output$corr_leaf_dormancy_precip)
output$corr_leaf_dormancy_sun <- ifelse(is.na(output$corr_leaf_dormancy_sun), 0, output$corr_leaf_dormancy_sun)
output$corr_leaf_dormancy_temp <- ifelse(is.na(output$corr_leaf_dormancy_temp), 0, output$corr_leaf_dormancy_temp)


output$dorm_groups <- ifelse(output$corr_leaf_dormancy_precip < 0 &  output$corr_leaf_dormancy_sun > 0 &    output$corr_leaf_dormancy_temp > 0, "group1",
                             ifelse(output$corr_leaf_dormancy_precip == 0 & output$corr_leaf_dormancy_sun > 0 &    output$corr_leaf_dormancy_temp > 0, "group1",
                                    ifelse(output$corr_leaf_dormancy_precip == 0 & output$corr_leaf_dormancy_sun > 0 &    output$corr_leaf_dormancy_temp == 0, "group1",
                                           ifelse(output$corr_leaf_dormancy_precip == 0 & output$corr_leaf_dormancy_sun == 0 &   output$corr_leaf_dormancy_temp > 0, "group1",
                                                  ifelse(output$corr_leaf_dormancy_precip < 0 &  output$corr_leaf_dormancy_sun == 0 &   output$corr_leaf_dormancy_temp > 0, "group1",
                                                         ifelse(output$corr_leaf_dormancy_precip < 0 &  output$corr_leaf_dormancy_sun > 0 &    output$corr_leaf_dormancy_temp == 0, "group1",
                                                                ifelse(output$corr_leaf_dormancy_precip < 0 &  output$corr_leaf_dormancy_sun == 0 &   output$corr_leaf_dormancy_temp == 0, "group1",

                                                                       ifelse(output$corr_leaf_dormancy_precip > 0 &  output$corr_leaf_dormancy_sun < 0 &    output$corr_leaf_dormancy_temp < 0, "group2",
                                                                              ifelse(output$corr_leaf_dormancy_precip == 0 & output$corr_leaf_dormancy_sun < 0 &    output$corr_leaf_dormancy_temp < 0, "group2",
                                                                                     ifelse(output$corr_leaf_dormancy_precip == 0 & output$corr_leaf_dormancy_sun < 0 &    output$corr_leaf_dormancy_temp == 0, "group2",
                                                                                            ifelse(output$corr_leaf_dormancy_precip == 0 & output$corr_leaf_dormancy_sun == 0 &   output$corr_leaf_dormancy_temp < 0, "group2",
                                                                                                   ifelse(output$corr_leaf_dormancy_precip > 0 &  output$corr_leaf_dormancy_sun == 0 &   output$corr_leaf_dormancy_temp < 0, "group2",
                                                                                                          ifelse(output$corr_leaf_dormancy_precip > 0 &  output$corr_leaf_dormancy_sun < 0 &    output$corr_leaf_dormancy_temp == 0, "group2",
                                                                                                                 ifelse(output$corr_leaf_dormancy_precip > 0 &  output$corr_leaf_dormancy_sun == 0 &   output$corr_leaf_dormancy_temp == 0, "group2", NA))))))))))))))

# back to NA
output$corr_leaf_turnover_precip <- ifelse(output$corr_leaf_turnover_precip == 0, NA, output$corr_leaf_turnover_precip)
output$corr_leaf_turnover_sun <- ifelse(output$corr_leaf_turnover_sun == 0, NA, output$corr_leaf_turnover_sun)
output$corr_leaf_turnover_temp <- ifelse(output$corr_leaf_turnover_temp == 0, NA, output$corr_leaf_turnover_temp)
output$corr_leaf_dormancy_precip <- ifelse(output$corr_leaf_dormancy_precip == 0, NA, output$corr_leaf_dormancy_precip)
output$corr_leaf_dormancy_sun <- ifelse(output$corr_leaf_dormancy_sun == 0, NA, output$corr_leaf_dormancy_sun)
output$corr_leaf_dormancy_temp <- ifelse(output$corr_leaf_dormancy_temp == 0, NA, output$corr_leaf_dormancy_temp)

# write.table(output, "data/timeseries_correlations_phase2.csv",
#             quote = FALSE,
#             col.names = TRUE,
#             row.names = FALSE,
#             sep = ",")

output <- merge(overview, output, by = "species_full", all.x = TRUE)
# output <- merge(output, dorm_events, by = "species_full", all.x = TRUE)
# output <- merge(output, turn_events, by = "species_full", all.x = TRUE)



#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# corr summary figure
# only species that have events of dormancy or turnover
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------


# dormancy
df_dorm <- output %>%
  filter(total_nr_events_leaf_dormancy > 0)

df_dorm$precip_phase <- ifelse(df_dorm$total_nr_events_leaf_dormancy <5, "few-events",
                               ifelse(df_dorm$corr_leaf_dormancy_precip_timing == "0" & df_dorm$corr_leaf_dormancy_precip < 0, "no-lag-neg",
                                      ifelse(df_dorm$corr_leaf_dormancy_precip_timing == "0" & df_dorm$corr_leaf_dormancy_precip > 0, "no-lag-pos",
                                             ifelse(df_dorm$corr_leaf_dormancy_precip_timing %in% c("-1","-2","-3") & df_dorm$corr_leaf_dormancy_precip < 0, "lag-neg",
                                                    ifelse(df_dorm$corr_leaf_dormancy_precip_timing %in% c("-1","-2","-3") & df_dorm$corr_leaf_dormancy_precip > 0, "lag-pos",
                                                           NA)))))
df_dorm$precip_phase <- ifelse(is.na(df_dorm$precip_phase) & df_dorm$total_nr_events_leaf_dormancy >= 5, "h-no-corr", df_dorm$precip_phase)


df_dorm$insol_phase <- ifelse(df_dorm$total_nr_events_leaf_dormancy <5, "few-events",
                              ifelse(df_dorm$corr_leaf_dormancy_sun_timing == "0" & df_dorm$corr_leaf_dormancy_sun < 0, "no-lag-neg",
                                     ifelse(df_dorm$corr_leaf_dormancy_sun_timing == "0" & df_dorm$corr_leaf_dormancy_sun > 0, "no-lag-pos",
                                            ifelse(df_dorm$corr_leaf_dormancy_sun_timing %in% c("-1","-2","-3") & df_dorm$corr_leaf_dormancy_sun < 0, "lag-neg",
                                                   ifelse(df_dorm$corr_leaf_dormancy_sun_timing %in% c("-1","-2","-3") & df_dorm$corr_leaf_dormancy_sun > 0, "lag-pos",
                                                          NA)))))
df_dorm$insol_phase <- ifelse(is.na(df_dorm$insol_phase) & df_dorm$total_nr_events_leaf_dormancy >= 5, "h-no-corr", df_dorm$insol_phase)

df_dorm$tmax_phase <- ifelse(df_dorm$total_nr_events_leaf_dormancy <5, "few-events",
                             ifelse(df_dorm$corr_leaf_dormancy_temp_timing == "0" & df_dorm$corr_leaf_dormancy_temp < 0, "no-lag-neg",
                                    ifelse(df_dorm$corr_leaf_dormancy_temp_timing == "0" & df_dorm$corr_leaf_dormancy_temp > 0, "no-lag-pos",
                                           ifelse(df_dorm$corr_leaf_dormancy_temp_timing %in% c("-1","-2","-3") & df_dorm$corr_leaf_dormancy_temp < 0, "lag-neg",
                                                  ifelse(df_dorm$corr_leaf_dormancy_temp_timing %in% c("-1","-2","-3") & df_dorm$corr_leaf_dormancy_temp > 0, "lag-pos",
                                                         NA)))))
df_dorm$tmax_phase <- ifelse(is.na(df_dorm$tmax_phase) & df_dorm$total_nr_events_leaf_dormancy >= 5, "h-no-corr", df_dorm$tmax_phase)


# turnover
df_turn <- output %>%
  filter(total_nr_events_leaf_turnover > 0)

df_turn$precip_phase <- ifelse(df_turn$total_nr_events_leaf_turnover <5, "few-events",
                               ifelse(df_turn$corr_leaf_turnover_precip_timing == "0" & df_turn$corr_leaf_turnover_precip < 0, "no-lag-neg",
                                      ifelse(df_turn$corr_leaf_turnover_precip_timing == "0" & df_turn$corr_leaf_turnover_precip > 0, "no-lag-pos",
                                             ifelse(df_turn$corr_leaf_turnover_precip_timing %in% c("-1","-2","-3") & df_turn$corr_leaf_turnover_precip < 0, "lag-neg",
                                                    ifelse(df_turn$corr_leaf_turnover_precip_timing %in% c("-1","-2","-3") & df_turn$corr_leaf_turnover_precip > 0, "lag-pos",
                                                           NA)))))
df_turn$precip_phase <- ifelse(is.na(df_turn$precip_phase) & df_turn$total_nr_events_leaf_turnover >= 5, "h-no-corr", df_turn$precip_phase)


df_turn$insol_phase <- ifelse(df_turn$total_nr_events_leaf_turnover <5, "few-events",
                              ifelse(df_turn$corr_leaf_turnover_sun_timing == "0" & df_turn$corr_leaf_turnover_sun < 0, "no-lag-neg",
                                     ifelse(df_turn$corr_leaf_turnover_sun_timing == "0" & df_turn$corr_leaf_turnover_sun > 0, "no-lag-pos",
                                            ifelse(df_turn$corr_leaf_turnover_sun_timing %in% c("-1","-2","-3") & df_turn$corr_leaf_turnover_sun < 0, "lag-neg",
                                                   ifelse(df_turn$corr_leaf_turnover_sun_timing %in% c("-1","-2","-3") & df_turn$corr_leaf_turnover_sun > 0, "lag-pos",
                                                          NA)))))
df_turn$insol_phase <- ifelse(is.na(df_turn$insol_phase) & df_turn$total_nr_events_leaf_turnover >= 5, "h-no-corr", df_turn$insol_phase)

df_turn$tmax_phase <- ifelse(df_turn$total_nr_events_leaf_turnover <5, "few-events",
                             ifelse(df_turn$corr_leaf_turnover_temp_timing == "0" & df_turn$corr_leaf_turnover_temp < 0, "no-lag-neg",
                                    ifelse(df_turn$corr_leaf_turnover_temp_timing == "0" & df_turn$corr_leaf_turnover_temp > 0, "no-lag-pos",
                                           ifelse(df_turn$corr_leaf_turnover_temp_timing %in% c("-1","-2","-3") & df_turn$corr_leaf_turnover_temp < 0, "lag-neg",
                                                  ifelse(df_turn$corr_leaf_turnover_temp_timing %in% c("-1","-2","-3") & df_turn$corr_leaf_turnover_temp > 0, "lag-pos",
                                                         NA)))))
df_turn$tmax_phase <- ifelse(is.na(df_turn$tmax_phase) & df_turn$total_nr_events_leaf_turnover >= 5, "h-no-corr", df_turn$tmax_phase)

#-----------------------------
# evergreen - dormancy
#-----------------------------
df_ever_dorm <- df_dorm %>%
  filter(grepl("evergreen",deciduousness))

# for each climate variable
# get % of species with cross-correlations in different phases and signs
counts_ever_dorm <- length(df_ever_dorm$species_full)
ever_dorm_prec <- df_ever_dorm %>%
  group_by(precip_phase) %>%
  dplyr::summarise(value = length(corr_leaf_dormancy_precip)/counts_ever_dorm*100) %>%
  dplyr::rename(relation = precip_phase) %>%
  add_column(variable = "precipitation")
ever_dorm_sun <- df_ever_dorm %>%
  group_by(insol_phase) %>%
  dplyr::summarise(value = length(corr_leaf_dormancy_sun)/counts_ever_dorm*100) %>%
  dplyr::rename(relation = insol_phase) %>%
  add_column(variable = "sunhours")
ever_dorm_temp <- df_ever_dorm %>%
  group_by(tmax_phase) %>%
  dplyr::summarise(value = length(corr_leaf_dormancy_temp)/counts_ever_dorm*100) %>%
  dplyr::rename(relation = tmax_phase) %>%
  add_column(variable = "tmax")
# rbind all
ever_dorm_summary <- rbind(ever_dorm_prec, ever_dorm_sun, ever_dorm_temp)
# % of species with 'few-events' and 'no-corr' will be added
# 1/2 at both the negative and positive correlation ends
# small hack to divide by 2 with seperate labels for relation
ever_dorm_summary$value <- ifelse(ever_dorm_summary$relation %in% c("few-events","h-no-corr"),
                                  ever_dorm_summary$value / 2, ever_dorm_summary$value)
hack <- ever_dorm_summary %>%
  filter(relation %in% "few-events") %>%
  mutate(relation = "few-events-neg")
hack2 <- ever_dorm_summary %>%
  filter(relation %in% "h-no-corr") %>%
  mutate(relation = "h-no-corr-neg")
ever_dorm_summary <- rbind(ever_dorm_summary, hack, hack2)
# ever_dorm_summary <- ever_dorm_summary %>%
#   add_column(class = "1evergeen-dormancy")
ever_dorm_summary <- ever_dorm_summary %>%
  add_column(class = "1evergeen",
             pheno = "dormancy")
#-----------------------------
# evergreen - turnover
#-----------------------------
df_ever_turn <- df_turn %>%
  filter(grepl("evergreen",deciduousness))

# for each climate variable
# get % of species with cross-correlations in different phases and signs
counts_ever_turn <- length(df_ever_turn$species_full)
ever_turn_prec <- df_ever_turn %>%
  group_by(precip_phase) %>%
  dplyr::summarise(value = length(corr_leaf_turnover_precip)/counts_ever_turn*100) %>%
  dplyr::rename(relation = precip_phase) %>%
  add_column(variable = "precipitation")
ever_turn_sun <- df_ever_turn %>%
  group_by(insol_phase) %>%
  dplyr::summarise(value = length(corr_leaf_turnover_sun)/counts_ever_turn*100) %>%
  dplyr::rename(relation = insol_phase) %>%
  add_column(variable = "sunhours")
ever_turn_temp <- df_ever_turn %>%
  group_by(tmax_phase) %>%
  dplyr::summarise(value = length(corr_leaf_turnover_temp)/counts_ever_turn*100) %>%
  dplyr::rename(relation = tmax_phase) %>%
  add_column(variable = "tmax")
# rbind all
ever_turn_summary <- rbind(ever_turn_prec, ever_turn_sun, ever_turn_temp)
# % of species with 'few-events' and 'no-corr' will be added
# 1/2 at both the negative and positive correlation ends
# small hack to divide by 2 with seperate labels for relation
ever_turn_summary$value <- ifelse(ever_turn_summary$relation %in% c("few-events","h-no-corr"),
                                  ever_turn_summary$value / 2, ever_turn_summary$value)
hack <- ever_turn_summary %>%
  filter(relation %in% "few-events") %>%
  mutate(relation = "few-events-neg")
hack2 <- ever_turn_summary %>%
  filter(relation %in% "h-no-corr") %>%
  mutate(relation = "h-no-corr-neg")
ever_turn_summary <- rbind(ever_turn_summary, hack, hack2)
# ever_turn_summary <- ever_turn_summary %>%
#   add_column(class = "1evergeen-turnover")
ever_turn_summary <- ever_turn_summary %>%
  add_column(class = "1evergeen",
             pheno = "turnover")


#-----------------------------
# deciduous - dormancy
#-----------------------------
df_dec_dorm <- df_dorm %>%
  filter(grepl("deciduous",deciduousness))


# for each climate variable
# get % of species with cross-correlations in different phases and signs
counts_dec_dorm <- length(df_dec_dorm$species_full)
dec_dorm_prec <- df_dec_dorm %>%
  group_by(precip_phase) %>%
  dplyr::summarise(value = length(corr_leaf_dormancy_precip)/counts_dec_dorm*100) %>%
  dplyr::rename(relation = precip_phase) %>%
  add_column(variable = "precipitation")
dec_dorm_sun <- df_dec_dorm %>%
  group_by(insol_phase) %>%
  dplyr::summarise(value = length(corr_leaf_dormancy_sun)/counts_dec_dorm*100) %>%
  dplyr::rename(relation = insol_phase) %>%
  add_column(variable = "sunhours")
dec_dorm_temp <- df_dec_dorm %>%
  group_by(tmax_phase) %>%
  dplyr::summarise(value = length(corr_leaf_dormancy_temp)/counts_dec_dorm*100) %>%
  dplyr::rename(relation = tmax_phase) %>%
  add_column(variable = "tmax")
# rbind all
dec_dorm_summary <- rbind(dec_dorm_prec, dec_dorm_sun, dec_dorm_temp)
# % of species with 'few-events' and 'no-corr' will be added
# 1/2 at both the negative and positive correlation ends
# small hack to divide by 2 with seperate labels for relation
dec_dorm_summary$value <- ifelse(dec_dorm_summary$relation %in% c("few-events","h-no-corr"),
                                  dec_dorm_summary$value / 2, dec_dorm_summary$value)
hack <- dec_dorm_summary %>%
  filter(relation %in% "few-events") %>%
  mutate(relation = "few-events-neg")
hack2 <- dec_dorm_summary %>%
  filter(relation %in% "h-no-corr") %>%
  mutate(relation = "h-no-corr-neg")
dec_dorm_summary <- rbind(dec_dorm_summary, hack, hack2)
# dec_dorm_summary <- dec_dorm_summary %>%
#   add_column(class = "2deciduous-dormancy")
dec_dorm_summary <- dec_dorm_summary %>%
  add_column(class = "deciduous",
             pheno = "dormancy")


#-----------------------------
# deciduous - turnover
#-----------------------------
df_dec_turn <- df_turn %>%
  filter(grepl("deciduous",deciduousness))

# for each climate variable
# get % of species with cross-correlations in different phases and signs
counts_dec_turn <- length(df_dec_turn$species_full)
dec_turn_prec <- df_dec_turn %>%
  group_by(precip_phase) %>%
  dplyr::summarise(value = length(corr_leaf_turnover_precip)/counts_dec_turn*100) %>%
  dplyr::rename(relation = precip_phase) %>%
  add_column(variable = "precipitation")
dec_turn_sun <- df_dec_turn %>%
  group_by(insol_phase) %>%
  dplyr::summarise(value = length(corr_leaf_turnover_sun)/counts_dec_turn*100) %>%
  dplyr::rename(relation = insol_phase) %>%
  add_column(variable = "sunhours")
dec_turn_temp <- df_dec_turn %>%
  group_by(tmax_phase) %>%
  dplyr::summarise(value = length(corr_leaf_turnover_temp)/counts_dec_turn*100) %>%
  dplyr::rename(relation = tmax_phase) %>%
  add_column(variable = "tmax")
# rbind all
dec_turn_summary <- rbind(dec_turn_prec, dec_turn_sun, dec_turn_temp)
# % of species with 'few-events' and 'no-corr' will be added
# 1/2 at both the negative and positive correlation ends
# small hack to divide by 2 with seperate labels for relation
dec_turn_summary$value <- ifelse(dec_turn_summary$relation %in% c("few-events","h-no-corr"),
                                  dec_turn_summary$value / 2, dec_turn_summary$value)
hack <- dec_turn_summary %>%
  filter(relation %in% "few-events") %>%
  mutate(relation = "few-events-neg")
hack2 <- dec_turn_summary %>%
  filter(relation %in% "h-no-corr") %>%
  mutate(relation = "h-no-corr-neg")
dec_turn_summary <- rbind(dec_turn_summary, hack, hack2)
# dec_turn_summary <- dec_turn_summary %>%
#   add_column(class = "2deciduous-turnover")
dec_turn_summary <- dec_turn_summary %>%
  add_column(class = "deciduous",
             pheno = "turnover")


corr_summary <- rbind(ever_dorm_summary, ever_turn_summary, dec_dorm_summary, dec_turn_summary)


##################################################################
new_labels <- c("1evergeen" = "", "deciduous" = "")

p_all <- ggplot(corr_summary,
               aes(x = variable,
                   y = ifelse(relation %in% c("no-lag-neg","lag-neg","few-events-neg","h-no-corr-neg"), -value, value),
                   fill = relation)) +
  geom_col() +
  scale_y_continuous(limits = c(-73,68),
                     breaks = c(-50,-25,0,25,50),
                     sec.axis = dup_axis(name = "",
                                         breaks = c(-25,25),
                                         labels = c("negative correlations", "positive correlations"))) +
  scale_x_discrete(limits = rev(levels(as.factor(corr_summary$variable)))) +
  coord_flip() +
  # scale_fill_manual(values = c("no-lag-pos" = "#018571",
  #                              "no-lag-neg" = "#018571",
  #                              "lag-pos" = "#80cdc1",
  #                              "lag-neg" = "#80cdc1",
  #                              "few-events" = "grey90",
  #                              "few-events-neg" = "grey90",
  #                              "h-no-corr" = "#dfc27d",
  #                              "h-no-corr-neg" = "#dfc27d")) +
  scale_fill_manual(values = c("grey90", "grey90","#dfc27d", "#dfc27d", "#80cdc1", "#80cdc1", "#018571","#018571"),
                    breaks = c("no-lag-neg","lag-neg","h-no-corr-neg","few-events-neg"),
                    labels = c(" in-phase   "," lag   "," no correlation   "," too few events   ")) +
  geom_hline(yintercept = 0, color =c("grey60"), size = 0.5) +
  labs(y = "% species with neg. or pos. correlations",
       x = "")  +
  theme(text = element_text(family = "Lato", color = "#22211d"),
        panel.grid.major.x = element_line(colour = "grey60"),#, size = 0.3),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(size = 11),
        strip.background = element_rect(fill="white"),
        strip.placement = "outside",
        # axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(),
        axis.title.y = element_text(),
        axis.title.x = element_text(),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm")
  ) +
  facet_wrap( ~ class + pheno,
              ncol = 1,
             labeller = labeller(class = new_labels),
             strip.position = "left")
# p_all

ann_text <- data.frame(x = 3,
                        y = 53,
                       relation = "few-events-neg",
                        class = c("1evergeen","1evergeen","deciduous","deciduous"),
                        pheno = c("dormancy", "turnover","dormancy", "turnover"),
                        notes = c(paste("n = ",counts_ever_dorm),
                                  paste("n = ",counts_ever_turn),
                                  paste("n = ",counts_dec_dorm),
                                  paste("n = ",counts_dec_turn)))
p1 <- p_all +
  geom_text(data = ann_text,
            aes(x, y, label = notes),
            size = 3,
            hjust = 0,
            family = "Lato", color = "#22211d")
# p1
p2 <- grid.arrange(arrangeGrob(p1,
                                right = textGrob("deciduous                               evergreen",
                                                 gp=gpar(family = "Lato", color = "#22211d",fontsize=11),
                                                 rot = 90,
                                                 hjust = 0.4,
                                                 vjust = -57)))





