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
#----- source required functions -----------------------------------------------#
source("R/timeline_gap_fill.R")
source("R/consec_years.R")
source("R/climate_ccf_function.R")
#-------------------------------------------------------------------------------#


#----------------------------------------------------------------------
#--------   Phenology data - species correction Meise   ---------------
#----------------------------------------------------------------------
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
df$join_id <- paste0("R",df$image,"-",df$image_row)
metadata <- read.csv("data/phenology_archives_species_long_format_20200324.csv",
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

# sum events for each id, each year, across phenophases
# years with zero observations across phenophases are possibly not observed
empty_years <- data %>%
  group_by(species_full,join_id,year) %>%
  dplyr::summarise(check_empty_years = sum(value))
data <- merge(data, empty_years, by = c("join_id","species_full","year"), all.x = TRUE)
data <- data %>%
  filter(check_empty_years > 0)
#----------------------------------------------------------------------
# only select parameters you need, more clear structure to work with
data <- data %>%
  select(species_full,
         id,
         phenophase,
         year,
         week,
         value)
#----------------------------------------------------------------------
rm(df,metadata, empty_years)
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
overview <- overview %>%
  dplyr::select(species_full,
                deciduousness,
                # site_years_with_leaf_dormancy, # in summary figure, only species with events are included
                # site_years_with_leaf_turnover,
                total_nr_events_leaf_dormancy, # species with too few events (< 5) not included for cross correlation analysis
                total_nr_events_leaf_turnover)

species_list <- overview$species_full

#--- leaf turnover ----------------------------------------------------------------------------------------------------
# for selected species and phenophase: get extended timelines at ID level with 2 year-gaps filled with zero
timelines_id_turn <- two_year_gaps(data = data,
                              species_name = species_list,
                              pheno = "leaf_turnover")
# for fourier, no gaps in timelines allowed
# get longest consecutive timeline; at species level
timelines_sp_consec_turn <- consecutive_timeline_sp(data = timelines_id_turn,
                                               species_name = species_list,
                                               pheno = "leaf_turnover")
timelines_sp_consec_turn <- merge(overview, timelines_sp_consec_turn, by = "species_full", all.x = TRUE)
#--- leaf dormancy ----------------------------------------------------------------------------------------------------
timelines_id_dorm <- two_year_gaps(data = data,
                                   species_name = species_list,
                                   pheno = "leaf_dormancy")
timelines_sp_consec_dorm <- consecutive_timeline_sp(data = timelines_id_dorm,
                                                    species_name = species_list,
                                                    pheno = "leaf_dormancy")
timelines_sp_consec_dorm <- merge(overview, timelines_sp_consec_dorm, by = "species_full", all.x = TRUE)
#----------------------------------------------------------------------------------------------------------------------

test_event_count <- function(
  data = data,
  species_name = "Afzelia bipindensis"){
  el_output <- data.frame()

  for (j in 1:length(species_name)){
    data_subset <- data %>%
      filter(species_full %in% species_name[j])
    # sort dataframe according to date
    data_subset <- data_subset %>%
      dplyr::arrange(date)
    # get first differences
    diff_values <- diff(data_subset$scaled_value)
    # get matching info
    start <- data_subset[which(diff_values == 1) + 1,]
    a <- length(start$scaled_value)
    el_sp <- data.frame(species_full = species_name[j],
                        nr_events = length(start$scaled_value))
    el_output <- rbind(el_output,el_sp)
  }
  return(el_output)
}

dorm_events <- test_event_count(data = timelines_sp_consec_dorm,
                         species_name = unique(timelines_sp_consec_dorm$species_full))
dorm_events <- dorm_events %>%
  dplyr::rename(nr_dorm_events_consec = nr_events)

turn_events <- test_event_count(data = timelines_sp_consec_turn,
                              species_name = unique(timelines_sp_consec_turn$species_full))
turn_events <- turn_events %>%
  dplyr::rename(nr_turn_events_consec = nr_events)

#----------------------------------------------------------------------------------------------------------------------
#--- cross correlations climate - phenology timeseries ----------------------------------------------------------------
#--- phenology timeseries with too few events (< 5) are not analysed within the function (set to NA) ------------------
#----------------------------------------------------------------------------------------------------------------------
# leaf turnover
crosscorr <- climate_ccf(data = timelines_sp_consec_turn,
                         climate = climate,
                         climate.variable = "precip",
                         species_name = species_list,
                         pheno = "leaf_turnover")
crosscorr2 <- climate_ccf(data = timelines_sp_consec_turn,
                          climate = climate,
                          climate.variable = "sun",
                          species_name = species_list,
                          pheno = "leaf_turnover")
crosscorr3 <- climate_ccf(data = timelines_sp_consec_turn,
                          climate = climate,
                          climate.variable = "temp",
                          species_name = species_list,
                          pheno = "leaf_turnover")
# leaf dormancy
crosscorr4 <- climate_ccf(data = timelines_sp_consec_dorm,
                          climate = climate,
                          climate.variable = "precip",
                          species_name = species_list,
                          pheno = "leaf_dormancy")
crosscorr5 <- climate_ccf(data = timelines_sp_consec_dorm,
                          climate = climate,
                          climate.variable = "sun",
                          species_name = species_list,
                          pheno = "leaf_dormancy")
crosscorr6 <- climate_ccf(data = timelines_sp_consec_dorm,
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

# write.table(output, "data/timeseries_correlations_phase2.csv",
#             quote = FALSE,
#             col.names = TRUE,
#             row.names = FALSE,
#             sep = ",")



output <- merge(overview, output, by = "species_full", all.x = TRUE)
output <- merge(output, dorm_events, by = "species_full", all.x = TRUE)
output <- merge(output, turn_events, by = "species_full", all.x = TRUE)

test <- output$total_nr_events_leaf_dormancy - output$nr_dorm_events_consec
test2 <- output$total_nr_events_leaf_turnover - output$nr_turn_events_consec
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# corr summary figure
# only species that have events of dormancy or turnover
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

test <- df_dorm %>%
  dplyr::select(species_full,
                total_nr_events_leaf_dormancy,
                nr_dorm_events_consec,
                precip_phase,
                insol_phase)

# dormancy
df_dorm <- output %>%
  filter(nr_dorm_events_consec > 0)

df_dorm$precip_phase <- ifelse(df_dorm$nr_dorm_events_consec <5, "few-events",
                               ifelse(df_dorm$corr_leaf_dormancy_precip_timing == "0" & df_dorm$corr_leaf_dormancy_precip < 0, "no-lag-neg",
                                      ifelse(df_dorm$corr_leaf_dormancy_precip_timing == "0" & df_dorm$corr_leaf_dormancy_precip > 0, "no-lag-pos",
                                             ifelse(df_dorm$corr_leaf_dormancy_precip_timing %in% c("-1","-2","-3") & df_dorm$corr_leaf_dormancy_precip < 0, "lag-neg",
                                                    ifelse(df_dorm$corr_leaf_dormancy_precip_timing %in% c("-1","-2","-3") & df_dorm$corr_leaf_dormancy_precip > 0, "lag-pos",
                                                           NA)))))
df_dorm$precip_phase <- ifelse(is.na(df_dorm$precip_phase) & df_dorm$nr_dorm_events_consec >= 5, "h-no-corr", df_dorm$precip_phase)


df_dorm$insol_phase <- ifelse(df_dorm$nr_dorm_events_consec <5, "few-events",
                              ifelse(df_dorm$corr_leaf_dormancy_sun_timing == "0" & df_dorm$corr_leaf_dormancy_sun < 0, "no-lag-neg",
                                     ifelse(df_dorm$corr_leaf_dormancy_sun_timing == "0" & df_dorm$corr_leaf_dormancy_sun > 0, "no-lag-pos",
                                            ifelse(df_dorm$corr_leaf_dormancy_sun_timing %in% c("-1","-2","-3") & df_dorm$corr_leaf_dormancy_sun < 0, "lag-neg",
                                                   ifelse(df_dorm$corr_leaf_dormancy_sun_timing %in% c("-1","-2","-3") & df_dorm$corr_leaf_dormancy_sun > 0, "lag-pos",
                                                          NA)))))
df_dorm$insol_phase <- ifelse(is.na(df_dorm$insol_phase) & df_dorm$nr_dorm_events_consec >= 5, "h-no-corr", df_dorm$insol_phase)

df_dorm$tmax_phase <- ifelse(df_dorm$nr_dorm_events_consec <5, "few-events",
                             ifelse(df_dorm$corr_leaf_dormancy_temp_timing == "0" & df_dorm$corr_leaf_dormancy_temp < 0, "no-lag-neg",
                                    ifelse(df_dorm$corr_leaf_dormancy_temp_timing == "0" & df_dorm$corr_leaf_dormancy_temp > 0, "no-lag-pos",
                                           ifelse(df_dorm$corr_leaf_dormancy_temp_timing %in% c("-1","-2","-3") & df_dorm$corr_leaf_dormancy_temp < 0, "lag-neg",
                                                  ifelse(df_dorm$corr_leaf_dormancy_temp_timing %in% c("-1","-2","-3") & df_dorm$corr_leaf_dormancy_temp > 0, "lag-pos",
                                                         NA)))))
df_dorm$tmax_phase <- ifelse(is.na(df_dorm$tmax_phase) & df_dorm$nr_dorm_events_consec >= 5, "h-no-corr", df_dorm$tmax_phase)


# turnover
df_turn <- output %>%
  filter(nr_turn_events_consec > 0)

df_turn$precip_phase <- ifelse(df_turn$nr_turn_events_consec <5, "few-events",
                               ifelse(df_turn$corr_leaf_turnover_precip_timing == "0" & df_turn$corr_leaf_turnover_precip < 0, "no-lag-neg",
                                      ifelse(df_turn$corr_leaf_turnover_precip_timing == "0" & df_turn$corr_leaf_turnover_precip > 0, "no-lag-pos",
                                             ifelse(df_turn$corr_leaf_turnover_precip_timing %in% c("-1","-2","-3") & df_turn$corr_leaf_turnover_precip < 0, "lag-neg",
                                                    ifelse(df_turn$corr_leaf_turnover_precip_timing %in% c("-1","-2","-3") & df_turn$corr_leaf_turnover_precip > 0, "lag-pos",
                                                           NA)))))
df_turn$precip_phase <- ifelse(is.na(df_turn$precip_phase) & df_turn$nr_turn_events_consec >= 5, "h-no-corr", df_turn$precip_phase)


df_turn$insol_phase <- ifelse(df_turn$nr_turn_events_consec <5, "few-events",
                              ifelse(df_turn$corr_leaf_turnover_sun_timing == "0" & df_turn$corr_leaf_turnover_sun < 0, "no-lag-neg",
                                     ifelse(df_turn$corr_leaf_turnover_sun_timing == "0" & df_turn$corr_leaf_turnover_sun > 0, "no-lag-pos",
                                            ifelse(df_turn$corr_leaf_turnover_sun_timing %in% c("-1","-2","-3") & df_turn$corr_leaf_turnover_sun < 0, "lag-neg",
                                                   ifelse(df_turn$corr_leaf_turnover_sun_timing %in% c("-1","-2","-3") & df_turn$corr_leaf_turnover_sun > 0, "lag-pos",
                                                          NA)))))
df_turn$insol_phase <- ifelse(is.na(df_turn$insol_phase) & df_turn$nr_turn_events_consec >= 5, "h-no-corr", df_turn$insol_phase)

df_turn$tmax_phase <- ifelse(df_turn$nr_turn_events_consec <5, "few-events",
                             ifelse(df_turn$corr_leaf_turnover_temp_timing == "0" & df_turn$corr_leaf_turnover_temp < 0, "no-lag-neg",
                                    ifelse(df_turn$corr_leaf_turnover_temp_timing == "0" & df_turn$corr_leaf_turnover_temp > 0, "no-lag-pos",
                                           ifelse(df_turn$corr_leaf_turnover_temp_timing %in% c("-1","-2","-3") & df_turn$corr_leaf_turnover_temp < 0, "lag-neg",
                                                  ifelse(df_turn$corr_leaf_turnover_temp_timing %in% c("-1","-2","-3") & df_turn$corr_leaf_turnover_temp > 0, "lag-pos",
                                                         NA)))))
df_turn$tmax_phase <- ifelse(is.na(df_turn$tmax_phase) & df_turn$nr_turn_events_consec >= 5, "h-no-corr", df_turn$tmax_phase)

## time lag of dormancy and turnover together for insol and tmax
# test <- df_dorm %>%
#   filter(insol_phase %in% "lag-neg")
# test2 <- df_turn %>%
#   filter(insol_phase %in% "lag-neg")
# a <- c(test$corr_leaf_dormancy_sun_timing,test2$corr_leaf_turnover_sun_timing)
# mean(as.numeric(a))
# sd(as.numeric(a))
# test3 <- df_dorm %>%
#   filter(tmax_phase %in% "lag-neg")
# test4 <- df_turn %>%
#   filter(tmax_phase %in% "lag-neg")
# b <- c(test3$corr_leaf_dormancy_temp_timing,test4$corr_leaf_turnover_temp_timing)
# mean(as.numeric(b))
# sd(as.numeric(b))


#-----------------------------
# evergreen - dormancy
#-----------------------------
df_ever_dorm <- df_dorm %>%
  filter(grepl("evergreen",deciduousness))

ed_precip <- as.data.frame(tapply(df_ever_dorm$corr_leaf_dormancy_precip, list(df_ever_dorm$precip_phase), length)) #ed = ever dorm
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

ed_insol <- as.data.frame(tapply(df_ever_dorm$corr_leaf_dormancy_sun, list(df_ever_dorm$insol_phase), length)) #ed = ever dorm
counts_ever_dorm <- length(df_ever_dorm$species_full)
colnames(ed_insol) <- "value"
ed_insol$value <- ed_insol$value / counts_ever_dorm *100
ed_insol$relation <- rownames(ed_insol)
ed_insol$variable <- "sunhours"

ed_tmax <- as.data.frame(tapply(df_ever_dorm$corr_leaf_dormancy_temp, list(df_ever_dorm$tmax_phase), length)) #ed = ever dorm
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
  scale_y_continuous(limits = c(-72,72),
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
        plot.margin = unit(c(0.5,0,-0.1,0.5),"cm")
  )



#-----------------------------
# evergreen - turnover
#-----------------------------
df_ever_turn <- df_turn %>%
  filter(grepl("evergreen",deciduousness))

et_precip <- as.data.frame(tapply(df_ever_turn$corr_leaf_turnover_precip, list(df_ever_turn$precip_phase), length)) #et = ever dorm
counts_ever_turn <- length(df_ever_turn$species_full)
colnames(et_precip) <- "value"
et_precip$value <- et_precip$value / counts_ever_turn *100
et_precip$relation <- rownames(et_precip)
et_precip$variable <- "precipitation"

et_insol <- as.data.frame(tapply(df_ever_turn$corr_leaf_turnover_sun, list(df_ever_turn$insol_phase), length)) #et = ever dorm
counts_ever_turn <- length(df_ever_turn$species_full)
colnames(et_insol) <- "value"
et_insol$value <- et_insol$value / counts_ever_turn *100
et_insol$relation <- rownames(et_insol)
et_insol$variable <- "sunhours"

et_tmax <- as.data.frame(tapply(df_ever_turn$corr_leaf_turnover_temp, list(df_ever_turn$tmax_phase), length)) #et = ever dorm
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
  scale_y_continuous(limits = c(-72,72),
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

dd_precip <- as.data.frame(tapply(df_dec_dorm$corr_leaf_dormancy_precip, list(df_dec_dorm$precip_phase), length)) #dd = dec dorm
counts_dec_dorm <- length(df_dec_dorm$species_full)
colnames(dd_precip) <- "value"
dd_precip$value <- dd_precip$value / counts_dec_dorm *100
dd_precip$relation <- rownames(dd_precip)
dd_precip$variable <- "precipitation"

dd_insol <- as.data.frame(tapply(df_dec_dorm$corr_leaf_dormancy_sun, list(df_dec_dorm$insol_phase), length)) #dd = dec dorm
counts_dec_dorm <- length(df_dec_dorm$species_full)
colnames(dd_insol) <- "value"
dd_insol$value <- dd_insol$value / counts_dec_dorm *100
dd_insol$relation <- rownames(dd_insol)
dd_insol$variable <- "sunhours"

dd_tmax <- as.data.frame(tapply(df_dec_dorm$corr_leaf_dormancy_temp, list(df_dec_dorm$tmax_phase), length)) #dd = dec dorm
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
  scale_y_continuous(limits = c(-72,72),
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

dt_precip <- as.data.frame(tapply(df_dec_turn$corr_leaf_turnover_precip, list(df_dec_turn$precip_phase), length)) #dt = dec turn
counts_dec_turn <- length(df_dec_turn$species_full)
colnames(dt_precip) <- "value"
dt_precip$value <- dt_precip$value / counts_dec_turn *100
dt_precip$relation <- rownames(dt_precip)
dt_precip$variable <- "precipitation"

dt_insol <- as.data.frame(tapply(df_dec_turn$corr_leaf_turnover_sun, list(df_dec_turn$insol_phase), length)) #dt = dec turn
counts_dec_turn <- length(df_dec_turn$species_full)
colnames(dt_insol) <- "value"
dt_insol$value <- dt_insol$value / counts_dec_turn *100
dt_insol$relation <- rownames(dt_insol)
dt_insol$variable <- "sunhours"

dt_tmax <- as.data.frame(tapply(df_dec_turn$corr_leaf_turnover_temp, list(df_dec_turn$tmax_phase), length)) #dt = dec turn
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
  scale_y_continuous(limits = c(-72,72),
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


p_all <- grid.arrange(arrangeGrob(p_ed, p_et, heights = c(1,0.7),
                                  left = textGrob("Evergreen", gp=gpar(fontsize=12), rot = 90, hjust = 0.7)), #, hjust = 0.05, vjust = 2
                      arrangeGrob(p_dd, p_dt, heights = c(0.5,1),
                                  left = textGrob("Deciduous", gp=gpar(fontsize=12), rot = 90, hjust = 0.02)),
                      ncol = 1,
                      heights = c(0.87,1))

# pdf("~/Desktop/figure3_corr_phase2.pdf",6,4) # 5,10)
# plot(p_all)
# dev.off()





