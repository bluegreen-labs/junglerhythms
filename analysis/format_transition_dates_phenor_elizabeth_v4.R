#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
graphics.off()
#----- load required packages --------------------------------------------------#
library(zoo)
library(tidyverse)
library(circular)
library(ggplot2)
library(gridExtra)
library(scales)
#----- source required functions -----------------------------------------------#
source("R/event_length.R")
source("R/timeline_gap_fill.R")
#-------------------------------------------------------------------------------#


#----------------------------------------------------------------------
#-------- input that you can change   ---------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------
minimum_events = 1    # species will only be included in this analysis
# when they have a minimum set of events
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#--------   Phenology data - species correction Meise   ---------------
#----------------------------------------------------------------------
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
df$join_id <- paste0("R",df$image,"-",df$image_row)

metadata <- read.csv("data/phenology_archives_species_long_format_20200324.csv",
                     header = TRUE, sep = ",")
metadata$join_id <- paste(metadata$image,metadata$row, sep = "-")
# merge the two tables based upon the unique join ID
data <- merge(df, metadata, by = c("join_id"), all.x = TRUE)
data$species_full <- paste(data$genus_Meise, data$species_Meise)

# remove column id.x and rename id.y to id (--> in id.y, empty ids are renamed to EK1, EK2, etc...)
data$id <- as.character(data$id.y)
data = data[,!(names(data) %in% c("id.x","id.y"))]

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

# #----------------------------------------------------------------------
# #-------- get the species list ----------------------------------------
# #----------------------------------------------------------------------
overview <- read.csv("data/species_meta_data_phase2.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)

overview <- overview %>%
  select(species_full,
         deciduousness,
         basal_area_site,
         sd_intrasp_onset_leaf_dormancy_weeks,
         sd_intrasp_onset_leaf_turnover_weeks)

#----------------------------------------------------------------------
#-- for this species list at ID level: --------------------------------
#-- get full range timelines with 2-year gap filled -------------------
#----------------------------------------------------------------------
species_list <- overview$species_full
timelines_dorm <- two_year_gaps(data = data,
                                species_name = species_list,
                                pheno = "leaf_dormancy")
timelines_turn <- two_year_gaps(data = data,
                                species_name = species_list,
                                pheno = "leaf_turnover")
data_timeline <- rbind(timelines_dorm, timelines_turn)
#----------------------------------------------------------------------
rm(timelines_dorm, timelines_turn)
#----------------------------------------------------------------------
data_timeline <- merge(data_timeline, overview, by = "species_full", all.x = TRUE)

dec_timeline <- data_timeline %>%
  filter(grepl("deciduous",deciduousness))
dec_sp <- unique(dec_timeline$species_full)
ever_timeline <- data_timeline %>%
  filter(grepl("evergreen",deciduousness))
ever_sp <- unique(ever_timeline$species_full)
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------- This is onset of event
#---------- turnover --------------------------------------
#---------- DECIDUOUS
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# use function event_length to get the starting week of an event
transition_dates_turnover_dec <- event_length(data = dec_timeline,
                                              species_name = dec_sp,
                                              pheno = "leaf_turnover")

# weeks to degrees
transition_dates_turnover_dec <- transition_dates_turnover_dec %>%
  filter(!is.na(year_start)) %>%
  mutate(degree = (week_start-1) * 360/48) %>%# -1 so that week 1 in janurari is degree 0
  mutate(duration = (phenophase_length) * 360/48)


# summary of onset dates per species
onset_turnover_dec <- transition_dates_turnover_dec %>%
  group_by(species_full) %>%
  dplyr::summarise(
    mean_degree = mean.circular(
      circular(degree, units = "degrees")),
    median_degree = median.circular(
      circular(degree, units = "degrees")),
    lower = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees"))$mu.ci[1],
    upper = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees"))$mu.ci[2],
    nr_events = length(week_start),
    duration = mean(duration))

# timing flushing added to match columns dormancy
onset_turnover_dec$flushing_degree <- NA

# merge with census and remove rows (species) not included in the Yangambi mixed forest census
onset_turnover_dec <- merge(onset_turnover_dec, overview, by = "species_full", all.x = TRUE)

# # sort by BA
# # set y value by sorted BA number
# onset_turnover_dec <- onset_turnover_dec %>%
#   arrange(basal_area_site) %>%
#   mutate(y_value = (1:length(species_full)))

# sort by intra-species variability
# set y value by sorted SD number
# first give year-round uncertainty to NA
onset_turnover_dec$sd_intrasp_onset_leaf_turnover_weeks <- ifelse(is.na(onset_turnover_dec$sd_intrasp_onset_leaf_turnover_weeks),52, onset_turnover_dec$sd_intrasp_onset_leaf_turnover_weeks)
onset_turnover_dec <- onset_turnover_dec %>%
  arrange(desc(sd_intrasp_onset_leaf_turnover_weeks)) %>%
  mutate(y_value = (1:length(species_full)))


# rescaling between 0 and 360 degrees
onset_turnover_dec$mean_rescaled <- ifelse(onset_turnover_dec$mean_degree < 0, onset_turnover_dec$mean_degree +360, onset_turnover_dec$mean_degree)
onset_turnover_dec$median_rescaled <- ifelse(onset_turnover_dec$median_degree < 0, onset_turnover_dec$median_degree +360, onset_turnover_dec$median_degree)
onset_turnover_dec$upper_rescaled <- ifelse (onset_turnover_dec$upper < 0, onset_turnover_dec$upper + 360, onset_turnover_dec$upper)
onset_turnover_dec$lower_rescaled <- ifelse (onset_turnover_dec$lower < 0, onset_turnover_dec$lower + 360, onset_turnover_dec$lower)
onset_turnover_dec$flushing_rescaled <- NA

# in order to plot the CI if CI goes over zero, two segments needed: lower - 360 and 0 - upper
onset_turnover_dec_CItwosegments <- onset_turnover_dec %>%
  filter(mean_rescaled < lower_rescaled | mean_rescaled > upper_rescaled)
onset_turnover_dec_CItwosegments$seg1_xend <- 360  # end of segment 1 --> lower:360
onset_turnover_dec_CItwosegments$seg2_x <- 0       # beginning of segment 2 --> 0:upper
# for CI not going over zero, the segment can just be lower - upper
onset_turnover_dec_CIonesegments <- onset_turnover_dec %>%
  filter(mean_rescaled >= lower_rescaled & mean_rescaled <= upper_rescaled)

onset_turnover_dec_CIonesegments$lower_rescaled <- ifelse(onset_turnover_dec_CIonesegments$sd_intrasp_onset_leaf_turnover_weeks == 52, 0, onset_turnover_dec_CIonesegments$lower_rescaled)
onset_turnover_dec_CIonesegments$upper_rescaled <- ifelse(onset_turnover_dec_CIonesegments$sd_intrasp_onset_leaf_turnover_weeks == 52, 360, onset_turnover_dec_CIonesegments$upper_rescaled)
onset_turnover_dec_CItwosegments$lower_rescaled <- ifelse(onset_turnover_dec_CItwosegments$sd_intrasp_onset_leaf_turnover_weeks == 52, 0, onset_turnover_dec_CItwosegments$lower_rescaled)
onset_turnover_dec_CItwosegments$upper_rescaled <- ifelse(onset_turnover_dec_CItwosegments$sd_intrasp_onset_leaf_turnover_weeks == 52, 0, onset_turnover_dec_CItwosegments$upper_rescaled)
#------------------------------------------------------------------------



#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------- This is onset of event
#---------- DORMANCY --------------------------------------
#---------- DECIDUOUS
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
transition_dates_dormancy_dec <- event_length(data = dec_timeline,
                                              species_name = dec_sp,
                                              pheno = "leaf_dormancy")
# weeks to degrees
transition_dates_dormancy_dec <- transition_dates_dormancy_dec %>%
  filter(!is.na(year_start)) %>%
  mutate(degree = (week_start-1) * 360/48) %>%# -1 so that week 1 in janurari is degree 0
  mutate(duration = (phenophase_length) * 360/48)

# summary of onset dates per species
onset_dormancy_dec <- transition_dates_dormancy_dec %>%
  group_by(species_full) %>%
  dplyr::summarise(
    mean_degree = mean.circular(
      circular(degree, units = "degrees")),
    median_degree = median.circular(
      circular(degree, units = "degrees")),
    lower = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees"))$mu.ci[1],
    upper = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees"))$mu.ci[2],
    nr_events = length(week_start),
    duration = mean(duration))
# timing flushing = degree + duration + (7.5 = 1 week (1*360/48))
onset_dormancy_dec$flushing_degree <- onset_dormancy_dec$median_degree + onset_dormancy_dec$duration + 7.5

# merge with census and remove rows (species) not included in the Yangambi mixed forest census
onset_dormancy_dec <- merge(onset_dormancy_dec, overview, by = "species_full", all.x = TRUE)

# # sort by BA
# # set y value by sorted BA number
# onset_dormancy_dec <- onset_dormancy_dec %>%
#   arrange(basal_area_site) %>%
#   mutate(y_value = (1:length(species_full) + max(onset_turnover_dec$y_value) + 10))

# sort by intra-species variability
# set y value by sorted SD number
# first give year-round uncertainty to NA
onset_dormancy_dec$sd_intrasp_onset_leaf_dormancy_weeks <- ifelse(is.na(onset_dormancy_dec$sd_intrasp_onset_leaf_dormancy_weeks),52, onset_dormancy_dec$sd_intrasp_onset_leaf_dormancy_weeks)
onset_dormancy_dec <- onset_dormancy_dec %>%
  arrange(desc(sd_intrasp_onset_leaf_dormancy_weeks)) %>%
  mutate(y_value = (1:length(species_full) + max(onset_turnover_dec$y_value) + 10))

# rescaling between 0 and 360 degrees
onset_dormancy_dec$mean_rescaled <- ifelse(onset_dormancy_dec$mean_degree < 0, onset_dormancy_dec$mean_degree +360, onset_dormancy_dec$mean_degree)
onset_dormancy_dec$median_rescaled <- ifelse(onset_dormancy_dec$median_degree < 0, onset_dormancy_dec$median_degree +360, onset_dormancy_dec$median_degree)
onset_dormancy_dec$upper_rescaled <- ifelse (onset_dormancy_dec$upper < 0, onset_dormancy_dec$upper + 360, onset_dormancy_dec$upper)
onset_dormancy_dec$lower_rescaled <- ifelse (onset_dormancy_dec$lower < 0, onset_dormancy_dec$lower + 360, onset_dormancy_dec$lower)
onset_dormancy_dec$flushing_rescaled <- ifelse(onset_dormancy_dec$flushing_degree < 0, onset_dormancy_dec$flushing_degree +360, onset_dormancy_dec$flushing_degree)

# in order to plot the CI if CI goes over zero, two segments needed: lower - 360 and 0 - upper
onset_dormancy_dec_CItwosegments <- onset_dormancy_dec %>%
  filter(mean_rescaled < lower_rescaled | mean_rescaled > upper_rescaled)
onset_dormancy_dec_CItwosegments$seg1_xend <- 360  # end of segment 1 --> lower:360
onset_dormancy_dec_CItwosegments$seg2_x <- 0       # beginning of segment 2 --> 0:upper
# for CI not going over zero, the segment can just be lower - upper
onset_dormancy_dec_CIonesegments <- onset_dormancy_dec %>%
  filter(mean_rescaled >= lower_rescaled & mean_rescaled <= upper_rescaled)

onset_dormancy_dec_CIonesegments$lower_rescaled <- ifelse(onset_dormancy_dec_CIonesegments$sd_intrasp_onset_leaf_dormancy_weeks == 52, 0, onset_dormancy_dec_CIonesegments$lower_rescaled)
onset_dormancy_dec_CIonesegments$upper_rescaled <- ifelse(onset_dormancy_dec_CIonesegments$sd_intrasp_onset_leaf_dormancy_weeks == 52, 360, onset_dormancy_dec_CIonesegments$upper_rescaled)
onset_dormancy_dec_CItwosegments$lower_rescaled <- ifelse(onset_dormancy_dec_CItwosegments$sd_intrasp_onset_leaf_dormancy_weeks == 52, 0, onset_dormancy_dec_CItwosegments$lower_rescaled)
onset_dormancy_dec_CItwosegments$upper_rescaled <- ifelse(onset_dormancy_dec_CItwosegments$sd_intrasp_onset_leaf_dormancy_weeks == 52, 0, onset_dormancy_dec_CItwosegments$upper_rescaled)
#------------------------------------------------------------------------


# merge tunover and dormancy for deciduous species
onset_turnover_dec$phase <- "turnover"
onset_dormancy_dec$phase <- "dormancy"
onset_turnover_dec$nr_events <- as.numeric(onset_turnover_dec$nr_events)
onset_turnover_dec$minimum_events <- ifelse(onset_turnover_dec$nr_events > 4, "turn_over","turn_under")
onset_dormancy_dec$nr_events <- as.numeric(onset_dormancy_dec$nr_events)
onset_dormancy_dec$minimum_events <- ifelse(onset_dormancy_dec$nr_events > 4, "dorm_over","dorm_under")
onset_dec <- rbind(onset_turnover_dec, onset_dormancy_dec)


onset_turnover_dec_CIonesegments$phase <- "turnover"
onset_dormancy_dec_CIonesegments$phase <- "dormancy"
onset_dec_CIonesegments <- rbind(onset_turnover_dec_CIonesegments, onset_dormancy_dec_CIonesegments)
onset_turnover_dec_CItwosegments$phase <- "turnover"
onset_dormancy_dec_CItwosegments$phase <- "dormancy"
onset_dec_CItwosegments <- rbind(onset_turnover_dec_CItwosegments, onset_dormancy_dec_CItwosegments)



#
onset_dec_allphases <- onset_dec[,(names(onset_dec) %in% c("species_full",
                                                           "y_value",
                                                           "median_rescaled",
                                                           "phase",
                                                           "minimum_events"))]
onset_dec_flushing <- onset_dec[,(names(onset_dec) %in% c("species_full",
                                                          "y_value",
                                                          "flushing_rescaled"))]
colnames(onset_dec_flushing)[3] <- "median_rescaled"
onset_dec_flushing <- onset_dec_flushing[!(is.na(onset_dec_flushing$median_rescaled)),]
onset_dec_flushing$phase <- "flushing"
onset_dec_flushing$minimum_events <- "flush_all"

onset_dec_allphases <- rbind(onset_dec_allphases, onset_dec_flushing)





#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------- This is onset of event
#---------- turnover --------------------------------------
#---------- EVERGREEN
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# use function event_length to get the starting week of an event
transition_dates_turnover_ever <- event_length(data = ever_timeline,
                                               species_name = ever_sp,
                                               pheno = "leaf_turnover")

# weeks to degrees
transition_dates_turnover_ever <- transition_dates_turnover_ever %>%
  filter(!is.na(year_start)) %>%
  mutate(degree = (week_start-1) * 360/48) %>%# -1 so that week 1 in janurari is degree 0
  mutate(duration = (phenophase_length) * 360/48)

# summary of onset dates per species
onset_turnover_ever <- transition_dates_turnover_ever %>%
  group_by(species_full) %>%
  dplyr::summarise(
    mean_degree = mean.circular(
      circular(degree, units = "degrees")),
    median_degree = median.circular(
      circular(degree, units = "degrees")),
    lower = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees"))$mu.ci[1],
    upper = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees"))$mu.ci[2],
    nr_events = length(week_start),
    duration = mean(duration))

# timing flushing to link to column names dormancy
onset_turnover_ever$flushing_degree <- NA

# merge with census and remove rows (species) not included in the Yangambi mixed forest census
onset_turnover_ever <- merge(onset_turnover_ever, overview, by = "species_full", all.x = TRUE)

# # sort by BA
# # set y value by sorted BA number
# onset_turnover_ever <- onset_turnover_ever %>%
#   arrange(basal_area_site) %>%
#   mutate(y_value = (1:length(species_full)))

# sort by intra-species variability
# set y value by sorted SD number
# first give year-round uncertainty to NA
onset_turnover_ever$sd_intrasp_onset_leaf_turnover_weeks <- ifelse(is.na(onset_turnover_ever$sd_intrasp_onset_leaf_turnover_weeks),52, onset_turnover_ever$sd_intrasp_onset_leaf_turnover_weeks)
onset_turnover_ever <- onset_turnover_ever %>%
  arrange(desc(sd_intrasp_onset_leaf_turnover_weeks)) %>%
  mutate(y_value = (1:length(species_full)))

# rescaling between 0 and 360 degrees
onset_turnover_ever$mean_rescaled <- ifelse(onset_turnover_ever$mean_degree < 0, onset_turnover_ever$mean_degree +360, onset_turnover_ever$mean_degree)
onset_turnover_ever$median_rescaled <- ifelse(onset_turnover_ever$median_degree < 0, onset_turnover_ever$median_degree +360, onset_turnover_ever$median_degree)
onset_turnover_ever$upper_rescaled <- ifelse (onset_turnover_ever$upper < 0, onset_turnover_ever$upper + 360, onset_turnover_ever$upper)
onset_turnover_ever$lower_rescaled <- ifelse (onset_turnover_ever$lower < 0, onset_turnover_ever$lower + 360, onset_turnover_ever$lower)
onset_turnover_ever$flushing_rescaled <- NA

# in order to plot the CI if CI goes over zero, two segments needed: lower - 360 and 0 - upper
onset_turnover_ever_CItwosegments <- onset_turnover_ever %>%
  filter(mean_rescaled < lower_rescaled | mean_rescaled > upper_rescaled)
onset_turnover_ever_CItwosegments$seg1_xend <- 360  # end of segment 1 --> lower:360
onset_turnover_ever_CItwosegments$seg2_x <- 0       # beginning of segment 2 --> 0:upper
# for CI not going over zero, the segment can just be lower - upper
onset_turnover_ever_CIonesegments <- onset_turnover_ever %>%
  filter(mean_rescaled >= lower_rescaled & mean_rescaled <= upper_rescaled)

onset_turnover_ever_CIonesegments$lower_rescaled <- ifelse(onset_turnover_ever_CIonesegments$sd_intrasp_onset_leaf_turnover_weeks == 52, 0, onset_turnover_ever_CIonesegments$lower_rescaled)
onset_turnover_ever_CIonesegments$upper_rescaled <- ifelse(onset_turnover_ever_CIonesegments$sd_intrasp_onset_leaf_turnover_weeks == 52, 360, onset_turnover_ever_CIonesegments$upper_rescaled)
onset_turnover_ever_CItwosegments$lower_rescaled <- ifelse(onset_turnover_ever_CItwosegments$sd_intrasp_onset_leaf_turnover_weeks == 52, 0, onset_turnover_ever_CItwosegments$lower_rescaled)
onset_turnover_ever_CItwosegments$upper_rescaled <- ifelse(onset_turnover_ever_CItwosegments$sd_intrasp_onset_leaf_turnover_weeks == 52, 0, onset_turnover_ever_CItwosegments$upper_rescaled)


#------------------------------------------------------------------------



#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------- This is onset of event
#---------- DORMANCY --------------------------------------
#---------- EVERGREEN
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# use function event_length to get the starting week of an event
transition_dates_dormancy_ever <- event_length(data = ever_timeline,
                                               species_name = ever_sp,
                                               pheno = "leaf_dormancy")

# weeks to degrees
transition_dates_dormancy_ever <- transition_dates_dormancy_ever %>%
  filter(!is.na(year_start)) %>%
  mutate(degree = (week_start-1) * 360/48) %>%# -1 so that week 1 in janurari is degree 0
  mutate(duration = (phenophase_length) * 360/48)

# summary of onset dates per species
onset_dormancy_ever <- transition_dates_dormancy_ever %>%
  group_by(species_full) %>%
  dplyr::summarise(
    mean_degree = mean.circular(
      circular(degree, units = "degrees")),
    median_degree = median.circular(
      circular(degree, units = "degrees")),
    lower = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees"))$mu.ci[1],
    upper = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees"))$mu.ci[2],
    nr_events = length(week_start),
    duration = mean(duration))
# timing flushing = degree + duration + (7.5 = 1 week (1*360/48))
onset_dormancy_ever$flushing_degree <- onset_dormancy_ever$median_degree + onset_dormancy_ever$duration + 7.5

# merge with census and remove rows (species) not included in the Yangambi mixed forest census
onset_dormancy_ever <- merge(onset_dormancy_ever, overview, by = "species_full", all.x = TRUE)

# # sort by BA
# # set y value by sorted BA number
# onset_dormancy_ever <- onset_dormancy_ever %>%
#   arrange(basal_area_site) %>%
#   mutate(y_value = (1:length(species_full) + max(onset_turnover_ever$y_value) + 8))

# sort by intra-species variability
# set y value by sorted SD number
# first give year-round uncertainty to NA
onset_dormancy_ever$sd_intrasp_onset_leaf_dormancy_weeks <- ifelse(is.na(onset_dormancy_ever$sd_intrasp_onset_leaf_dormancy_weeks),52, onset_dormancy_ever$sd_intrasp_onset_leaf_dormancy_weeks)
onset_dormancy_ever <- onset_dormancy_ever %>%
  arrange(desc(sd_intrasp_onset_leaf_dormancy_weeks)) %>%
  mutate(y_value = (1:length(species_full) + max(onset_turnover_ever$y_value) + 8))

# rescaling between 0 and 360 degrees
onset_dormancy_ever$mean_rescaled <- ifelse(onset_dormancy_ever$mean_degree < 0, onset_dormancy_ever$mean_degree +360, onset_dormancy_ever$mean_degree)
onset_dormancy_ever$median_rescaled <- ifelse(onset_dormancy_ever$median_degree < 0, onset_dormancy_ever$median_degree +360, onset_dormancy_ever$median_degree)
onset_dormancy_ever$upper_rescaled <- ifelse (onset_dormancy_ever$upper < 0, onset_dormancy_ever$upper + 360, onset_dormancy_ever$upper)
onset_dormancy_ever$lower_rescaled <- ifelse (onset_dormancy_ever$lower < 0, onset_dormancy_ever$lower + 360, onset_dormancy_ever$lower)
onset_dormancy_ever$flushing_rescaled <- ifelse(onset_dormancy_ever$flushing_degree < 0, onset_dormancy_ever$flushing_degree +360, onset_dormancy_ever$flushing_degree)

# in order to plot the CI if CI goes over zero, two segments needed: lower - 360 and 0 - upper
onset_dormancy_ever_CItwosegments <- onset_dormancy_ever %>%
  filter(mean_rescaled < lower_rescaled | mean_rescaled > upper_rescaled)
onset_dormancy_ever_CItwosegments$seg1_xend <- 360  # end of segment 1 --> lower:360
onset_dormancy_ever_CItwosegments$seg2_x <- 0       # beginning of segment 2 --> 0:upper
# for CI not going over zero, the segment can just be lower - upper
onset_dormancy_ever_CIonesegments <- onset_dormancy_ever %>%
  filter(mean_rescaled >= lower_rescaled & mean_rescaled <= upper_rescaled)

onset_dormancy_ever_CIonesegments$lower_rescaled <- ifelse(onset_dormancy_ever_CIonesegments$sd_intrasp_onset_leaf_dormancy_weeks == 52, 0, onset_dormancy_ever_CIonesegments$lower_rescaled)
onset_dormancy_ever_CIonesegments$upper_rescaled <- ifelse(onset_dormancy_ever_CIonesegments$sd_intrasp_onset_leaf_dormancy_weeks == 52, 360, onset_dormancy_ever_CIonesegments$upper_rescaled)
onset_dormancy_ever_CItwosegments$lower_rescaled <- ifelse(onset_dormancy_ever_CItwosegments$sd_intrasp_onset_leaf_dormancy_weeks == 52, 0, onset_dormancy_ever_CItwosegments$lower_rescaled)
onset_dormancy_ever_CItwosegments$upper_rescaled <- ifelse(onset_dormancy_ever_CItwosegments$sd_intrasp_onset_leaf_dormancy_weeks == 52, 0, onset_dormancy_ever_CItwosegments$upper_rescaled)

#------------------------------------------------------------------------


# merge tunover and dormancy for deciduous species
onset_turnover_ever$phase <- "turnover"
onset_dormancy_ever$phase <- "dormancy"
onset_turnover_ever$nr_events <- as.numeric(onset_turnover_ever$nr_events)
onset_turnover_ever$minimum_events <- ifelse(onset_turnover_ever$nr_events > 4, "turn_over","turn_under")
onset_dormancy_ever$nr_events <- as.numeric(onset_dormancy_ever$nr_events)
onset_dormancy_ever$minimum_events <- ifelse(onset_dormancy_ever$nr_events > 4, "dorm_over","dorm_under")
onset_ever <- rbind(onset_turnover_ever, onset_dormancy_ever)


onset_turnover_ever_CIonesegments$phase <- "turnover"
onset_dormancy_ever_CIonesegments$phase <- "dormancy"
onset_ever_CIonesegments <- rbind(onset_turnover_ever_CIonesegments, onset_dormancy_ever_CIonesegments)
onset_turnover_ever_CItwosegments$phase <- "turnover"
onset_dormancy_ever_CItwosegments$phase <- "dormancy"
onset_ever_CItwosegments <- rbind(onset_turnover_ever_CItwosegments, onset_dormancy_ever_CItwosegments)

#
onset_ever_allphases <- onset_ever[,(names(onset_ever) %in% c("y_value",
                                                              "median_rescaled",
                                                              "phase",
                                                              "minimum_events"))]
onset_ever_flushing <- onset_ever[,(names(onset_ever) %in% c("y_value",
                                                             "flushing_rescaled"))]
colnames(onset_ever_flushing)[2] <- "median_rescaled"
onset_ever_flushing <- onset_ever_flushing[!(is.na(onset_ever_flushing$median_rescaled)),]
onset_ever_flushing$phase <- "flushing"
onset_ever_flushing$minimum_events <- "flush_all"

onset_ever_allphases <- rbind(onset_ever_allphases, onset_ever_flushing)

#------------------------------------------------------------------------
# FIGURE DECIDUOUS
#------------------------------------------------------------------------
p_onset_dec <- ggplot() +
  annotate("rect", xmin = 330, xmax = 360, ymin = -20, ymax = max(onset_dec$y_value) , alpha = .2) + #Dec
  annotate("rect", xmin = 0, xmax = 60, ymin = -20, ymax = max(onset_dec$y_value) , alpha = .2) + # jan - feb
  annotate("rect", xmin = 150, xmax = 210, ymin = -20, ymax = max(onset_dec$y_value) , alpha = .2) + # jun-jul

  annotate("text", x = 0, y = -8, label = "LD", col = "grey50") +
  annotate("text", x = 90, y = -8, label = "SW", col = "grey50") +
  annotate("text", x = 180, y = -8, label = "SD", col = "grey50") +
  annotate("text", x = 270, y = -8, label = "LW", col = "grey50") +

  geom_segment(data = onset_dec_CIonesegments,
               aes(x = lower_rescaled,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = onset_dec_CItwosegments,
               aes(x = lower_rescaled,
                   xend = seg1_xend,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = onset_dec_CItwosegments,
               aes(x = seg2_x,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_point(data = onset_dec_allphases,
             aes(x = median_rescaled, y = y_value, shape = phase)) + #, col = minimum_events)) +
  # aes(x = median_rescaled, y = y_value, shape = minimum_events)) +
  scale_shape_manual(values = c(19,4,1),
                     name = "onset of:",
                     labels = c("dormancy", "flushing","turnover")) +
  # scale_color_manual(values = c("turn_over" = "#d8b365", "turn_under" = "black",
  #                               "dorm_over" = "#8c510a", "dorm_under" = "black",
  #                               "flush_all" = "black")) +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0,359,30),
                     labels = month.abb) +
  scale_y_continuous(limits = c(-20,max(onset_dec$y_value))) +  #staring y-axis at '-20' so that points are squeezed in the center
  # minor_breaks = seq(1, 48, 1), breaks = seq(1, 48, 1)) +
  coord_polar() +
  labs(x="",
       y="",
       title = paste("(b) deciduous")) +
  theme(panel.grid.major.x = element_line(colour = "grey75",
                                          size = 0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),# element_line(colour = "grey90", size = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11),
        strip.text = element_text(face = "italic", size = 13),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(vjust = -13,
                                  hjust = 0.05,
                                  size = 13),
        legend.position = c(0.5,-0.01),#"bottom",
        # legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.key = element_rect(fill = "white"),
        plot.margin=unit(c(0,0,0,0),"cm")) +
  guides(shape = guide_legend(title.position = "left",
                              byrow = TRUE,
                              nrow = 1))
p_onset_dec
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# FIGURE EVERGREEN
#------------------------------------------------------------------------
p_onset_ever <- ggplot() +
  annotate("rect", xmin = 330, xmax = 360, ymin = -20, ymax = max(onset_ever$y_value) , alpha = .2) + #Dec
  annotate("rect", xmin = 0, xmax = 60, ymin = -20, ymax = max(onset_ever$y_value) , alpha = .2) + # jan - feb
  annotate("rect", xmin = 150, xmax = 210, ymin = -20, ymax = max(onset_ever$y_value) , alpha = .2) + # jun-jul

  annotate("text", x = 0, y = -8, label = "LD", col = "grey50") +
  annotate("text", x = 90, y = -8, label = "SW", col = "grey50") +
  annotate("text", x = 180, y = -8, label = "SD", col = "grey50") +
  annotate("text", x = 270, y = -8, label = "LW", col = "grey50") +

  geom_segment(data = onset_ever_CIonesegments,
               aes(x = lower_rescaled,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = onset_ever_CItwosegments,
               aes(x = lower_rescaled,
                   xend = seg1_xend,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = onset_ever_CItwosegments,
               aes(x = seg2_x,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_point(data = onset_ever_allphases,
             aes(x = median_rescaled, y = y_value, shape = phase)) +#, col = minimum_events)) +
             # aes(x = median_rescaled, y = y_value, shape = minimum_events)) +
  scale_shape_manual(values = c(19,4,1),
                     name = "onset of:",
                     labels = c("dormancy", "flushing","turnover")) +
  # scale_shape_manual(values = c("turn_over" = 15, "turn_under" = 0,
  #                               "dorm_over" = 19, "dorm_under" = 1,
  #                               "flush_all" = 4)) +
  # scale_color_manual(values = c("turn_over" = "#8c510a", "turn_under" = "black",
  #                               "dorm_over" = "#8c510a", "dorm_under" = "black",
  #                               "flush_all" = "black")) +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0,359,30),
                     labels = month.abb) +
  scale_y_continuous(limits = c(-20,max(onset_ever$y_value))) +  #staring y-axis at '-20' so that points are squeezed in the center
  # minor_breaks = seq(1, 48, 1), breaks = seq(1, 48, 1)) +
  coord_polar() +
  labs(x="",
       y="",
       title = paste("(a) evergreen")) +
  theme(panel.grid.major.x = element_line(colour = "grey75",
                                          size = 0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11),
        strip.text = element_text(face = "italic", size = 13),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(vjust = -13,
                                  hjust = 0.05,
                                  size = 13),
        legend.position = "none",
        plot.margin=unit(c(0,0,-1,0),"cm")
  )
p_onset_ever
#------------------------------------------------------------------------

# p_onset_ever <- ggplot_gtable(ggplot_build(p_onset_ever))
# p_onset_dec <- ggplot_gtable(ggplot_build(p_onset_dec))
# p_onset_ever$widths <-p_onset_dec$widths

p_onset <- grid.arrange(p_onset_ever, p_onset_dec, heights = c(1,1))


# pdf("~/Desktop/figure2_test3.pdf",4.95,10)
# plot(p_onset)
# dev.off()

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------
# FIGURE DECIDUOUS
#------------------------------------------------------------------------
dec_dorm_CIone <- onset_dec_CIonesegments %>%
  filter(phase %in% c("dormancy","flushing"))
dec_dorm_CItwo <- onset_dec_CItwosegments %>%
  filter(phase %in% c("dormancy","flushing"))
dec_dorm <- onset_dec_allphases %>%
  filter(phase %in% c("dormancy","flushing"))

p_dec_dorm <- ggplot() +
  annotate("rect", xmin = 330, xmax = 360, ymin = min(dec_dorm$y_value)-10, ymax = max(dec_dorm$y_value) , alpha = .2) + #Dec
  annotate("rect", xmin = 0, xmax = 60, ymin = min(dec_dorm$y_value)-10, ymax = max(dec_dorm$y_value) , alpha = .2) + # jan - feb
  annotate("rect", xmin = 150, xmax = 210, ymin = min(dec_dorm$y_value)-10, ymax = max(dec_dorm$y_value) , alpha = .2) + # jun-jul

  geom_segment(data = dec_dorm_CIone,
               aes(x = lower_rescaled,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = dec_dorm_CItwo,
               aes(x = lower_rescaled,
                   xend = seg1_xend,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = dec_dorm_CItwo,
               aes(x = seg2_x,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_point(data = dec_dorm,
             aes(x = median_rescaled, y = y_value, shape = phase)) +
  scale_shape_manual(values = c(19,4,1),
                     name = "onset of:",
                     labels = c("dormancy", "flushing","turnover")) +

  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0,359,30),
                     labels = month.abb) +
  scale_y_continuous(limits = c(min(dec_dorm$y_value)-10,max(dec_dorm$y_value))) +  #staring y-axis at '-20' so that points are squeezed in the center
  # minor_breaks = seq(1, 48, 1), breaks = seq(1, 48, 1)) +
  coord_polar() +
  labs(x="",
       y="",
       title = paste("(b) deciduous - dorm")) +
  theme(panel.grid.major.x = element_line(colour = "grey75",
                                          size = 0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),# element_line(colour = "grey90", size = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11),
        strip.text = element_text(face = "italic", size = 13),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(size = 13),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.key = element_rect(fill = "white"),
        plot.margin=unit(c(0,0,0,0),"cm"))
p_dec_dorm
#------------------------------------------------------------------------
dec_turn_CIone <- onset_dec_CIonesegments %>%
  filter(phase == "turnover")
dec_turn_CItwo <- onset_dec_CItwosegments %>%
  filter(phase == "turnover")
dec_turn <- onset_dec_allphases %>%
  filter(phase == "turnover")

p_dec_turn <- ggplot() +
  annotate("rect", xmin = 330, xmax = 360, ymin = min(dec_turn$y_value)-10, ymax = max(dec_turn$y_value) , alpha = .2) + #Dec
  annotate("rect", xmin = 0, xmax = 60, ymin = min(dec_turn$y_value)-10, ymax = max(dec_turn$y_value) , alpha = .2) + # jan - feb
  annotate("rect", xmin = 150, xmax = 210, ymin = min(dec_turn$y_value)-10, ymax = max(dec_turn$y_value) , alpha = .2) + # jun-jul

  geom_segment(data = dec_turn_CIone,
               aes(x = lower_rescaled,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = dec_turn_CItwo,
               aes(x = lower_rescaled,
                   xend = seg1_xend,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = dec_turn_CItwo,
               aes(x = seg2_x,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_point(data = dec_turn,
             aes(x = median_rescaled, y = y_value, shape = phase)) +
  scale_shape_manual(values = 1) +

  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0,359,30),
                     labels = month.abb) +
  scale_y_continuous(limits = c(min(dec_turn$y_value)-10,max(dec_turn$y_value))) +  #staring y-axis at '-20' so that points are squeezed in the center
  # minor_breaks = seq(1, 48, 1), breaks = seq(1, 48, 1)) +
  coord_polar() +
  labs(x="",
       y="",
       title = paste("(b) deciduous - turn")) +
  theme(panel.grid.major.x = element_line(colour = "grey75",
                                          size = 0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),# element_line(colour = "grey90", size = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11),
        strip.text = element_text(face = "italic", size = 13),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(size = 13),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.key = element_rect(fill = "white"),
        plot.margin=unit(c(0,0,0,0),"cm"))
p_dec_turn
#---------------
#------------------------------------------------------------------------
# FIGURE EVERGREEN
#------------------------------------------------------------------------
ever_dorm_CIone <- onset_ever_CIonesegments %>%
  filter(phase %in% c("dormancy","flushing"))
ever_dorm_CItwo <- onset_ever_CItwosegments %>%
  filter(phase %in% c("dormancy","flushing"))
ever_dorm <- onset_ever_allphases %>%
  filter(phase %in% c("dormancy","flushing"))

p_ever_dorm <- ggplot() +
  annotate("rect", xmin = 330, xmax = 360, ymin = min(ever_dorm$y_value)-10, ymax = max(ever_dorm$y_value) , alpha = .2) + #ever
  annotate("rect", xmin = 0, xmax = 60, ymin = min(ever_dorm$y_value)-10, ymax = max(ever_dorm$y_value) , alpha = .2) + # jan - feb
  annotate("rect", xmin = 150, xmax = 210, ymin = min(ever_dorm$y_value)-10, ymax = max(ever_dorm$y_value) , alpha = .2) + # jun-jul

  geom_segment(data = ever_dorm_CIone,
               aes(x = lower_rescaled,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = ever_dorm_CItwo,
               aes(x = lower_rescaled,
                   xend = seg1_xend,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = ever_dorm_CItwo,
               aes(x = seg2_x,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_point(data = ever_dorm,
             aes(x = median_rescaled, y = y_value, shape = phase)) +
  scale_shape_manual(values = c(19,4,1),
                     name = "onset of:",
                     labels = c("dormancy", "flushing","turnover")) +

  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0,359,30),
                     labels = month.abb) +
  scale_y_continuous(limits = c(min(ever_dorm$y_value)-10,max(ever_dorm$y_value))) +  #staring y-axis at '-20' so that points are squeezed in the center
  # minor_breaks = seq(1, 48, 1), breaks = seq(1, 48, 1)) +
  coord_polar() +
  labs(x="",
       y="",
       title = paste("(b) evergreen - dorm")) +
  theme(panel.grid.major.x = element_line(colour = "grey75",
                                          size = 0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),# element_line(colour = "grey90", size = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11),
        strip.text = element_text(face = "italic", size = 13),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(size = 13),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.key = element_rect(fill = "white"),
        plot.margin=unit(c(0,0,0,0),"cm"))
p_ever_dorm
#------------------------------------------------------------------------
ever_turn_CIone <- onset_ever_CIonesegments %>%
  filter(phase == "turnover")
ever_turn_CItwo <- onset_ever_CItwosegments %>%
  filter(phase == "turnover")
ever_turn <- onset_ever_allphases %>%
  filter(phase == "turnover")

p_ever_turn <- ggplot() +
  annotate("rect", xmin = 330, xmax = 360, ymin = min(ever_turn$y_value)-10, ymax = max(ever_turn$y_value) , alpha = .2) + #ever
  annotate("rect", xmin = 0, xmax = 60, ymin = min(ever_turn$y_value)-10, ymax = max(ever_turn$y_value) , alpha = .2) + # jan - feb
  annotate("rect", xmin = 150, xmax = 210, ymin = min(ever_turn$y_value)-10, ymax = max(ever_turn$y_value) , alpha = .2) + # jun-jul

  geom_segment(data = ever_turn_CIone,
               aes(x = lower_rescaled,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = ever_turn_CItwo,
               aes(x = lower_rescaled,
                   xend = seg1_xend,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = ever_turn_CItwo,
               aes(x = seg2_x,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_point(data = ever_turn,
             aes(x = median_rescaled, y = y_value, shape = phase)) +
  scale_shape_manual(values = 1) +

  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0,359,30),
                     labels = month.abb) +
  scale_y_continuous(limits = c(min(ever_turn$y_value)-10,max(ever_turn$y_value))) +  #staring y-axis at '-20' so that points are squeezed in the center
  # minor_breaks = seq(1, 48, 1), breaks = seq(1, 48, 1)) +
  coord_polar() +
  labs(x="",
       y="",
       title = paste("(b) evergreen - turn")) +
  theme(panel.grid.major.x = element_line(colour = "grey75",
                                          size = 0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),# element_line(colour = "grey90", size = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11),
        strip.text = element_text(face = "italic", size = 13),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(size = 13),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.key = element_rect(fill = "white"),
        plot.margin=unit(c(0,0,0,0),"cm"))
p_ever_turn
#---------------

p_all <- grid.arrange(p_ever_dorm, p_ever_turn, p_dec_dorm, p_dec_turn, heights = c(1,1))
