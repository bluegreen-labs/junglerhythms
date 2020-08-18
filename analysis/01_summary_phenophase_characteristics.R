#----- reset your R session. ------------------------------------------
rm(list=ls())
# graphics.off()
#----- load required packages -----------------------------------------
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(circular)
#----- source required functions --------------------------------------
source("R/event_length.R")
source("R/summary_stats.R")
source("R/timeline_gap_fill.R")
#----------------------------------------------------------------------

# Suppress summarise info
# because `summarise()` ungrouping output (override with `.groups` argument)
# comes up a lot in the consol, which is annoying
options(dplyr.summarise.inform = FALSE)

#----------------------------------------------------------------------
#--------   Phenology data --------------------------------------------
#----------------------------------------------------------------------
data <- readRDS("data/jungle_rhythms_data_manuscript_leaf_repro.rds")
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#-------- census data  ------------------------------------------------
#-------- get species-specific basal area, abundance, max dbh
#----------------------------------------------------------------------
census <- read.csv("data/YGB_ForestPlotsNET_corrected_indet.csv",
                   header = TRUE,
                   sep = ",",
                   stringsAsFactors = FALSE)
census <- census %>%
  dplyr::rename("species_full" = Species)
census$species_full <- ifelse(census$species_full == "Unknown", NA, census$species_full)

# remove individuals without C1DBH4,
# these are new recruits for census2
census <- census[!(is.na(census$C1DBH4)),]
# and calculate basal_area for each individual
# C1DBH4 unit is mm -> /1000 for m
census$basal_area = pi*(census$C1DBH4/2000)^2

# only use mixed plots
# get %basal area, abundance and max dbh for each species at site level (1-ha stats)
census <- census %>%
  filter(grepl("MIX",Plot)) %>%
  group_by(species_full) %>%
  dplyr::summarise(basal_area_site = sum(basal_area)/5,
                   abundance = length(basal_area)/5, #/5, because 5 1-ha plots
                   dbh_max_mm_census = max(C1DBH4))
total_basal_area <- sum(census$basal_area_site)
census$basal_area_percentage <- census$basal_area_site/total_basal_area *100
census <- census %>%
  dplyr::select(species_full,
                basal_area_percentage,
                abundance,
                dbh_max_mm_census)
#----------------------------------------------------------------------
rm(total_basal_area)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#-- for this species list at ID level: --------------------------------
#-- get full range timelines ------------------------------------------
#----------------------------------------------------------------------
species_list <- unique(data$species_full)

timelines_dorm <- missing_year_gaps(data = data,
                                    species_name = species_list,
                                    pheno = "leaf_dormancy",
                                    gapfill_missingyears = 0)
timelines_turn <- missing_year_gaps(data = data,
                                    species_name = species_list,
                                    pheno = "leaf_turnover",
                                    gapfill_missingyears = 0)
data_timeline <- rbind(timelines_dorm, timelines_turn)
#----------------------------------------------------------------------
rm(timelines_dorm, timelines_turn)
#----------------------------------------------------------------------


#-------------------------------------------------------------------------------
# summary statistics:
# start year, end year,
# number of observation years,
# number of site years with events,
# nr of individuals
#-------------------------------------------------------------------------------
stats_LD <- overview_stats(data = data_timeline,
                           species_name = species_list,
                           pheno = "leaf_dormancy")
stats_LT <- overview_stats(data = data_timeline,
                           species_name = species_list,
                           pheno = "leaf_turnover")
output_stats <- rbind(stats_LD, stats_LT)
rm(stats_LD, stats_LT)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Average duration of an event
#-------------------------------------------------------------------------------
# use function event_length to get the starting week of each event (individual level)
event_LD <- event_length(data = data_timeline,
                         species_name = species_list,
                         pheno = "leaf_dormancy")
event_LT <- event_length(data = data_timeline,
                         species_name = species_list,
                         pheno = "leaf_turnover")
event_indiv <- rbind(event_LD, event_LT)

# mean and sd of length phenophase (species level)
# rescaling 48-week year to 52-week year
output_event_duration <- event_indiv %>%
  group_by(species_full, phenophase) %>%
  dplyr::summarise(mean_phenophase_length_weeks = mean(phenophase_length, na.rm = TRUE)/48 *52,
                   sd_phenophase_length_weeks = sd(phenophase_length, na.rm = TRUE)/48 *52,
                   total_nr_events = length(phenophase_length[!(is.na(phenophase_length))]))
#-------------------------------------------------------------------------------
rm(event_LD, event_LT)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#------------ Synchrony index among individuals across years -------------------
#-------------------------------------------------------------------------------
# calculate the synchrony between events as mean and SD of onset (= week_start)
# then reformat those taking the year - date line into consideration
transition_dates <- event_indiv %>%
  filter(!is.na(week_start)) %>%
  mutate(degree = week_start * 360/48)

onset_event <- transition_dates %>%
  group_by(species_full, phenophase) %>%
  dplyr::summarise(
    mean_degree = mean.circular(circular(degree, units = "degrees")),
    median_degree = median.circular(circular(degree, units = "degrees")),
    sd_rad = sd.circular(circular(degree, units = "degrees")), # eventhough units in degree, sd calculated in radians # https://stats.stackexchange.com/questions/185361/how-do-i-interpret-the-standard-deviation-of-a-directional-dataset
    nr_events_onset = length(week_start)) %>%
  mutate(sd_degree = deg(sd_rad))
# only sd for species with more than 3 events
onset_event$sd_degree <- ifelse(onset_event$nr_events_onset >= 3, onset_event$sd_degree, NA)
# mean: rescaling between 0 and 360 degrees
onset_event$mean_rescaled <- ifelse(onset_event$mean_degree < 0, onset_event$mean_degree +360, onset_event$mean_degree)
onset_event$median_rescaled <- ifelse(onset_event$median_degree < 0, onset_event$median_degree +360, onset_event$median_degree)
# rescaling degrees to weeks
# and only keep columns you need
output_onset_event <- onset_event %>%
  mutate(mean_intrasp_onset_weeks = mean_rescaled /360 *52,
         median_intrasp_onset_weeks = median_rescaled /360 *52,
         sd_intrasp_onset_weeks = sd_degree/360*52) %>%
  select(species_full,
         phenophase,
         mean_intrasp_onset_weeks,
         median_intrasp_onset_weeks,
         sd_intrasp_onset_weeks)
#-------------------------------------------------------------------------------
rm(onset_event, event_indiv)



# #-------------------------------------------------------------------------------
# #------------ Synchrony index between species as -------------------
# #----- average pairwise distance between median onset date of species
# #-------------------------------------------------------------------------------
# # you need to merge with census here, because else distance is calculated to all other species in full phenology dataset
# census_overview <- census[,(names(census) %in% c("species_full","basal_area_site"))]
# distance_events_LD <- merge(onset_LD, census_overview, by = "species_full", all.x = TRUE)
# distance_events_LD <- distance_events_LD[!(is.na(distance_events_LD$basal_area_site)),]
# distance_events_LD$median_rad <- rad(distance_events_LD$median_degree)
# distance_events_LD <- distance_events_LD[!(is.na(distance_events_LD$median_rad)),]
# a <- as.matrix(dist(distance_events_LD$median_rad), labels = TRUE)
# rownames(a) <- distance_events_LD[['species_full']]
# colnames(a) <- rownames(a)
# b <- deg(a)
# c <- ifelse(b > 180, 360-b, b)
# d <- as.data.frame(rowMeans(c))
# d$species_full <- rownames(d)
# colnames(d)[1] <- "mean_distance_onset_leaf_dormancy_weeks"
# d <- d %>%
#   mutate(mean_distance_onset_leaf_dormancy_weeks = mean_distance_onset_leaf_dormancy_weeks /360 *52)
# #---------- now with a mimimum of 5 events
# onset_LD_minfreq <- transition_dates_LD %>%
#   group_by(species_full) %>%
#   filter(length(week_start) >= 5) %>%    # only species with more than 5 events
#   summarise(median_degree = median.circular(circular(degree, units = "degrees")))
# distance_events_minfreq_LD <- merge(onset_LD_minfreq, census_overview, by = "species_full", all.x = TRUE)
# distance_events_minfreq_LD <- distance_events_minfreq_LD[!(is.na(distance_events_minfreq_LD$basal_area_site)),]
# distance_events_minfreq_LD$median_rad <- rad(distance_events_minfreq_LD$median_degree)
# distance_events_minfreq_LD <- distance_events_minfreq_LD[!(is.na(distance_events_minfreq_LD$median_rad)),]
# a_minfreq <- as.matrix(dist(distance_events_minfreq_LD$median_rad), labels = TRUE)
# rownames(a_minfreq) <- distance_events_minfreq_LD[['species_full']]
# colnames(a_minfreq) <- rownames(a_minfreq)
# b_minfreq <- deg(a_minfreq)
# c_minfreq <- ifelse(b_minfreq > 180, 360-b_minfreq, b_minfreq)
# d_minfreq <- as.data.frame(rowMeans(c_minfreq))
# d_minfreq$species_full <- rownames(d_minfreq)
# colnames(d_minfreq)[1] <- "mean_distance_onset_minfreq_leaf_dormancy_weeks"
# d_minfreq <- d_minfreq %>%
#   mutate(mean_distance_onset_minfreq_leaf_dormancy_weeks = mean_distance_onset_minfreq_leaf_dormancy_weeks /360 *52)
# #-----------
# onset_LD <- merge(onset_LD, d, by = "species_full", all.x = TRUE)
# onset_LD <- merge(onset_LD, d_minfreq, by = "species_full", all.x = TRUE)
# # remove columns you don't need
# onset_LD = onset_LD[,!(names(onset_LD) %in% c("mean_degree","sd_rad","sd_degree","nr_events_onset","mean_rescaled","median_degree","median_rescaled","median_rad"))]
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
#------------ Synchrony indexes within individuals across years ----------------
#-------------------------------------------------------------------------------
#-- Synchrony index within individuals across years ----------------------------
#-- will be writen out as separate file
onset_ind <- transition_dates %>%
  group_by(id, species_full, phenophase) %>%
  filter(length(week_start) >= 3) %>%    # only species with more than 3 events
  dplyr::summarise(sd_rad = sd.circular(circular(degree, units = "degrees")), # eventhough units in degree, sd calculated in radians
                   nr_events_onset = length(week_start)) %>%
  mutate(sd_degree = deg(sd_rad))
# rescaling degrees to weeks
output_synchrony_ind <- onset_ind %>%
  mutate(onset_sd_indiv_weeks = sd_degree/360*52) %>%
  select(species_full,
         phenophase,
         id,
         onset_sd_indiv_weeks,
         nr_events_onset)

#-- Synchrony index within species across years and individuals ----------------
onset_sp <- transition_dates %>%
  group_by(species_full, phenophase) %>%
  filter(length(week_start) >= 3) %>%    # only species with more than 3 events
  dplyr::summarise(sd_rad = sd.circular(circular(degree, units = "degrees"))) %>% # eventhough units in degree, sd calculated in radians
  mutate(sd_degree = deg(sd_rad))
# rescaling degrees to weeks
output_synchrony_sp <- onset_sp %>%
  mutate(onset_sd_species_weeks = sd_degree/360*52) %>%
  select(species_full,
         phenophase,
         onset_sd_species_weeks)
#-------------------------------------------------------------------------------
rm(onset_ind, onset_sp, transition_dates)
#-------------------------------------------------------------------------------


#----------------------------------------------------------------------
#-------- read in trait data from Steven Janssens ---------------------
#-------- get species-specific basal area at plot and site level ------
#----------------------------------------------------------------------
traits <- read.csv("data/Dataset_traits_African_trees.csv",
                   header = TRUE,
                   sep = ",",
                   stringsAsFactors = FALSE)
traits <- traits %>%
  dplyr::rename("mating_system" = Mating_system,
                "ecology" = Ecology,
                "height_m_literature" = Height_lit_m,
                "dbh_max_cm_literature" = Dmax_lit_cm)

# redefine deciousness of some species (because classified as 'evergreen or deciduous'; 'sometimes deciduous')
# redefine based on the archive phenology data
#----------------------------------------------------------------------
# clear patterns:
traits$deciduousness <- ifelse(traits$species_full %in% "Pericopsis elata", "deciduous*",traits$deciduousness)
traits$deciduousness <- ifelse(traits$species_full %in% "Trichilia welwitschii", "evergreen*",traits$deciduousness)
traits$deciduousness <- ifelse(traits$species_full %in% "Copaifera mildbraedii", "deciduous*",traits$deciduousness)
traits$deciduousness <- ifelse(traits$species_full %in% "Tridesmostemon omphalocarpoides", "evergreen*",traits$deciduousness)
traits$deciduousness <- ifelse(traits$species_full %in% "Omphalocarpum lecomteanum", "evergreen*",traits$deciduousness)
traits$deciduousness <- ifelse(traits$species_full %in% "Fernandoa adolfi-friderici", "deciduous*",traits$deciduousness)
traits$deciduousness <- ifelse(traits$species_full %in% "Tabernaemontana crassa", "evergreen*",traits$deciduousness)
traits$deciduousness <- ifelse(traits$species_full %in% "Trichilia tessmannii", "deciduous*",traits$deciduousness)
# not so sure, limited data:
traits$deciduousness <- ifelse(traits$species_full %in% "Trichilia gilletii", "evergreen* (?)",traits$deciduousness)
# not so sure, unclear phenological data:
traits$deciduousness <- ifelse(traits$species_full %in% "Radlkofera calodendron", "evergreen* (?)",traits$deciduousness)
traits$deciduousness <- ifelse(traits$species_full %in% "Gilletiodendron mildbraedii", "evergreen* (?)",traits$deciduousness)
## two stars, in literature found as evergreen or (sometimes) deciduous
## selected a class based on the actual data
traits$deciduousness <- ifelse(traits$species_full %in% "Celtis mildbraedii", "evergreen**",traits$deciduousness)
traits$deciduousness <- ifelse(traits$species_full %in% "Combretum lokele", "deciduous**",traits$deciduousness)
traits$deciduousness <- ifelse(traits$species_full %in% "Homalium africanum", "evergreen**",traits$deciduousness)
traits$deciduousness <- ifelse(traits$species_full %in% "Quassia silvestris", "evergreen**",traits$deciduousness)
# not so sure, unclear phenological data
traits$deciduousness <- ifelse(traits$species_full %in% "Homalium longistylum", "evergreen** (?)",traits$deciduousness)
traits$deciduousness <- ifelse(traits$species_full %in% "Irvingia gabonensis", "deciduous** (?)",traits$deciduousness)



#--------------------------------------------------------------------
#--------------------------------------------------------------------
# merge everything
#--------------------------------------------------------------------
#--------------------------------------------------------------------
final <- merge(output_stats, output_event_duration, by = c("species_full", "phenophase"), all.x = TRUE)
final <- merge(final, output_onset_event, by = c("species_full", "phenophase"), all.x = TRUE)
final <- merge(final, output_synchrony_sp, by = c("species_full", "phenophase"), all.x = TRUE)

final <- merge(final, census, by = c("species_full"), all.x = TRUE)
final <- merge(final, traits, by = c("species_full"), all.x = TRUE)

# # to get some statistics on full database (total species, individuals, siteyear etc.)
# test <- final %>%
#   filter(phenophase == "leaf_dormancy")
# sum(test$nr_indiv)
# sum(test$site_years)
# sum(test$basal_area_percentage)
# metadata <- read.csv("data/phenology_archives_species_long_format_20200324.csv",
#                      header = TRUE, sep = ",")
# metadata$species_full <- paste(metadata$genus_Meise, metadata$species_Meise)
# metadata = metadata[,(names(metadata) %in% c("species_full",
#                                               "genus_Meise",
#                                               "family_Meise"))]
# metadata <- metadata[!duplicated(metadata), ]
# test <- merge(test, metadata, by = "species_full", all.x = TRUE)
# length(unique(test$genus_Meise)) # unknown included still
# length(unique(test$family_Meise)) # unknown included still


#--------------------------------------------------------------------
#--------------------------------------------------------------------
# write to file
#--------------------------------------------------------------------
#--------------------------------------------------------------------
write.table(final, "data/summ_species_pheno_characteristics.csv",
          quote = FALSE,
          col.names = TRUE,
          row.names = FALSE,
          sep = ",")

write.table(output_synchrony_ind, "data/synchrony_individuals.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")
#--------------------------------------------------------------------




