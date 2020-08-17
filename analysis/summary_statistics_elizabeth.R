#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
graphics.off()
#----- load required packages --------------------------------------------------#
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(circular)
#----- source required functions -----------------------------------------------#
source("R/event_length.R")
source("R/timeline_gap_fill.R")
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


#-------- read in census data  ----------------------------------------
#-------- get species-specific basal area at plot and site level ------
#----------------------------------------------------------------------
# read in census data and convert basal area
census <- read.csv("data/YGB_ForestPlotsNET_corrected_indet.csv",
                   header = TRUE,
                   sep = ",",
                   stringsAsFactors = FALSE)
census <- census %>%
  dplyr::rename("species_full" = Species)
census$species_full <- ifelse(census$species_full == "Unknown", NA, census$species_full)

# remove individuals without C1DBH4, these are new recruits for census2
# and calculate basal_area for each individual
census <- census[!(is.na(census$C1DBH4)),]
census$basal_area = pi*(census$C1DBH4/2000)^2
# only keep mixed plots
# summarize at species level
census <- census %>%
  filter(grepl("MIX", Plot)) %>%
  group_by(species_full) %>%
  dplyr::summarise(basal_area_site = sum(basal_area)/5,
                   abundance = length(basal_area)/5,
                   dbh_max_mm_census = max(C1DBH4))
basal_area_total <- sum(census$basal_area_site)
census$basal_area_percentage <- census$basal_area_site/basal_area_total *100
#----------------------------------------------------------------------
rm(basal_area_total)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#-- get a species list for analysis to reduce computing time ----------
#-- only species in the census + at species level (not genus) ---------
#----------------------------------------------------------------------
# remove those only at genus level
census_sp <- census %>%
  filter(!grepl("sp\\.",species_full))
# species lists different datasets
census_species <- unique(census_sp$species_full)
phen_species <- unique(data$species_full)
species_list <- intersect(census_species, phen_species)
#----------------------------------------------------------------------
rm(census_species, phen_species, census_sp)
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#-- for this species list at ID level: --------------------------------
#-- get full range timelines with 2-year gap filled -------------------
#----------------------------------------------------------------------
# species_list <- unique(data$species_full)

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


#-------------------------------------------------------------------------------
# summary statistics: start year, end year, site years, site years with events, nr of individuals
#-------------------------------------------------------------------------------
stats_LD <- overview_stats(data = data_timeline,
                           species_name = species_list,
                           pheno = "leaf_dormancy")
stats_LD <- stats_LD %>%
  dplyr::rename("site_years_with_leaf_dormancy" = site_years_with_phenophase,
                "ratio_site_years_with_leaf_dormancy" = ratio_site_years_with_phenophase)
#-------------------------------------------------------------------------------
stats_LT <- overview_stats(data = data_timeline,
                           species_name = species_list,
                           pheno = "leaf_turnover")
stats_LT <- stats_LT %>%
  dplyr::rename("site_years_with_leaf_turnover" = site_years_with_phenophase,
                "ratio_site_years_with_leaf_turnover" = ratio_site_years_with_phenophase)
#-------------------------------------------------------------------------------
# merge and drop merge column
stats_LD <- stats_LD %>%
  select(-phenophase)
stats_LT <- stats_LT %>%
  select(-phenophase)
output_stats <- merge(stats_LD, stats_LT, by = c("species_full","nr_indiv","start_year","end_year","site_years"), all.x = TRUE)
rm(stats_LD, stats_LT)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#------------ Average duration of an event -------------------------------------
#-------------------------------------------------------------------------------
event_LD <- event_length(data = data_timeline,
                         species_name = species_list,
                         pheno = "leaf_dormancy")
# mean and sd of length phenophase
# rescaling 48-week year to 52-week year
event_duration_LD <- event_LD %>%
  group_by(species_full) %>%
  dplyr::summarise(mean_duration_leaf_dormancy_weeks = mean(phenophase_length, na.rm = TRUE),
                   sd_duration_leaf_dormancy_weeks = sd(phenophase_length, na.rm = TRUE),
                   total_nr_events_leaf_dormancy = length(phenophase_length[!(is.na(phenophase_length))]))
event_duration_LD <- event_duration_LD %>%
  mutate(mean_duration_leaf_dormancy_weeks = mean_duration_leaf_dormancy_weeks /48 *52,
         sd_duration_leaf_dormancy_weeks = sd_duration_leaf_dormancy_weeks /48 *52)
#-------------------------------------------------------------------------------
event_LT <- event_length(data = data_timeline,
                         species_name = species_list,
                         pheno = "leaf_turnover")
event_duration_LT <- event_LT %>%
  group_by(species_full) %>%
  dplyr::summarise(mean_duration_leaf_turnover_weeks = mean(phenophase_length, na.rm = TRUE),
                   sd_duration_leaf_turnover_weeks = sd(phenophase_length, na.rm = TRUE),
                   total_nr_events_leaf_turnover = length(phenophase_length[!(is.na(phenophase_length))]))
event_duration_LT <- event_duration_LT %>%
  mutate(mean_duration_leaf_turnover_weeks = mean_duration_leaf_turnover_weeks /48 *52,
         sd_duration_leaf_turnover_weeks = sd_duration_leaf_turnover_weeks /48 *52)
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
#------------ Synchrony index among individuals across years -------------------
#-------------------------------------------------------------------------------
# calculate the synchrony between events as mean and SD of onset
# then reformat those taking the year - date line into consideration
transition_dates_LD <- event_LD %>%
  filter(!is.na(week_start)) %>%
  mutate(degree = week_start * 360/48)

onset_LD <- transition_dates_LD %>%
  group_by(species_full) %>%
  dplyr::summarise(
    mean_degree = mean.circular(circular(degree, units = "degrees")),
    median_degree = median.circular(circular(degree, units = "degrees")),
    sd_rad = sd.circular(circular(degree, units = "degrees")), # eventhough units in degree, sd calculated in radians # https://stats.stackexchange.com/questions/185361/how-do-i-interpret-the-standard-deviation-of-a-directional-dataset
    nr_events_onset = length(week_start)) %>%
  mutate(sd_degree = deg(sd_rad))
# only sd for species with more than 5 events
onset_LD$sd_degree <- ifelse(onset_LD$nr_events_onset >= 3, onset_LD$sd_degree, NA)
# mean: rescaling between 0 and 360 degrees
onset_LD$mean_rescaled <- ifelse(onset_LD$mean_degree < 0, onset_LD$mean_degree +360, onset_LD$mean_degree)
onset_LD$median_rescaled <- ifelse(onset_LD$median_degree < 0, onset_LD$median_degree +360, onset_LD$median_degree)
# rescaling degrees to weeks
# and only keep columns you need
onset_LD <- onset_LD %>%
  mutate(mean_intrasp_onset_leaf_dormancy_weeks = mean_rescaled /360 *52,
         median_intrasp_onset_leaf_dormancy_weeks = median_rescaled /360 *52,
         sd_intrasp_onset_leaf_dormancy_weeks = sd_degree/360*52) %>%
  select(species_full,
         mean_intrasp_onset_leaf_dormancy_weeks,
         median_intrasp_onset_leaf_dormancy_weeks,
         sd_intrasp_onset_leaf_dormancy_weeks)
#-------------------------------------------------------------------------------
transition_dates_LT <- event_LT %>%
  filter(!is.na(week_start)) %>%
  mutate(degree = week_start * 360/48)

onset_LT <- transition_dates_LT %>%
  group_by(species_full) %>%
  dplyr::summarise(
    mean_degree = mean.circular(circular(degree, units = "degrees")),
    median_degree = median.circular(circular(degree, units = "degrees")),
    sd_rad = sd.circular(circular(degree, units = "degrees")),
    nr_events_onset = length(week_start)) %>%
  mutate(sd_degree = deg(sd_rad))
# only sd for species with more than 5 events
onset_LT$sd_degree <- ifelse(onset_LT$nr_events_onset >= 3, onset_LT$sd_degree, NA)
# mean: rescaling between 0 and 360 degrees
onset_LT$mean_rescaled <- ifelse(onset_LT$mean_degree < 0, onset_LT$mean_degree +360, onset_LT$mean_degree)
onset_LT$median_rescaled <- ifelse(onset_LT$median_degree < 0, onset_LT$median_degree +360, onset_LT$median_degree)
# rescaling degrees to weeks
# and only keep columns you need
onset_LT <- onset_LT %>%
  mutate(mean_intrasp_onset_leaf_turnover_weeks = mean_rescaled /360 *52,
         median_intrasp_onset_leaf_turnover_weeks = median_rescaled /360 *52,
         sd_intrasp_onset_leaf_turnover_weeks = sd_degree/360*52) %>%
  select(species_full,
         mean_intrasp_onset_leaf_turnover_weeks,
         median_intrasp_onset_leaf_turnover_weeks,
         sd_intrasp_onset_leaf_turnover_weeks)
#-------------------------------------------------------------------------------


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
#------------ Synchrony index within individuals across years ------------------
#------------ And then the average SIind of a species --------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
onset_ind_LD <- transition_dates_LD %>%
  group_by(id, species_full) %>%
  filter(length(week_start) >= 3) %>%    # only species with more than 4 events
  dplyr::summarise(mean_degree = mean.circular(circular(degree, units = "degrees")),
    sd_rad = sd.circular(circular(degree, units = "degrees")), # eventhough units in degree, sd calculated in radians
    nr_events_onset = length(week_start)) %>%
  mutate(sd_degree = deg(sd_rad))
# mean: rescaling between 0 and 360 degrees
onset_ind_LD$mean_rescaled <- ifelse(onset_ind_LD$mean_degree < 0, onset_ind_LD$mean_degree +360, onset_ind_LD$mean_degree)
# rescaling degrees to weeks
onset_ind_LD <- onset_ind_LD %>%
  mutate(onset_mean_weeks = mean_rescaled /360 *52,
         onset_sd_weeks = sd_degree/360*52)

synchrony_ind_LD <- onset_ind_LD %>%
  group_by(species_full) %>%
  dplyr::summarise(mean_synchrony_individuals_onset_leaf_dormancy_weeks = mean(onset_sd_weeks,na.rm=TRUE),
                   sd_synchrony_individuals_onset_leaf_dormancy_weeks = sd(onset_sd_weeks,na.rm=TRUE),
                   mean_nr_events_within_individuals_leaf_dormancy = mean(nr_events_onset))
#-------------------------------------------------------------------------------
onset_ind_LT <- transition_dates_LT %>%
  group_by(id, species_full) %>%
  filter(length(week_start) >= 3) %>%    # only species with more than 4 events
  dplyr::summarise(mean_degree = mean.circular(circular(degree, units = "degrees")),
                   sd_rad = sd.circular(circular(degree, units = "degrees")), # eventhough units in degree, sd calculated in radians
                   nr_events_onset = length(week_start)) %>%
  mutate(sd_degree = deg(sd_rad))
# mean: rescaling between 0 and 360 degrees
onset_ind_LT$mean_rescaled <- ifelse(onset_ind_LT$mean_degree < 0, onset_ind_LT$mean_degree +360, onset_ind_LT$mean_degree)
# rescaling degrees to weeks
onset_ind_LT <- onset_ind_LT %>%
  mutate(onset_mean_weeks = mean_rescaled /360 *52,
         onset_sd_weeks = sd_degree/360*52)

synchrony_ind_LT <- onset_ind_LT %>%
  group_by(species_full) %>%
  dplyr::summarise(mean_synchrony_individuals_onset_leaf_turnover_weeks = mean(onset_sd_weeks,na.rm=TRUE),
                   sd_synchrony_individuals_onset_leaf_turnover_weeks = sd(onset_sd_weeks,na.rm=TRUE),
                   mean_nr_events_within_individuals_leaf_turnover = mean(nr_events_onset))
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Sychrony at the individual level
# will be writen out as seperate file
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
onset_ind_LD <- onset_ind_LD %>%
  select(species_full,
         id,
         onset_sd_weeks,
         nr_events_onset) %>%
  mutate(phenophase = "leaf_dormancy")

onset_ind_LT <- onset_ind_LT %>%
  select(species_full,
         id,
         onset_sd_weeks,
         nr_events_onset) %>%
  mutate(phenophase = "leaf_turnover")

onset_ind <- rbind(onset_ind_LD, onset_ind_LT)
#-------------------------------------------------------------------------------
rm(onset_ind_LD, onset_ind_LT)
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
final <- merge(output_stats, census, by = "species_full", all.x = TRUE)
final <- merge(final, traits, by = "species_full", all.x = TRUE)
final <- merge(final, event_duration_LD, by = "species_full", all.x = TRUE)
final <- merge(final, event_duration_LT, by = "species_full", all.x = TRUE)
final <- merge(final, onset_LD, by = "species_full", all.x = TRUE)
final <- merge(final, onset_LT, by = "species_full", all.x = TRUE)
final <- merge(final, synchrony_ind_LD, by = "species_full", all.x = TRUE)
final <- merge(final, synchrony_ind_LT, by = "species_full", all.x = TRUE)


# # to get some statistics on full database (total species, individuals, siteyear etc.)
# sum(final$nr_indiv)
# sum(final$site_years)
# sum(final$basal_area_percentage)
# metadata <- read.csv("data/phenology_archives_species_long_format_20200324.csv",
#                      header = TRUE, sep = ",")
# metadata$species_full <- paste(metadata$genus_Meise, metadata$species_Meise)
# metadata = metadata[,(names(metadata) %in% c("species_full",
#                                               "genus_Meise",
#                                               "family_Meise"))]
# metadata <- metadata[!duplicated(metadata), ]
# final <- merge(final, metadata, by = "species_full", all.x = TRUE)
# length(unique(final$genus_Meise)) # unknown included still
# length(unique(final$family_Meise)) # unknown included still
#
#
# onset_ind <- merge(onset_ind, census, by = "species_full", all.x = TRUE)
# onset_ind <- merge(onset_ind, traits, by = "species_full", all.x = TRUE)
#--------------------------------------------------------------------
#--------------------------------------------------------------------
# write to file
#--------------------------------------------------------------------
#--------------------------------------------------------------------
write.table(final, "data/species_meta_data_phase2.csv",
          quote = FALSE,
          col.names = TRUE,
          row.names = FALSE,
          sep = ",")

write.table(onset_ind, "data/synchrony_individuals_phase2.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")
#--------------------------------------------------------------------




