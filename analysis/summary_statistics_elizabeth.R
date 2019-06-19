#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
# graphics.off()
#----- load required packages --------------------------------------------------#
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(circular)
#----- source required functions -----------------------------------------------#
source("R/event_length.R")
#-------------------------------------------------------------------------------#

#----------------------------------------------------------------------
#--------   Phenology data - species correction Meise   ---------------
#----------------------------------------------------------------------
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
df <- df[which(df$value != 0),]
df$join_id <- paste0("R",df$image,"-",df$image_row)

metadata <- read.csv("data/phenology_archives_species_long_format_20190619.csv",
                     header = TRUE, sep = ",")
metadata$join_id <- paste(metadata$image,metadata$row, sep = "-")

# test merge the two tables based upon the unique join ID
df <- merge(df, metadata, by = c("join_id"), all.x = TRUE)
df$species_full <- paste(df$genus_Meise, df$species_Meise)
# remove column id.x and rename id.y to id (--> in id.y, empty ids are renamed to EK1, EK2, etc...)
df = df[,!(names(df) %in% "id.x")]
df <- df %>%
  rename("id" = id.y)
df$id <- as.character(df$id)
#----------------------------------------------------------------------

# # read in census data and convert basal area
# census <- read.csv2("data/YGB_ForestPlotsNET.csv",
#                      header = TRUE,
#                      sep = ",",
#                      stringsAsFactors = FALSE)
#
# # calculate mean basal area across all species across all
# # mixed plots
# census <- census %>%
#   filter(grepl("MIX",Plot)) %>%
#   group_by(Plot, Species) %>%
#   summarize(basal_area = sum(pi*(C2DBH4/2000)^2, na.rm = TRUE)) %>%
#   group_by(Species) %>%
#   summarize(basal_area = mean(basal_area))
#
# # replace empty slots with Unknown for consistency
# df$genus <- ifelse(df$genus == "","Unknown",df$genus)
# df$species <- ifelse(df$species == "","Unknown",df$species)

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
  rename("species_full" = Species)
census$species_full <- ifelse(census$species_full == "Unknown", NA, census$species_full)

# remove individuals without C1DBH4, these are new recruits for census2 + understory if understory.remove = TRUE
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

# read in trait data from Steven Janssens
traits <- read.csv("data/Dataset_traits_African_trees.csv",
                    header = TRUE,
                    sep = ",",
                    stringsAsFactors = FALSE)
traits <- traits %>%
  rename("mating_system" = Mating_system,
         "ecology" = Ecology,
         "height_m_literature" = Height_lit_m,
         "dbh_max_cm_literature" = Dmax_lit_cm)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# summary statistics: start year, end year, site years, site years with events, nr of individuals
# this is now per id
# counting unique year - id combinations for either
# the full or reduced (phenophase based) subsets
phenophase_selected <- c("leaf_dormancy")
out_LD <- lapply(unique(df$species_full), function(species_selected){

  # list of all observed data regardless of phenophase
  # get then number of unique years
  ss_full <- subset(df, species_full == species_selected,
                    select = c('year','id'))
  nr_years_full <- length(unique(paste(ss_full$year, ss_full$id)))
  nr_indiv <- length(unique(ss_full$id))

  start_year <- min(ss_full$year)
  end_year <- max(ss_full$year)

  # subset based on phenophase of interest + species
  # get the number of unique years (with observations)
  ss_phen <- subset(df,
                    species_full == species_selected &
                      phenophase == phenophase_selected,
                    select = c('year','id'))
  nr_years_phen <- length(unique(paste(ss_phen$year, ss_phen$id)))

  # return the ratio
  ratio <- nr_years_phen/nr_years_full

  return(data.frame(
    "species_full" = species_selected,
    "nr_indiv" = nr_indiv,
    "start_year" = start_year,
    "end_year"= end_year,
    "site_years" = nr_years_full,
    "site_years_with_leaf_dormancy" = nr_years_phen,
    "ratio_site_years_with_leaf_dormancy" = ratio))
})

# bind everything row wise
out_LD <- do.call("rbind", out_LD)
#-------------------------------------------------------------------------------
phenophase_selected <- c("leaf_turnover")
out_LT <- lapply(unique(df$species_full), function(species_selected){

  # list of all observed data regardless of phenophase
  # get then number of unique years
  ss_full <- subset(df, species_full == species_selected,
                    select = c('year','id'))
  nr_years_full <- length(unique(paste(ss_full$year, ss_full$id)))
  nr_indiv <- length(unique(ss_full$id))

  start_year <- min(ss_full$year)
  end_year <- max(ss_full$year)

  # subset based on phenophase of interest + species
  # get the number of unique years (with observations)
  ss_phen <- subset(df,
                    species_full == species_selected &
                      phenophase == phenophase_selected,
                    select = c('year','id'))
  nr_years_phen <- length(unique(paste(ss_phen$year, ss_phen$id)))

  # return the ratio
  ratio <- nr_years_phen/nr_years_full

  return(data.frame(
    "species_full" = species_selected,
    "nr_indiv" = nr_indiv,
    "start_year" = start_year,
    "end_year"= end_year,
    "site_years" = nr_years_full,
    "site_years_with_leaf_turnover" = nr_years_phen,
    "ratio_site_years_with_leaf_turnover" = ratio))
})

# bind everything row wise
out_LT <- do.call("rbind", out_LT)

# merge and drop merge column
out <- merge(out_LD, out_LT, by = c("species_full","nr_indiv","start_year","end_year","site_years"), all.x = TRUE)
rm(out_LD)
rm(out_LT)
#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
#------------ Average duration of an event -------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
event_LD <- df %>%
  filter(phenophase == "leaf_dormancy") %>%
  group_by(species_full, id) %>%
  do(event_length(.))

event_duration_LD <- event_LD %>%
  group_by(species_full) %>%
  dplyr::summarise(mean_duration_leaf_dormancy_weeks = mean(phenophase_length),
                   sd_duration_leaf_dormancy_weeks = sd(phenophase_length),
                   total_nr_events_leaf_dormancy = length(phenophase_length))
# rescaling 48-week year to 52-week year
event_duration_LD <- event_duration_LD %>%
  mutate(mean_duration_leaf_dormancy_weeks = mean_duration_leaf_dormancy_weeks /48 *52,
         sd_duration_leaf_dormancy_weeks = sd_duration_leaf_dormancy_weeks /48 *52)
#-------------------------------------------------------------------------------
event_LT <- df %>%
  filter(phenophase == "leaf_turnover") %>%
  group_by(species_full, id) %>%
  do(event_length(.))

event_duration_LT <- event_LT %>%
  group_by(species_full) %>%
  dplyr::summarise(mean_duration_leaf_turnover_weeks = mean(phenophase_length),
                   sd_duration_leaf_turnover_weeks = sd(phenophase_length),
                   total_nr_events_leaf_turnover = length(phenophase_length))
# rescaling 48-week year to 52-week year
event_duration_LT <- event_duration_LT %>%
  mutate(mean_duration_leaf_turnover_weeks = mean_duration_leaf_turnover_weeks /48 *52,
         sd_duration_leaf_turnover_weeks = sd_duration_leaf_turnover_weeks /48 *52)
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
#------------ Synchrony index among individuals across years -------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# calculate the synchrony between events as mean and SD of onset
# then reformat those taking the year - date line into consideration
transition_dates_LD <- event_LD %>%
  mutate(degree = week_start * 360/48) #%>%

onset_LD <- transition_dates_LD %>%
  group_by(species_full) %>%
  # filter(length(week_start) > 5) %>%    # only species with more than 10 events
  summarise(
    mean_degree = mean.circular(
      circular(degree, units = "degrees")
    ),
    median_degree = median.circular(
      circular(degree, units = "degrees")
    ),
    sd_rad = sd.circular(
      circular(degree, units = "degrees") # eventhough units in degree, sd calculated in radians
    ),                                    # https://stats.stackexchange.com/questions/185361/how-do-i-interpret-the-standard-deviation-of-a-directional-dataset
    nr_events_onset = length(week_start)
  )
# sd: rad back to degrees
onset_LD$sd_degree <- deg(onset_LD$sd_rad)
# mean: rescaling between 0 and 360 degrees
onset_LD$mean_rescaled <- ifelse(onset_LD$mean_degree < 0, onset_LD$mean_degree +360, onset_LD$mean_degree)
onset_LD$median_rescaled <- ifelse(onset_LD$median_degree < 0, onset_LD$median_degree +360, onset_LD$median_degree)
# rescaling degrees to weeks
onset_LD <- onset_LD %>%
  mutate(mean_intrasp_onset_leaf_dormancy_weeks = mean_rescaled /360 *52,
         median_intrasp_onset_leaf_dormancy_weeks = median_rescaled /360 *52,
         sd_intrasp_onset_leaf_dormancy_weeks = sd_degree/360*52)
#-------------------------------------------------------------------------------
#------------ Synchrony index between species as -------------------
#----- average pairwise distance between median onset date of species
#-------------------------------------------------------------------------------
# you need to merge with census here, because else distance is calculated to all other species in full phenology dataset
census_overview <- census[,(names(census) %in% c("species_full","basal_area_site"))]
distance_events_LD <- merge(onset_LD, census_overview, by = "species_full", all.x = TRUE)
distance_events_LD <- distance_events_LD[!(is.na(distance_events_LD$basal_area_site)),]
distance_events_LD$median_rad <- rad(distance_events_LD$median_degree)
distance_events_LD <- distance_events_LD[!(is.na(distance_events_LD$median_rad)),]
a <- as.matrix(dist(distance_events_LD$median_rad), labels = TRUE)
rownames(a) <- distance_events_LD[['species_full']]
colnames(a) <- rownames(a)
b <- deg(a)
c <- ifelse(b > 180, 360-b, b)
d <- as.data.frame(rowMeans(c))
d$species_full <- rownames(d)
colnames(d)[1] <- "mean_distance_onset_leaf_dormancy_weeks"
d <- d %>%
  mutate(mean_distance_onset_leaf_dormancy_weeks = mean_distance_onset_leaf_dormancy_weeks /360 *52)
#---------- now with a mimimum of 5 events
onset_LD_minfreq <- transition_dates_LD %>%
  group_by(species_full) %>%
  filter(length(week_start) >= 5) %>%    # only species with more than 5 events
  summarise(median_degree = median.circular(circular(degree, units = "degrees")))
distance_events_minfreq_LD <- merge(onset_LD_minfreq, census_overview, by = "species_full", all.x = TRUE)
distance_events_minfreq_LD <- distance_events_minfreq_LD[!(is.na(distance_events_minfreq_LD$basal_area_site)),]
distance_events_minfreq_LD$median_rad <- rad(distance_events_minfreq_LD$median_degree)
distance_events_minfreq_LD <- distance_events_minfreq_LD[!(is.na(distance_events_minfreq_LD$median_rad)),]
a_minfreq <- as.matrix(dist(distance_events_minfreq_LD$median_rad), labels = TRUE)
rownames(a_minfreq) <- distance_events_minfreq_LD[['species_full']]
colnames(a_minfreq) <- rownames(a_minfreq)
b_minfreq <- deg(a_minfreq)
c_minfreq <- ifelse(b_minfreq > 180, 360-b_minfreq, b_minfreq)
d_minfreq <- as.data.frame(rowMeans(c_minfreq))
d_minfreq$species_full <- rownames(d_minfreq)
colnames(d_minfreq)[1] <- "mean_distance_onset_minfreq_leaf_dormancy_weeks"
d_minfreq <- d_minfreq %>%
  mutate(mean_distance_onset_minfreq_leaf_dormancy_weeks = mean_distance_onset_minfreq_leaf_dormancy_weeks /360 *52)
#-----------
onset_LD <- merge(onset_LD, d, by = "species_full", all.x = TRUE)
onset_LD <- merge(onset_LD, d_minfreq, by = "species_full", all.x = TRUE)
# remove columns you don't need
onset_LD = onset_LD[,!(names(onset_LD) %in% c("mean_degree","sd_rad","sd_degree","nr_events_onset","mean_rescaled","median_degree","median_rescaled","median_rad"))]
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
#------------ Synchrony index among individuals across years -------------------
#------------Leaf turnover---------------------------------------------
#-------------------------------------------------------------------------------
# calculate the synchrony between events as mean and SD of onset
# then reformat those taking the year - date line into consideration
transition_dates_LT <- event_LT %>%
  mutate(degree = week_start * 360/48) #%>%

onset_LT <- transition_dates_LT %>%
  group_by(species_full) %>%
  # filter(length(week_start) > 5) %>%    # only species with more than 10 events
  summarise(
    mean_degree = mean.circular(
      circular(degree, units = "degrees")
    ),
    median_degree = median.circular(
      circular(degree, units = "degrees")
    ),
    sd_rad = sd.circular(
      circular(degree, units = "degrees") # eventhough units in degree, sd calculated in radians
    ),                                    # https://stats.stackexchange.com/questions/185361/how-do-i-interpret-the-standard-deviation-of-a-directional-dataset
    nr_events_onset = length(week_start)
  )
# sd: rad back to degrees
onset_LT$sd_degree <- deg(onset_LT$sd_rad)
# mean: rescaling between 0 and 360 degrees
onset_LT$mean_rescaled <- ifelse(onset_LT$mean_degree < 0, onset_LT$mean_degree +360, onset_LT$mean_degree)
onset_LT$median_rescaled <- ifelse(onset_LT$median_degree < 0, onset_LT$median_degree +360, onset_LT$median_degree)
# rescaling degrees to weeks
onset_LT <- onset_LT %>%
  mutate(mean_intrasp_onset_leaf_turnover_weeks = mean_rescaled /360 *52,
         median_intrasp_onset_leaf_turnover_weeks = median_rescaled /360 *52,
         sd_intrasp_onset_leaf_turnover_weeks = sd_degree/360*52)
#-------------------------------------------------------------------------------
#------------ Synchrony index between species as -------------------
#----- average pairwise distance between median onset date of species
#-------------------------------------------------------------------------------
# you need to merge with census here, because else distance is calculated to all other species in full phenology dataset
census_overview <- census[,(names(census) %in% c("species_full","basal_area_site"))]
distance_events_LT <- merge(onset_LT, census_overview, by = "species_full", all.x = TRUE)
distance_events_LT <- distance_events_LT[!(is.na(distance_events_LT$basal_area_site)),]
distance_events_LT$median_rad <- rad(distance_events_LT$median_degree)
distance_events_LT <- distance_events_LT[!(is.na(distance_events_LT$median_rad)),]
a <- as.matrix(dist(distance_events_LT$median_rad), labels = TRUE)
rownames(a) <- distance_events_LT[['species_full']]
colnames(a) <- rownames(a)
b <- deg(a)
c <- ifelse(b > 180, 360-b, b)
d <- as.data.frame(rowMeans(c))
d$species_full <- rownames(d)
colnames(d)[1] <- "mean_distance_onset_leaf_turnover_weeks"
d <- d %>%
  mutate(mean_distance_onset_leaf_turnover_weeks = mean_distance_onset_leaf_turnover_weeks /360 *52)
#---------- now with a mimimum of 5 events
onset_LT_minfreq <- transition_dates_LT %>%
  group_by(species_full) %>%
  filter(length(week_start) >= 5) %>%    # only species with more than 5 events
  summarise(median_degree = median.circular(circular(degree, units = "degrees")))
distance_events_minfreq_LT <- merge(onset_LT_minfreq, census_overview, by = "species_full", all.x = TRUE)
distance_events_minfreq_LT <- distance_events_minfreq_LT[!(is.na(distance_events_minfreq_LT$basal_area_site)),]
distance_events_minfreq_LT$median_rad <- rad(distance_events_minfreq_LT$median_degree)
distance_events_minfreq_LT <- distance_events_minfreq_LT[!(is.na(distance_events_minfreq_LT$median_rad)),]
a_minfreq <- as.matrix(dist(distance_events_minfreq_LT$median_rad), labels = TRUE)
rownames(a_minfreq) <- distance_events_minfreq_LT[['species_full']]
colnames(a_minfreq) <- rownames(a_minfreq)
b_minfreq <- deg(a_minfreq)
c_minfreq <- ifelse(b_minfreq > 180, 360-b_minfreq, b_minfreq)
d_minfreq <- as.data.frame(rowMeans(c_minfreq))
d_minfreq$species_full <- rownames(d_minfreq)
colnames(d_minfreq)[1] <- "mean_distance_onset_minfreq_leaf_turnover_weeks"
d_minfreq <- d_minfreq %>%
  mutate(mean_distance_onset_minfreq_leaf_turnover_weeks = mean_distance_onset_minfreq_leaf_turnover_weeks /360 *52)
#-----------
onset_LT <- merge(onset_LT, d, by = "species_full", all.x = TRUE)
onset_LT <- merge(onset_LT, d_minfreq, by = "species_full", all.x = TRUE)
# remove columns you don't need
onset_LT = onset_LT[,!(names(onset_LT) %in% c("mean_degree","sd_rad","sd_degree","nr_events_onset","mean_rescaled","median_degree","median_rescaled","median_rad"))]
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#------------ Synchrony index within individuals across years ------------------
#------------ And then the average SIind of a species --------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
onset_ind_LD <- transition_dates_LD %>%
  group_by(id, species_full) %>%
  # filter(length(week_start) > 5) %>%    # only species with more than 10 events
  summarise(
    mean_degree = mean.circular(
      circular(degree, units = "degrees")
    ),
    sd_rad = sd.circular(
      circular(degree, units = "degrees") # eventhough units in degree, sd calculated in radians
    ),                                    # https://stats.stackexchange.com/questions/185361/how-do-i-interpret-the-standard-deviation-of-a-directional-dataset
    nr_events_onset = length(week_start)
  )
# sd: rad back to degrees
onset_ind_LD$sd_degree <- deg(onset_ind_LD$sd_rad)
# mean: rescaling between 0 and 360 degrees
onset_ind_LD$mean_rescaled <- ifelse(onset_ind_LD$mean_degree < 0, onset_ind_LD$mean_degree +360, onset_ind_LD$mean_degree)
# rescaling degrees to weeks
onset_ind_LD <- onset_ind_LD %>%
  mutate(onset_mean_weeks = mean_rescaled /360 *52,
         onset_sd_weeks = sd_degree/360*52)

synchrony_ind_LD <- onset_ind_LD %>%
  group_by(species_full) %>%
  summarise(mean_synchrony_individuals_onset_leaf_dormancy_weeks = mean(onset_sd_weeks,na.rm=TRUE),
            sd_synchrony_individuals_onset_leaf_dormancy_weeks = sd(onset_sd_weeks,na.rm=TRUE),
            mean_nr_events_within_individuals_leaf_dormancy = mean(nr_events_onset))
#-------------------------------------------------------------------------------
onset_ind_LT <- transition_dates_LT %>%
  group_by(id, species_full) %>%
  # filter(length(week_start) > 5) %>%    # only species with more than 10 events
  summarise(
    mean_degree = mean.circular(
      circular(degree, units = "degrees")
    ),
    sd_rad = sd.circular(
      circular(degree, units = "degrees") # eventhough units in degree, sd calculated in radians
    ),                                    # https://stats.stackexchange.com/questions/185361/how-do-i-interpret-the-standard-deviation-of-a-directional-dataset
    nr_events_onset = length(week_start)
  )
# sd: rad back to degrees
onset_ind_LT$sd_degree <- deg(onset_ind_LT$sd_rad)
# mean: rescaling between 0 and 360 degrees
onset_ind_LT$mean_rescaled <- ifelse(onset_ind_LT$mean_degree < 0, onset_ind_LT$mean_degree +360, onset_ind_LT$mean_degree)
# rescaling degrees to weeks
onset_ind_LT <- onset_ind_LT %>%
  mutate(onset_mean_weeks = mean_rescaled /360 *52,
         onset_sd_weeks = sd_degree/360*52)

synchrony_ind_LT <- onset_ind_LT %>%
  group_by(species_full) %>%
  summarise(mean_synchrony_individuals_onset_leaf_turnover_weeks = mean(onset_sd_weeks,na.rm=TRUE),
            sd_synchrony_individuals_onset_leaf_turnover_weeks = sd(onset_sd_weeks,na.rm=TRUE),
            mean_nr_events_within_individuals_leaf_turnover = mean(nr_events_onset))
#-------------------------------------------------------------------------------




#--------------------------------------------------------------------
#--------------------------------------------------------------------
# merge everything
#--------------------------------------------------------------------
#--------------------------------------------------------------------

final <- merge(out, census, by = "species_full", all.x = TRUE)
final <- merge(final, traits, by = "species_full", all.x = TRUE)
final <- merge(final, event_duration_LD, by = "species_full", all.x = TRUE)
final <- merge(final, event_duration_LT, by = "species_full", all.x = TRUE)
final <- merge(final, onset_LD, by = "species_full", all.x = TRUE)
final <- merge(final, onset_LT, by = "species_full", all.x = TRUE)
final <- merge(final, synchrony_ind_LD, by = "species_full", all.x = TRUE)
final <- merge(final, synchrony_ind_LT, by = "species_full", all.x = TRUE)

# remove rows (species) not included in the Yangambi mixed forest census
final <- final[!(is.na(final$basal_area_site)),]

#--------------------------------------------------------------------
# Sychrony at the individual level
onset_ind_LD = onset_ind_LD[,(names(onset_ind_LD) %in% c("id","species_full","onset_sd_weeks","nr_events_onset"))]
onset_ind_LD$phenophase <- "leaf_dormancy"
onset_ind_LT = onset_ind_LT[,(names(onset_ind_LT) %in% c("id","species_full","onset_sd_weeks","nr_events_onset"))]
onset_ind_LT$phenophase <- "leaf_turnover"
onset_ind <- rbind(onset_ind_LD, onset_ind_LT)
onset_ind <- merge(onset_ind, census, by = "species_full", all.x = TRUE)
onset_ind <- merge(onset_ind, traits, by = "species_full", all.x = TRUE)
# remove rows (species) not included in the Yangambi mixed forest census
onset_ind <- onset_ind[!(is.na(onset_ind$basal_area_site)),]


#--------------------------------------------------------------------
#--------------------------------------------------------------------
# write to file
#--------------------------------------------------------------------
#--------------------------------------------------------------------
write.table(final, "data/species_meta_data.csv",
          quote = FALSE,
          col.names = TRUE,
          row.names = FALSE,
          sep = ",")

write.table(onset_ind, "data/synchrony_individuals.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")
#--------------------------------------------------------------------



# # sum the basal area of the observed species
# # if >30 good coverage
# # message(sum(final$BA, na.rm = TRUE))
#
# hist1 <- ggplot(final, aes(site_years)) +
#   geom_histogram(bins = 15) +
#   labs(x = "Site Years",
#        y = "Frequency") +
#   theme_minimal() +
#   theme(text = element_text(size=20))
#
# hist2 <- ggplot(final, aes(nr_indiv)) +
#   geom_histogram(bins = 15) +
#   labs(x = "Individuals per species",
#        y = "Frequency") +
#   theme_minimal() +
#   theme(text = element_text(size=20))
#
# # ggsave("site_year_histogram.png", hist1)
# # ggsave("individuals_histogram.png", hist2)


