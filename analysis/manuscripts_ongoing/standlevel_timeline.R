#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
# graphics.off()
#----- load required packages --------------------------------------------------#
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(stats)
library(Hmisc)
library(scales)
#----- source required files  --------------------------------------------------#
# source("analysis/remote_sensing_plot.R")
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
# calculate sum basal area across all species at plotlevel and at site level
census_plot <- census %>%
  filter(grepl("MIX",Plot)) %>%
  group_by(Plot, species_full) %>%
  dplyr::summarise(basal_area_plot = sum(basal_area))
census_site <- census %>%
  filter(grepl("MIX",Plot)) %>%
  group_by(species_full) %>%
  dplyr::summarise(basal_area_site = sum(basal_area))

# calculate total basal area, for weighted means later on
total_basal_area_site <- sum(census_site$basal_area_site)
total_basal_area_plot <- census_plot %>%
  group_by(Plot) %>%
  dplyr::summarise(total_basal_area_plot = sum(basal_area_plot))
#----------------------------------------------------------------------


# #----------------------------------------------------------------------
# #-------- load trait data ---------------------------------------------
# #----------------------------------------------------------------------
# traits <- read.csv("data/Dataset_traits_African_trees.csv",
#                         header = TRUE,
#                         sep = ",")
# traits$layer <- ifelse(traits$Dmax_lit_cm >= understory_dbh_max, "canopy", "understory")

#----------------------------------------------------------------------
#-------- merge census at site_leval data to phenology data  ----------
#-------- just to have the right species, without duplication at plot-level
#----------------------------------------------------------------------
data <- inner_join(data, census_site, by = "species_full")
#----------------------------------------------------------------------
#-- for this species list at ID level: --------------------------------
#-- get full range timelines with 2-year gap filled -------------------
#----------------------------------------------------------------------
test <- data %>%
  filter(!grepl("sp\\.",species_full))
species_list <- unique(test$species_full)
# timelines_dorm <- two_year_gaps(data = data,
#                                 species_name = species_list,
#                                 pheno = "leaf_dormancy")
timelines_turn <- two_year_gaps(data = data,
                                species_name = species_list,
                                pheno = "leaf_turnover")

#-----------------------------------------------------------------------
#------------ Leaf turnover  -------------------------------------------
#-----------------------------------------------------------------------
# stats_LT <- overview_stats(data = timelines_turn,
#                            species_name = species_list,
#                            pheno = "leaf_turnover")
#
# summary(stats_LT$end_year)

# average by date
data_LT <- timelines_turn %>%
  group_by(species_full, date) %>%
  dplyr::summarise(percent_value = mean(value, na.rm=TRUE) *100)

# data_LD <- data_LD %>%
#   filter(total_week >= minimum_siteyears)

data_LT_site <- inner_join(data_LT, census_site, by = c("species_full"))


final_LT_site <- data_LT_site %>%
  group_by(date) %>%
  summarise(ss = sum(percent_value*basal_area_site, na.rm = TRUE)/total_basal_area_site)


p_lin <- ggplot(final_LT_site, aes(x = date,
                             y = ss)) +
  geom_line() +
  theme_minimal() +
  labs(y = "% of individuals with events",
       x = "Year",
       color = "Event") +
  scale_x_date(date_breaks = "1 years",
               date_labels = "%Y",
               limits = as.Date(c('1947-01-01','1952-12-31')),
               expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,10),
                     breaks = c(0,5,10)) +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = -0.25),
        strip.text = element_text(hjust = 0, size = 13, face = "italic"),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.1,0,1,1),"cm")
  ) #+
  # facet_wrap( ~ species_full, ncol = 3)
p_lin




timelines_turn <- timelines_turn %>%
  filter(!is.na(value))

data_LT <- timelines_turn %>%
  group_by(species_full, week) %>%
  summarise(mean_week = mean(value),
            count_week = sum(value),
            total_week = length(value))

# filter out species with too few total siteyears across individuals to have a meaningfull average year
data_LT <- data_LT %>%
  filter(total_week >= minimum_siteyears)

# data_LT <- data_LT %>%
#   filter(mean_week >= minimum_event_frequency)
# data_LT$mean_week <- ifelse(data_LT$mean_week < minimum_event_frequency, 0, data_LT$mean_week)

# if(rescale_event){
#   data_LT <- data_LT %>%
#     group_by(species_full) %>%
#     mutate(mean_week = rescale(count_week, c(minimum_event_frequency,1)))
# }

data_LT_site <- inner_join(data_LT, census_site, by = c("species_full"))
data_LT_plot <- inner_join(data_LT, census_plot, by = c("species_full"))

# # check to see if total basal area per week is constant
# test <- data_LT_plot %>%
#   group_by(Plot, week) %>%
#   summarise(total_basal_area_week = sum(basal_area_plot)) #,

# ### reference website for confidence intervals
# ### https://stats.stackexchange.com/questions/224029/calculating-weighted-cis-and-interpretation
final_LT_site <- data_LT_site %>%
  group_by(week) %>%
  summarise(ss = sum(mean_week*basal_area_site, na.rm = TRUE)/total_basal_area_site *100) #,
# ss2 = wtd.mean(mean_week, basal_area_site, na.rm = TRUE)*100,
# weighted_sd = sqrt(wtd.var(mean_week, basal_area_site, na.rm = TRUE))*100,
# weighted_CI = ss2 + 1.96*weighted_sd)

final_LT_plot <- data_LT_plot %>%
  group_by(Plot, week) %>%
  summarise(ss = sum(mean_week*basal_area_plot, na.rm = TRUE))
final_LT_plot <- inner_join(final_LT_plot, total_basal_area_plot, by = c("Plot"))
final_LT_plot <- final_LT_plot %>%
  mutate(ss = ss/total_basal_area_plot*100)
#-----------------------------------------------------------------------






