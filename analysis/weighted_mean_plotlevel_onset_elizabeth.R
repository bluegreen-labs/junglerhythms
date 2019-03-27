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
library(circular)
#----- source required files  --------------------------------------------------#
source("analysis/remote_sensing_plot.R")
source("R/event_length.R")
#-------------------------------------------------------------------------------#


#-------------------------------------------------------------------------------#
# SETTINGS                                                                      #
# Change these setting according to the output you want.                        #
#-------------------------------------------------------------------------------#
# zero_events_remove: remove events = 0 from initial annotations                #
# understory.remove:  remove individuals in the census with understory_dbh_max  #
# understory_dbh_max: every dbh above this is assumed to be in the canopy       #
#                     also used for defining classes 'understory','canopy'      #
#                     from traits                                               #
# minimum_event_frequency: lower freq. of events removed from upscaling         #
#                          to standlevel                                        #
#         --> this is especially important when rescaling event frequency       #
#         --> to avoid low low freq. events to drive the signal                 #
# rescale_event: at species-level, event freq. will be rescaled                 #
#                between c(minimum_event_frequency,1)                           #
#-------------------------------------------------------------------------------#
zero_events_remove <- TRUE  # in this case TRUE because this reduces calculation time for function event_length
understory_remove <- TRUE
understory_dbh_max = 30 # unit = cm
minimum_event_frequency = 0 # range c(0,1)
rescale_event <- FALSE
minimum_events = 5 # from onset
#-------------------------------------------------------------------------------#



#----------------------------------------------------------------------
#--------   Phenology data - species correction Meise   ---------------
#----------------------------------------------------------------------
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
if(zero_events_remove){
  df <- df[which(df$value != 0),]
}

df$join_id <- paste0("R",df$image,"-",df$image_row)
metadata <- read.csv("data/phenology_archives_species_long_format_20190319.csv",
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
data$species_full <- as.character(data$species_full)
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

# remove trees in understory
if(understory_remove){
  census$C1DBH4 <- ifelse(census$C1DBH4 >= (understory_dbh_max*10), census$C1DBH4, NA) #*10 because units here is mm
}
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
  group_by(species_full) %>%
  dplyr::summarise(basal_area_site = sum(basal_area))
#----------------------------------------------------------------------


# #----------------------------------------------------------------------
# #-------- load trait data ---------------------------------------------
# #----------------------------------------------------------------------
# traits <- read.csv("data/Dataset_traits_African_trees.csv",
#                    header = TRUE,
#                    sep = ",")
# traits$layer <- ifelse(traits$Dmax_lit_cm >= understory_dbh_max, "canopy", "understory")

#----------------------------------------------------------------------
#-------- merge census & traits data to phenology data  ---------------
#----------------------------------------------------------------------
# data <- inner_join(data, census_plot, by = "species_full")
data <- inner_join(data, census_site, by = "species_full")
# data <- inner_join(data, traits, by = "species_full")

#----------------------------------------------------------------------
#-------- load trait data ---------------------------------------------
#----------------------------------------------------------------------
climate <- read.csv("~/Dropbox/Phenology_JR/Manuscript/Phenology_Leaf_draft/ClimData_test.csv",
                    header = TRUE,
                    sep = ",")
#----------------------------------------------------------------------


#-----------------------------------------------------------------------
#------------ Leaf turnover
#------------ median onset
#-----------------------------------------------------------------------
transition_dates <- data %>%
  filter(phenophase == "leaf_turnover") %>%
  group_by(species_full, id, phenophase) %>%
  do(event_length(.))
# weeks to degrees
transition_dates <- transition_dates %>%
  mutate(degree = (week_start-1) * 360/48) # -1 so that week 1 in janurari is degree 0
# summary of onset dates per species
onset <- transition_dates %>%
  group_by(species_full) %>%
  filter(length(week_start) >= minimum_events) %>%    # only species with a minimum of events
  summarise(median_degree = median.circular(
      circular(degree, units = "degrees")
    )
  )
# merge with census and remove rows (species) not included in the Yangambi mixed forest census
onset <- merge(onset, census_site, by = "species_full", all.x = TRUE)
onset <- onset[!(is.na(onset$basal_area_site)),]
# rescaling between 0 and 360 degrees
onset$median_rescaled <- ifelse(onset$median_degree < 0, onset$median_degree +360, onset$median_degree)
# degrees to weeks
onset <- onset %>%
  mutate(week = (median_rescaled*48/360)+1)

onset$freq = 1

#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#-----------something wrong with calculation of weighted mean------------------------------------------------
#-----------first check in weigthed_mean_plotlevel_elizabeth --------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
total_basal_area_site <- sum(onset$basal_area_site)


final_LT_site <- onset %>%
  group_by(week) %>%
  summarise(ss = sum(freq*basal_area_site, na.rm = TRUE)/total_basal_area_site)

# final_LT_plot <- data_LT_plot %>%
#   group_by(Plot, week) %>%
#   summarise(ss = sum(mean_week*basal_area_plot, na.rm = TRUE)/(sum(basal_area_plot))*100)
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------
transition_dates_LD <- data %>%
  filter(phenophase == "leaf_dormancy") %>%
  group_by(species_full, id, phenophase) %>%
  do(event_length(.))
# weeks to degrees
transition_dates_LD <- transition_dates_LD %>%
  mutate(degree = (week_start-1) * 360/48) # -1 so that week 1 in janurari is degree 0
# summary of onset dates per species
onset_LD <- transition_dates_LD %>%
  group_by(species_full) %>%
  filter(length(week_start) >= minimum_events) %>%    # only species with a minimum of events
  summarise(median_degree = median.circular(
    circular(degree, units = "degrees")
  )
  )
# merge with census and remove rows (species) not included in the Yangambi mixed forest census
onset_LD <- merge(onset_LD, census_site, by = "species_full", all.x = TRUE)
onset_LD <- onset_LD[!(is.na(onset_LD$basal_area_site)),]
# rescaling between 0 and 360 degrees
onset_LD$median_rescaled <- ifelse(onset_LD$median_degree < 0, onset_LD$median_degree +360, onset_LD$median_degree)
# degrees to weeks
onset_LD <- onset_LD %>%
  mutate(week = (median_rescaled*48/360)+1)

onset_LD$freq = 1

#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#-----------something wrong with calculation of weighted mean------------------------------------------------
#-----------first check in weigthed_mean_plotlevel_elizabeth --------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
total_basal_area_site_LD <- sum(onset_LD$basal_area_site)


final_LD_site <- onset_LD %>%
  group_by(week) %>%
  summarise(ss = sum(freq*basal_area_site, na.rm = TRUE)/total_basal_area_site_LD)

# final_LT_plot <- data_LT_plot %>%
#   group_by(Plot, week) %>%
#   summarise(ss = sum(mean_week*basal_area_plot, na.rm = TRUE)/(sum(basal_area_plot))*100)
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------
#
# #-----------------------------------------------------------------------
# #------------ Leaf dormancy  -------------------------------------------
# #-----------------------------------------------------------------------
# data_LD <- data %>%
#   filter(phenophase == "leaf_dormancy") %>%
#   # filter(layer == "canopy") %>%
#   # filter(Ecology == "shade") %>%
#   na.omit()
#
# data_LD <- data_LD %>%
#   group_by(species_full, week) %>%
#   summarise(mean_week = mean(value, na.rm = TRUE),
#             count_week = sum(value),
#             total_week = length(value))
#
# data_LD <- data_LD %>%
#   filter(mean_week >= minimum_event_frequency)
#
# if(rescale_event){
#   data_LD <- data_LD %>%
#     group_by(species_full) %>%
#     mutate(mean_week = rescale(count_week, c(minimum_event_frequency,1)))
# }
#
# data_LD_site <- inner_join(data_LD, census_site, by = c("species_full"))
# data_LD_plot <- inner_join(data_LD, census_plot, by = c("species_full"))
#
# final_LD_site <- data_LD_site %>%
#   group_by(week) %>%
#   summarise(ss = sum(mean_week*basal_area_site, na.rm = TRUE)/(sum(basal_area_site))*100)
#
# final_LD_plot <- data_LD_plot %>%
#   group_by(Plot, week) %>%
#   summarise(ss = sum(mean_week*basal_area_plot, na.rm = TRUE)/(sum(basal_area_plot))*100)
# #-----------------------------------------------------------------------



#-----------------------------------------------------------------------
#------------ plots -------  -------------------------------------------
#-----------------------------------------------------------------------
p_turnover <- ggplot() +
  # geom_point(data = final_LT_plot,
  #            aes(week, ss, shape = Plot),
  #            col="grey40") +
  # geom_smooth(data = final_LT_site,
  #             aes(week, ss), span = 0.2, se = FALSE, col = "red", size = 1.2) +
  geom_line(data = final_LT_site,
            aes(week, ss),
            col="red",
            size=1.2) +
  geom_point(data = final_LT_site,
             aes(week, ss),
             col="red",
             size=2) +
  scale_x_continuous(limits = c(1,49),
                     breaks = seq(1,48,4),
                     labels = month.abb) +
  scale_y_continuous(limits = c(0,0.3),
                     breaks = seq(0,0.3,0.1),
                     labels = scales::number_format(accuracy = 0.1)) +
  annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 2.8, alpha = .2) + # jan - febr
  annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 2.8, alpha = .2) + # jun - jul
  annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 2.8, alpha = .2) + # dec
  labs(y = "freq. leaf turnover",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,0,0.2,0.2),"cm")
  )
print(p_turnover)


p_dormancy <- ggplot() +
  # geom_point(data = final_LD_plot,
  #            aes(week, ss, shape = Plot),
  #            col="grey40") +
  # geom_smooth(data = final_LD_site,
  #             aes(week, ss), span = 0.2, se = FALSE, col = "red", size = 1.2) +
  geom_line(data = final_LD_site,
            aes(week, ss),
            col="red",
            size=1.2) +
  geom_point(data = final_LD_site,
             aes(week, ss),
             col="red",
             size=2) +
  scale_x_continuous(limits = c(1,49),
                     breaks = seq(1,48,4),
                     labels = month.abb) +
  scale_y_continuous(limits = c(0,0.3),
                     breaks = seq(0,0.3,0.1),
                     labels = scales::number_format(accuracy = 0.1)) +
  annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 2.8, alpha = .2) + # jan - febr
  annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 2.8, alpha = .2) + # jun - jul
  annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 2.8, alpha = .2) + # dec
  labs(y = "freq. leaf dormancy",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,0,0.2,0.2),"cm")
  )
print(p_dormancy)


p_precip <- ggplot(climate) +
  geom_col(aes(x = Month,
               y = prec),
           col = "grey70",
           fill = "grey70") +
  scale_x_continuous(limits = c(0.5,12.5),
                     breaks = seq(1,12,1),
                     labels = month.abb) +
  labs(y = "precip. (mm)",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0.2),"cm")
  )

p_par <- ggplot(climate) +
  geom_line(aes(x = Month,
                y = PAR_Hauser),
            size = 1.2) +
  scale_x_continuous(limits = c(0.5,12.5),
                     breaks = seq(1,12,1),
                     labels = month.abb) +
  labs(y = "PAR",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        # axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none",
        plot.margin=unit(c(0,0,0,0.2),"cm")
  )


grid.arrange(p_modis, p_dormancy, p_turnover, p_par,p_precip, heights = c(3,3,3.4,1,2)) #
# grid.arrange(p_dormancy,p_turnover,p_climate, heights = c(1,1,1))
