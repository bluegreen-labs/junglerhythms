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
zero_events_remove <- FALSE  # in this case FALSE because zero_events are used to calculate event frequency
# and zero_events are needed so that sum(basal_area) per week is constant
understory_remove <- FALSE
understory_dbh_max = 0 # unit = cm
minimum_siteyears = 10
minimum_event_frequency = 0 # range c(0,1)
# rescale_event <- FALSE
#-------------------------------------------------------------------------------#

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# when understory_remove = TRUE, total_basal_area is not correct
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#----------------------------------------------------------------------
#--------   Phenology data - species correction Meise   ---------------
#----------------------------------------------------------------------
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
if(zero_events_remove){
  df <- df[which(df$value != 0),]
}

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

# # remove trees in understory
# if(understory_remove){
#   census$C1DBH4 <- ifelse(census$C1DBH4 >= (understory_dbh_max*10), census$C1DBH4, NA) #*10 because units here is mm
# }

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
#-------- load climate data -------------------------------------------
#----------------------------------------------------------------------
climate <- read.csv("data/ClimData_monthly_avg.csv",
                    header = TRUE,
                    sep = ",")
#----------------------------------------------------------------------




#-----------------------------------------------------------------------
#------------ flowering  -------------------------------------------
#-----------------------------------------------------------------------
data_flower <- data %>%
  filter(phenophase == "flowers") %>%
  # filter(layer == "canopy") %>%
  # filter(Ecology == "shade") %>%
  na.omit()

data_flower <- data_flower %>%
  group_by(species_full, week) %>%
  summarise(mean_week = mean(value, na.rm = TRUE),
            count_week = sum(value),
            total_week = length(value))

# filter out species with too few total siteyears to have a meaningfull average year
data_flower <- data_flower %>%
  filter(total_week >= minimum_siteyears)

# data_flower <- data_flower %>%
#   filter(mean_week >= minimum_event_frequency)
data_flower$mean_week <- ifelse(data_flower$mean_week < minimum_event_frequency, 0, data_flower$mean_week)

# if(rescale_event){
#   data_flower <- data_flower %>%
#     group_by(species_full) %>%
#     mutate(mean_week = rescale(count_week, c(minimum_event_frequency,1)))
# }

data_flower_site <- inner_join(data_flower, census_site, by = c("species_full"))
data_flower_plot <- inner_join(data_flower, census_plot, by = c("species_full"))

# # check to see if total basal area per week is constant
# test <- data_flower_plot %>%
#   group_by(Plot, week) %>%
#   summarise(total_basal_area_week = sum(basal_area_plot)) #,

# ### reference website for confidence intervals
# ### https://stats.stackexchange.com/questions/224029/calculating-weighted-cis-and-interpretation
final_flower_site <- data_flower_site %>%
  group_by(week) %>%
  summarise(ss = sum(mean_week*basal_area_site, na.rm = TRUE)/total_basal_area_site *100) #,
# ss2 = wtd.mean(mean_week, basal_area_site, na.rm = TRUE)*100,
# weighted_sd = sqrt(wtd.var(mean_week, basal_area_site, na.rm = TRUE))*100,
# weighted_CI = ss2 + 1.96*weighted_sd)

final_flower_plot <- data_flower_plot %>%
  group_by(Plot, week) %>%
  summarise(ss = sum(mean_week*basal_area_plot, na.rm = TRUE))
final_flower_plot <- inner_join(final_flower_plot, total_basal_area_plot, by = c("Plot"))
final_flower_plot <- final_flower_plot %>%
  mutate(ss = ss/total_basal_area_plot*100)
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#------------ fruiting  -------------------------------------------
#-----------------------------------------------------------------------
data_fruit <- data %>%
  filter(phenophase == "fruit") %>%
  # filter(layer == "canopy") %>%
  # filter(Ecology == "shade") %>%
  na.omit()

data_fruit <- data_fruit %>%
  group_by(species_full, week) %>%
  summarise(mean_week = mean(value, na.rm = TRUE),
            count_week = sum(value),
            total_week = length(value))

# filter out species with too few total siteyears to have a meaningfull average year
data_fruit <- data_fruit %>%
  filter(total_week >= minimum_siteyears)

# data_fruit <- data_fruit %>%
#   filter(mean_week >= minimum_event_frequency)
data_fruit$mean_week <- ifelse(data_fruit$mean_week < minimum_event_frequency, 0, data_fruit$mean_week)

# if(rescale_event){
#   data_fruit <- data_fruit %>%
#     group_by(species_full) %>%
#     mutate(mean_week = rescale(count_week, c(minimum_event_frequency,1)))
# }

data_fruit_site <- inner_join(data_fruit, census_site, by = c("species_full"))
data_fruit_plot <- inner_join(data_fruit, census_plot, by = c("species_full"))

final_fruit_site <- data_fruit_site %>%
  group_by(week) %>%
  summarise(ss = sum(mean_week*basal_area_site, na.rm = TRUE)/total_basal_area_site *100) #,
# ss2 = wtd.mean(mean_week, basal_area_site, na.rm = TRUE)*100,
# weighted_sd = sqrt(wtd.var(mean_week, basal_area_site, na.rm = TRUE))*100,
# weighted_CI = ss2 + 1.96*weighted_sd)

final_fruit_plot <- data_fruit_plot %>%
  group_by(Plot, week) %>%
  summarise(ss = sum(mean_week*basal_area_plot, na.rm = TRUE))
final_fruit_plot <- inner_join(final_fruit_plot, total_basal_area_plot, by = c("Plot"))
final_fruit_plot <- final_fruit_plot %>%
  mutate(ss = ss/total_basal_area_plot*100)
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#------------ dispersal  -------------------------------------------
#-----------------------------------------------------------------------
data_dispersal <- data %>%
  filter(phenophase == "fruit_drop") %>%
  # filter(layer == "canopy") %>%
  # filter(Ecology == "shade") %>%
  na.omit()

data_dispersal <- data_dispersal %>%
  group_by(species_full, week) %>%
  summarise(mean_week = mean(value, na.rm = TRUE),
            count_week = sum(value),
            total_week = length(value))

# filter out species with too few total siteyears to have a meaningfull average year
data_dispersal <- data_dispersal %>%
  filter(total_week >= minimum_siteyears)

# data_dispersal <- data_dispersal %>%
#   filter(mean_week >= minimum_event_frequency)
data_dispersal$mean_week <- ifelse(data_dispersal$mean_week < minimum_event_frequency, 0, data_dispersal$mean_week)

# if(rescale_event){
#   data_dispersal <- data_dispersal %>%
#     group_by(species_full) %>%
#     mutate(mean_week = rescale(count_week, c(minimum_event_frequency,1)))
# }

data_dispersal_site <- inner_join(data_dispersal, census_site, by = c("species_full"))
data_dispersal_plot <- inner_join(data_dispersal, census_plot, by = c("species_full"))

# # check to see if total basal area per week is constant
# test <- data_dispersal_plot %>%
#   group_by(Plot, week) %>%
#   summarise(total_basal_area_week = sum(basal_area_plot)) #,

# ### reference website for confidence intervals
# ### https://stats.stackexchange.com/questions/224029/calculating-weighted-cis-and-interpretation
final_dispersal_site <- data_dispersal_site %>%
  group_by(week) %>%
  summarise(ss = sum(mean_week*basal_area_site, na.rm = TRUE)/total_basal_area_site *100) #,
# ss2 = wtd.mean(mean_week, basal_area_site, na.rm = TRUE)*100,
# weighted_sd = sqrt(wtd.var(mean_week, basal_area_site, na.rm = TRUE))*100,
# weighted_CI = ss2 + 1.96*weighted_sd)

final_dispersal_plot <- data_dispersal_plot %>%
  group_by(Plot, week) %>%
  summarise(ss = sum(mean_week*basal_area_plot, na.rm = TRUE))
final_dispersal_plot <- inner_join(final_dispersal_plot, total_basal_area_plot, by = c("Plot"))
final_dispersal_plot <- final_dispersal_plot %>%
  mutate(ss = ss/total_basal_area_plot*100)
#-----------------------------------------------------------------------




#-----------------------------------------------------------------------
#------------ plots -------  -------------------------------------------
#-----------------------------------------------------------------------

p_flowering <- ggplot() +
  geom_point(data = final_flower_plot,
             aes(week, ss, shape = Plot),
             col="grey40") +
  geom_smooth(data = final_flower_site,
              aes(week, ss), span = 0.2, se = FALSE, col = "black", size = 1.2) +
  # geom_line(data = final_flower_site,
  #           aes(week, ss),
  #           col="black",
  #           size=1.2) +
  geom_point(data = final_flower_site,
             aes(week, ss),
             col="black",
             size=2) +
  scale_x_continuous(limits = c(1,49),
                     breaks = seq(1,48,4),
                     labels = month.abb) +
  # scale_y_continuous(limits = c(0,2.5),
  #                    breaks = seq(0,2,0.5),
  #                    labels = scales::number_format(accuracy = 0.1)) +
  annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 2.2, alpha = .1) + # jan - febr
  annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 2.2, alpha = .1) + # jun - jul
  annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 2.2, alpha = .1) + # dec
  labs(y = "freq. canopy flowering",
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
        # axis.text.x = element_text(angle = 90, hjust = 1,vjust = 1.5,size = 10), # vjust to center the label
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 3),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0.5),"cm")
  )
# print(p_flowering)


p_fruit <- ggplot() +
  geom_point(data = final_fruit_plot,
             aes(week, ss, shape = Plot),
             col="grey40") +
  geom_smooth(data = final_fruit_site,
              aes(week, ss), span = 0.2, se = FALSE, col = "black", size = 1.2) +
  # geom_line(data = final_fruit_site,
  #           aes(week, ss),
  #           col="black",
  #           size=1.2) +
  geom_point(data = final_fruit_site,
             aes(week, ss),
             col="black",
             size=2) +
  scale_x_continuous(limits = c(1,49),
                     breaks = seq(1,48,4),
                     labels = month.abb) +
  # scale_y_continuous(limits = c(0,2.5),
  #                    breaks = seq(0,2,0.5),
  #                    labels = scales::number_format(accuracy = 0.1)) +
  annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 2.2, alpha = .1) + # jan - febr
  annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 2.2, alpha = .1) + # jun - jul
  annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 2.2, alpha = .1) + # dec
  labs(y = "freq. canopy fruiting",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 2,size = 10), # vjust to center the label
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 3),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0.5),"cm")
  )
# print(p_fruit)


p_dispersal <- ggplot() +
  geom_point(data = final_dispersal_plot,
             aes(week, ss, shape = Plot),
             col="grey40") +
  geom_smooth(data = final_dispersal_site,
              aes(week, ss), span = 0.2, se = FALSE, col = "black", size = 1.2) +
  # geom_line(data = final_dispersal_site,
  #           aes(week, ss),
  #           col="black",
  #           size=1.2) +
  geom_point(data = final_dispersal_site,
             aes(week, ss),
             col="black",
             size=2) +
  scale_x_continuous(limits = c(1,49),
                     breaks = seq(1,48,4),
                     labels = month.abb) +
  # scale_y_continuous(limits = c(0,2.5),
  #                    breaks = seq(0,2,0.5),
  #                    labels = scales::number_format(accuracy = 0.1)) +
  annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 2.2, alpha = .1) + # jan - febr
  annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 2.2, alpha = .1) + # jun - jul
  annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 2.2, alpha = .1) + # dec
  labs(y = "freq. canopy seed dispersal",
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
        # axis.text.x = element_text(angle = 90, hjust = 1,vjust = 1.5,size = 10), # vjust to center the label
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 3),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0.5),"cm")
  )
# print(p_dispersal)


p_precip <- ggplot(climate) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = 0, ymax = 330, alpha = .1) + # jan - febr
  annotate("rect", xmin = 5.5, xmax = 7.5, ymin = 0, ymax = 330, alpha = .1) + # jun - jul
  annotate("rect", xmin = 11.5, xmax = 12.5, ymin = 0, ymax = 330, alpha = .1) + # dec
  geom_col(aes(x = Month,
               y = prec_JR),
           col = "grey70",
           fill = "grey70") +
  geom_linerange(aes(x = Month, ymin = prec_JR_10, ymax = prec_JR_90),
                 col = "grey40") +
  scale_x_continuous(limits = c(0.5,12.5),
                     breaks = seq(1,12,1),
                     labels = month.abb) +
  scale_y_continuous(limits = c(0,330),
                     breaks = seq(0,300,50)) +
  labs(y = "precip. (mm)",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), #element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        # axis.text.x=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size = 10), # vjust to center the label
        axis.title.y = element_text(vjust = 3),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0.5),"cm")
  )

p_sun <- ggplot(climate) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = 120, ymax = 230, alpha = .1) + # jan - febr
  annotate("rect", xmin = 5.5, xmax = 7.5, ymin = 120, ymax = 230, alpha = .1) + # jun - jul
  annotate("rect", xmin = 11.5, xmax = 12.5, ymin = 120, ymax = 230, alpha = .1) + # dec
  geom_ribbon(aes(x = Month, ymin = insol_JR_10 ,ymax = insol_JR_90), fill="grey60", alpha="0.5") +
  geom_line(aes(x = Month,
                y = insol_JR),
            size = 1.2,
            col = "grey30") +
  geom_point(aes(x = Month,
                 y = insol_JR),
             col = "grey30") +
  scale_x_continuous(limits = c(0.5,12.5),
                     breaks = seq(1,12,1),
                     labels = month.abb) +
  scale_y_continuous(limits = c(120,230),
                     breaks = seq(120,230,40)) +
  labs(y = "sun (h)",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), # element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        # axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(vjust = 3),
        legend.position = "none",
        plot.margin=unit(c(0,0,0,0.5),"cm")
  )
# p_sun

p_tmax <- ggplot(climate) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = 28, ymax = 31.2, alpha = .1) + # jan - febr
  annotate("rect", xmin = 5.5, xmax = 7.5, ymin = 28, ymax = 31.2, alpha = .1) + # jun - jul
  annotate("rect", xmin = 11.5, xmax = 12.5, ymin = 28, ymax = 31.2, alpha = .1) + # dec
  geom_ribbon(aes(x = Month, ymin = tmax_JR_10 ,ymax = tmax_JR_90), fill="grey60", alpha="0.5") +
  geom_line(aes(x = Month,
                y = tmax_JR),
            size = 1.2,
            col = "grey30") +
  geom_point(aes(x = Month,
                 y = tmax_JR),
             col = "grey30") +
  scale_x_continuous(limits = c(0.5,12.5),
                     breaks = seq(1,12,1),
                     labels = month.abb) +
  scale_y_continuous(limits = c(28,31.2),
                     breaks = seq(28,31.2,1)) +
  labs(y = "tmax (°C)",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), # element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        # axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(vjust = 3),
        legend.position = "none",
        plot.margin=unit(c(0,0,0,0.5),"cm")
  )


#-----------------------------------------------------------------------
# combine plots
#-----------------------------------------------------------------------
# make the widths of the plots with different y-axis equal
# https://stackoverflow.com/questions/30402930/align-x-axes-of-box-plot-and-line-plot-using-ggplot
p_flowering <- ggplot_gtable(ggplot_build(p_flowering))
p_precip <- ggplot_gtable(ggplot_build(p_precip))
p_sun <- ggplot_gtable(ggplot_build(p_sun))
p_tmax <- ggplot_gtable(ggplot_build(p_tmax))
p_fruit <- ggplot_gtable(ggplot_build(p_fruit))
p_dispersal <- ggplot_gtable(ggplot_build(p_dispersal))

p_precip$widths <-p_flowering$widths
p_sun$widths <-p_flowering$widths
p_tmax$widths <-p_flowering$widths
p_fruit$widths <-p_flowering$widths
p_dispersal$widths <-p_flowering$widths

p_all <- grid.arrange(p_flowering, p_fruit,p_dispersal, p_tmax, p_sun, p_precip, heights = c(4,4,4,1,1,4))


pdf("~/Desktop/repro_plotlevel.pdf",4,8.5) # 5,10)
plot(p_all)
dev.off()


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# correlations between standlevel events and climate
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

final_flower_site$Month <- ifelse(final_flower_site$week %in% c(1,2,3,4),1,
                              ifelse(final_flower_site$week %in% c(5,6,7,8),2,
                                     ifelse(final_flower_site$week %in% c(9,10,11,12),3,
                                            ifelse(final_flower_site$week %in% c(13,14,15,16),4,
                                                   ifelse(final_flower_site$week %in% c(17,18,19,20),5,
                                                          ifelse(final_flower_site$week %in% c(21,22,23,24),6,
                                                                 ifelse(final_flower_site$week %in% c(25,26,27,28),7,
                                                                        ifelse(final_flower_site$week %in% c(29,30,31,32),8,
                                                                               ifelse(final_flower_site$week %in% c(33,34,35,36),9,
                                                                                      ifelse(final_flower_site$week %in% c(37,38,39,40),10,
                                                                                             ifelse(final_flower_site$week %in% c(41,42,43,44),11,
                                                                                                    ifelse(final_flower_site$week %in% c(45,46,47,48),12,
                                                                                                           NA))))))))))))

final_flower_month <- final_flower_site %>%
  group_by(Month) %>%
  summarise(flower_stand = mean(ss))

final_fruit_site$Month <- ifelse(final_fruit_site$week %in% c(1,2,3,4),1,
                              ifelse(final_fruit_site$week %in% c(5,6,7,8),2,
                                     ifelse(final_fruit_site$week %in% c(9,10,11,12),3,
                                            ifelse(final_fruit_site$week %in% c(13,14,15,16),4,
                                                   ifelse(final_fruit_site$week %in% c(17,18,19,20),5,
                                                          ifelse(final_fruit_site$week %in% c(21,22,23,24),6,
                                                                 ifelse(final_fruit_site$week %in% c(25,26,27,28),7,
                                                                        ifelse(final_fruit_site$week %in% c(29,30,31,32),8,
                                                                               ifelse(final_fruit_site$week %in% c(33,34,35,36),9,
                                                                                      ifelse(final_fruit_site$week %in% c(37,38,39,40),10,
                                                                                             ifelse(final_fruit_site$week %in% c(41,42,43,44),11,
                                                                                                    ifelse(final_fruit_site$week %in% c(45,46,47,48),12,
                                                                                                           NA))))))))))))

final_fruit_month <- final_fruit_site %>%
  group_by(Month) %>%
  summarise(fruit_stand = mean(ss))

climate.corr <- merge(climate, final_flower_month, by = c("Month"), all.x = TRUE)
climate.corr <- merge(climate.corr, final_fruit_month, by = c("Month"), all.x = TRUE)

# flowering
# cor.test(climate.corr$PAR_Ygb_Yoko_Hauser, climate.corr$flower_stand, method = 'pearson')
# cor.test(climate.corr$PAR_Ygb_Hauser, climate.corr$flower_stand, method = 'pearson')
# cor.test(climate.corr$insol_all, climate.corr$flower_stand, method = 'pearson')
cor.test(climate.corr$insol_JR, climate.corr$flower_stand, method = 'pearson')
# cor.test(climate.corr$prec_all, climate.corr$flower_stand, method = 'pearson')
cor.test(climate.corr$prec_JR, climate.corr$flower_stand, method = 'pearson')
cor.test(climate.corr$tmax_JR, climate.corr$flower_stand, method = 'pearson')


# fruit
# cor.test(climate.corr$PAR_Ygb_Yoko_Hauser, climate.corr$fruit_stand, method = 'pearson')
# cor.test(climate.corr$PAR_Ygb_Hauser, climate.corr$fruit_stand, method = 'pearson')
# cor.test(climate.corr$insol_all, climate.corr$fruit_stand, method = 'pearson')
cor.test(climate.corr$insol_JR, climate.corr$fruit_stand, method = 'pearson')
# cor.test(climate.corr$prec_all, climate.corr$fruit_stand, method = 'pearson')
cor.test(climate.corr$prec_JR, climate.corr$fruit_stand, method = 'pearson')
cor.test(climate.corr$tmax_JR, climate.corr$fruit_stand, method = 'pearson')

# climate inter - correlations
cor.test(climate.corr$insol_JR, climate.corr$prec_JR, method = 'pearson')
cor.test(climate.corr$insol_all, climate.corr$prec_all, method = 'pearson')

# MODIS
VI_s$date <- as.Date(VI_s$doy, origin="2010-12-31")
VI_s$Month <- format(as.Date(VI_s$date), "%m")
modis <- VI_s %>%
  group_by(Month)%>%
  dplyr::summarise(EVIm = mean(EVI))
modis$Month <- as.numeric(modis$Month)

climate.corr <- merge(climate.corr, modis, by = c("Month"), all.x = TRUE)

cor.test(climate.corr$EVIm, climate.corr$prec_all, method = 'pearson')
cor.test(climate.corr$EVIm, climate.corr$insol_all, method = 'pearson')
cor.test(climate.corr$EVIm, climate.corr$PAR_Ygb_Hauser, method = 'pearson')
cor.test(climate.corr$EVIm, climate.corr$fruit_stand, method = 'pearson')
cor.test(climate.corr$EVIm, climate.corr$flower_stand, method = 'pearson')
