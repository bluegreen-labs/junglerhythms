#----- reset your R session. -----------------------------------------
rm(list=ls())
# graphics.off()
#----- load required packages ----------------------------------------
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(stats)
library(Hmisc)
library(scales)
library(DescTools)
library(showtext)
font_add_google(
  "Lato",
  regular.wt = 300,
  bold.wt = 700)
#----- source required files  -----------------------------------------
source("analysis/08_remote_sensing_plot.R")
source("R/event_length.R")
source("R/standlevel_phen.R")
source("R/standlevel_phen_plotlevel.R")
#----------------------------------------------------------------------

# Suppress summarise info
# because `summarise()` ungrouping output (override with `.groups` argument)
# comes up a lot in the consol, which is annoying
options(dplyr.summarise.inform = FALSE)
#----- settings  ------------------------------------------------------
# minimum_siteyears: minimum observations years for a species to be included
#                   (cumulative years across individuals of a species)
# -> e.g. if only 1 observation year for a species, an event will have
#     a large impact (weight) on the outcome
# -> this is used in function 'standlevel_phen()'
minimum_siteyears = 10
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#--------   Phenology data --------------------------------------------
#----------------------------------------------------------------------
data <- readRDS("data/jungle_rhythms_data_manuscript_leaf_repro.rds")
species_list <- unique(data$species_full)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#-------- climate data ------------------------------------------------
#----------------------------------------------------------------------
climate <- read.csv("data/ClimData_monthly_avg.csv",
                    header = TRUE,
                    sep = ",")
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#-------- census data  ------------------------------------------------
#-------- get species-specific basal area at plot and site level ------
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
# get sum basal area across all species at plotlevel and at site level
census_plot <- census %>%
  filter(grepl("MIX",Plot)) %>%
  group_by(Plot, species_full) %>%
  dplyr::summarise(basal_area_plot = sum(basal_area))
census_site <- census %>%
  filter(grepl("MIX",Plot)) %>%
  group_by(species_full) %>%
  dplyr::summarise(basal_area_site = sum(basal_area))

# get total basal area, for weighted means later on
total_basal_area_site <- sum(census_site$basal_area_site)
total_basal_area_plot <- census_plot %>%
  group_by(Plot) %>%
  dplyr::summarise(total_basal_area_plot = sum(basal_area_plot))
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# use function standlevel_phen
# to get stand-level weighted means of phenological pattern
# based on total species list
#----------------------------------------------------------------------
standlevel_full <- standlevel_phen(data = data,
                                   census_site = census_site,
                                   total_basal_area_site = total_basal_area_site,
                                   species_list_dorm = species_list,
                                   species_list_turn = species_list,
                                   minimum_siteyears = minimum_siteyears)

#----------------------------------------------------------------------
# use function standlevel_phen_plotlevel
# weekly min/max for the different plots gives the range
# to get stand-level weighted means of phenological pattern
# of plot-level phenological patterns
#----------------------------------------------------------------------
standlevel_full_plots <- standlevel_phen_plotlevel(data = data,
                                                   census_plot = census_plot,
                                                   total_basal_area_plot = total_basal_area_plot,
                                                   species_list_dorm = species_list,
                                                   species_list_turn = species_list,
                                                   minimum_siteyears = minimum_siteyears)

standlevel_range <- standlevel_full_plots %>%
  group_by(week) %>%
  dplyr::summarise(turn_min = min(ss_turn),
                   turn_max = max(ss_turn),
                   dorm_min = min(ss_dorm),
                   dorm_max = max(ss_dorm))
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# manuscript - figure 4
#----------------------------------------------------------------------
#----------------------------------------------------------------------
p_combined_all <- ggplot() +
  annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 2.15, alpha = .1) + # jan - febr
  annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 2.15, alpha = .1) + # jun - jul
  annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 2.15, alpha = .1) + # dec
  annotate("text", x = 37.6, y = 1.55, label = "129 species\n91.2% BA", family = "Lato", color = "#22211d", hjust = 0) +

  # dormancy
  geom_ribbon(data = standlevel_range,
              aes(x = week, ymin = dorm_min, ymax = dorm_max), fill="#8c510a", alpha=0.3) +
  geom_smooth(data = standlevel_full,
              aes(week, ss_dorm, colour = "line1"), span = 0.15, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(data = standlevel_full,
             aes(week, ss_dorm),
             col="#8c510a",
             size=2) +
  # turnover
  geom_ribbon(data = standlevel_range,
              aes(x = week, ymin = turn_min, ymax = turn_max), fill="#018571", alpha=0.3) +
  geom_smooth(data = standlevel_full,
              aes(week, ss_turn, colour = "line2"), span = 0.15, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(data = standlevel_full,
             aes(week, ss_turn),
             col="#018571",
             shape = 1,
             stroke = 1.3) +
  scale_colour_manual(values = c("line1" = "#8c510a", "line2" = "#018571"),
                      labels = c(" senescence   "," turnover   ")) +
  scale_x_continuous(limits = c(1,49),
                     breaks = seq(1,48,4),
                     labels = month.abb) +
  scale_y_continuous(limits = c(0,2.15)) +
  labs(y = "% of canopy in phenophase",
       x = "") +
  theme_minimal() +
  theme(text = element_text(family = "Lato", color = "#22211d"),
        panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 3),
        legend.position = c(0.82,0.9),
        legend.title = element_blank(),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.text = element_text(size = 11),
        plot.margin = unit(c(0,0,0,0.5),"cm")
  )
# p_combined_all

p_precip <- ggplot(climate) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = 0, ymax = 330, alpha = .1) + # jan - febr
  annotate("rect", xmin = 5.5, xmax = 7.5, ymin = 0, ymax = 330, alpha = .1) + # jun - jul
  annotate("rect", xmin = 11.5, xmax = 12.5, ymin = 0, ymax = 330, alpha = .1) + # dec

  annotate("text", x = 1.5, y = 310, label = "LD", family = "Lato", col = "grey50") +
  annotate("text", x = 4, y = 310, label = "SW", family = "Lato", col = "grey50") +
  annotate("text", x = 6.5, y = 310, label = "SD", family = "Lato", col = "grey50") +
  annotate("text", x = 9.5, y = 310, label = "LW", family = "Lato", col = "grey50") +

  geom_col(aes(x = Month,
               y = prec_JR),
           col = "grey70",
           fill = "grey70") +
  geom_linerange(aes(x = Month, ymin = prec_JR_10, ymax = prec_JR_90),
                 col = "grey40") +
  scale_x_continuous(limits = c(0.5,12.5),
                     breaks = seq(0.5,11.5,1),#seq(1,12,1),
                     labels = month.abb) +
  scale_y_continuous(limits = c(0,330),
                     breaks = seq(0,300,50)) +
  labs(y = "precip. (mm)",
       x = "") +
  theme_minimal() +
  theme(text = element_text(family = "Lato", color = "#22211d"),
        panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        # axis.text.x=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 2,size = 10), # vjust to center the label
        axis.title.y = element_text(vjust = 3),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0.5),"cm")
  )

p_sun <- ggplot(climate) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = 120, ymax = 230, alpha = .1) + # jan - febr
  annotate("rect", xmin = 5.5, xmax = 7.5, ymin = 120, ymax = 230, alpha = .1) + # jun - jul
  annotate("rect", xmin = 11.5, xmax = 12.5, ymin = 120, ymax = 230, alpha = .1) + # dec
  geom_ribbon(aes(x = Month, ymin = insol_JR_10 ,ymax = insol_JR_90), fill="grey60", alpha=0.5) +
  geom_line(aes(x = Month,
                y = insol_JR),
            size = 1.2,
            col = "grey30") +
  geom_point(aes(x = Month,
                 y = insol_JR),
             col = "grey30") +
  scale_x_continuous(limits = c(0.5,12.5),
                     breaks = seq(0.5,11.5,1)) +
  scale_y_continuous(limits = c(120,230),
                     breaks = seq(120,230,40)) +
  labs(y = "sun (h)",
       x = "") +
  theme_minimal() +
  theme(text = element_text(family = "Lato", color = "#22211d"),
        panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
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
  geom_ribbon(aes(x = Month, ymin = tmax_JR_10 ,ymax = tmax_JR_90), fill="grey60", alpha=0.5) +
  geom_line(aes(x = Month,
                y = tmax_JR),
            size = 1.2,
            col = "grey30") +
  geom_point(aes(x = Month,
                 y = tmax_JR),
             col = "grey30") +
  scale_x_continuous(limits = c(0.5,12.5),
                     breaks = seq(0.5,11.5,1)) +
  scale_y_continuous(limits = c(28,31.2),
                     breaks = seq(28,31.2,1)) +
  labs(y = expression('t'[max]*' (°C)'),#"tmax (°C)",
       x = "") +
  theme_minimal() +
  theme(text = element_text(family = "Lato", color = "#22211d"),
        panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(vjust = 3),
        legend.position = "none",
        plot.margin=unit(c(0,0,0,0.5),"cm")
  )


#-----------------------------------------------------------------------
# combine plots
#-----------------------------------------------------------------------
p_modis <- ggplot_gtable(ggplot_build(p_modis))
p_precip <- ggplot_gtable(ggplot_build(p_precip))
p_sun <- ggplot_gtable(ggplot_build(p_sun))
p_tmax <- ggplot_gtable(ggplot_build(p_tmax))
p_combined_all <- ggplot_gtable(ggplot_build(p_combined_all))

p_modis$widths <-p_combined_all$widths
p_precip$widths <-p_combined_all$widths
p_sun$widths <-p_combined_all$widths
p_tmax$widths <-p_combined_all$widths

p_all <- grid.arrange(p_combined_all, p_modis, p_tmax, p_sun, p_precip, heights = c(4,4,1.2,1.2,4))#c(4,4,1,1,4))


# #-----------------------------------------------------------------------
# #-----------------------------------------------------------------------
# #-----------------------------------------------------------------------
# #-----------------------------------------------------------------------
# #-----------------------------------------------------------------------
# # correlations between standlevel events and climate
# #-----------------------------------------------------------------------
# #-----------------------------------------------------------------------
# #-----------------------------------------------------------------------
# #-----------------------------------------------------------------------
# #-----------------------------------------------------------------------
# #-----------------------------------------------------------------------
#
# combined <- standlevel_full
# # combined <- standlevel_ann
# # combined <- standlevel_corr
#
# combined$Month <- ifelse(combined$week %in% c(1,2,3,4),1,
#                          ifelse(combined$week %in% c(5,6,7,8),2,
#                                 ifelse(combined$week %in% c(9,10,11,12),3,
#                                        ifelse(combined$week %in% c(13,14,15,16),4,
#                                               ifelse(combined$week %in% c(17,18,19,20),5,
#                                                      ifelse(combined$week %in% c(21,22,23,24),6,
#                                                             ifelse(combined$week %in% c(25,26,27,28),7,
#                                                                    ifelse(combined$week %in% c(29,30,31,32),8,
#                                                                           ifelse(combined$week %in% c(33,34,35,36),9,
#                                                                                  ifelse(combined$week %in% c(37,38,39,40),10,
#                                                                                         ifelse(combined$week %in% c(41,42,43,44),11,
#                                                                                                ifelse(combined$week %in% c(45,46,47,48),12,
#                                                                                                       NA))))))))))))
#
# combined_month <- combined %>%
#   group_by(Month) %>%
#   summarise(dormancy_stand = mean(ss),
#             turnover_stand = mean(ss_turn),
#             senescence_stand = mean(ss_senescence),
#             flushing_stand = mean(ss_flush))
#
# climate.corr <- merge(climate, combined_month, by = c("Month"), all.x = TRUE)
#
#
#
# # turnover
# cor.test(climate.corr$insol_JR, climate.corr$turnover_stand, method = 'pearson')
# cor.test(climate.corr$prec_JR, climate.corr$turnover_stand, method = 'pearson')
# cor.test(climate.corr$tmax_JR, climate.corr$turnover_stand, method = 'pearson')
#
#
# # dormancy
# cor.test(climate.corr$insol_JR, climate.corr$dormancy_stand, method = 'pearson')
# cor.test(climate.corr$prec_JR, climate.corr$dormancy_stand, method = 'pearson')
# cor.test(climate.corr$tmax_JR, climate.corr$dormancy_stand, method = 'pearson')
#
# # flushin
# cor.test(climate.corr$insol_JR, climate.corr$flushing_stand, method = 'pearson')
# cor.test(climate.corr$prec_JR, climate.corr$flushing_stand, method = 'pearson')
# cor.test(climate.corr$tmax_JR, climate.corr$flushing_stand, method = 'pearson')
#
# # senescence
# cor.test(climate.corr$insol_JR, climate.corr$senescence_stand, method = 'pearson')
# cor.test(climate.corr$prec_JR, climate.corr$senescence_stand, method = 'pearson')
# cor.test(climate.corr$tmax_JR, climate.corr$senescence_stand, method = 'pearson')
#
# # climate inter - correlations
# cor.test(climate.corr$insol_JR, climate.corr$prec_JR, method = 'pearson')
# cor.test(climate.corr$insol_all, climate.corr$prec_all, method = 'pearson')
#
# # MODIS
# VI_s$date <- as.Date(VI_s$doy, origin="2010-12-31")
# VI_s$Month <- format(as.Date(VI_s$date), "%m")
# modis <- VI_s %>%
#   group_by(Month)%>%
#   dplyr::summarise(EVIm = mean(EVI))
# modis$Month <- as.numeric(modis$Month)
#
# climate.corr <- merge(climate.corr, modis, by = c("Month"), all.x = TRUE)
#
# cor.test(climate.corr$EVIm, climate.corr$prec_all, method = 'pearson')
# cor.test(climate.corr$EVIm, climate.corr$insol_all, method = 'pearson')
# cor.test(climate.corr$EVIm, climate.corr$PAR_Ygb_Hauser, method = 'pearson')
# cor.test(climate.corr$EVIm, climate.corr$tmax_JR, method = 'pearson')
#
# cor.test(climate.corr$EVIm, climate.corr$dormancy_stand, method = 'pearson')
# cor.test(climate.corr$EVIm, climate.corr$turnover_stand, method = 'pearson')
# cor.test(climate.corr$EVIm, climate.corr$senescence_stand, method = 'pearson')
# cor.test(climate.corr$EVIm, climate.corr$flushing_stand, method = 'pearson')
#
