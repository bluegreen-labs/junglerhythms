#----- reset your R session. -----------------------------------------
rm(list=ls())
# graphics.off()
#----- load required packages ----------------------------------------
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)
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
source("analysis/remote_sensing_plot.R")
source("R/event_length.R")
source("R/timeline_gap_fill.R")
source("R/standlevel_phen.R")
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
data <- readRDS("data/jungle_rhythms_data_cleaned.rds")
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
# census_plot <- census %>%
#   filter(grepl("MIX",Plot)) %>%
#   group_by(Plot, species_full) %>%
#   dplyr::summarise(basal_area_plot = sum(basal_area))
census_site <- census %>%
  filter(grepl("MIX",Plot)) %>%
  group_by(species_full) %>%
  dplyr::summarise(basal_area_site = sum(basal_area))

# get total basal area, for weighted means later on
total_basal_area_site <- sum(census_site$basal_area_site)
# total_basal_area_plot <- census_plot %>%
#   group_by(Plot) %>%
#   dplyr::summarise(total_basal_area_plot = sum(basal_area_plot))
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# use function standlevel_phen
# to get stand-level weighted means of phenological pattern
# based on subsetted species lists
#----------------------------------------------------------------------

overview_dorm <- read.csv("data/SI_table2_dormancy.csv",
                          header = TRUE,
                          sep = ",",
                          stringsAsFactors = FALSE)
overview_dorm$dorm_groups <- ifelse(overview_dorm$dorm_groups == "", NA, overview_dorm$dorm_groups)
# # look for species with climate correlations that do not fit in a group
# overview_dorm$dorm_corrs <- ifelse(overview_dorm$corr.precip != "" |
#                                    overview_dorm$corr.sun != "" |
#                                    overview_dorm$corr.tmax != "", "no_group", NA)
# overview_dorm$dorm_groups <- ifelse(is.na(overview_dorm$dorm_groups), overview_dorm$dorm_corrs, overview_dorm$dorm_groups)

# annual with different climate groups
dorm_ann1 <- overview_dorm %>%
  filter(cyclicity_dormancy == "annual") %>%
  filter(dorm_groups == "group1")
dorm_ann2 <- overview_dorm %>%
  filter(cyclicity_dormancy == "annual") %>%
  filter(dorm_groups == "group2")
# dorm_ann_nogroup <- overview_dorm %>%
#   filter(cyclicity_dormancy == "annual") %>%
#   filter(dorm_groups == "no_group")
# dorm_ann_noclim <- overview_dorm %>%
#   filter(cyclicity_dormancy == "annual") %>%
#   filter(is.na(dorm_groups))
dorm_ann_nogroup <- overview_dorm %>%
  filter(cyclicity_dormancy == "annual") %>%
  filter(!dorm_groups %in% c("group1","group2"))
# annual with different climate groups
dorm_group1 <- overview_dorm %>%
  filter(!cyclicity_dormancy %in% c("annual")) %>%
  filter(dorm_groups == "group1")
dorm_group2 <- overview_dorm %>%
  filter(!cyclicity_dormancy %in% c("annual")) %>%
  filter(dorm_groups == "group2")
# dorm_nogroup <- overview_dorm %>%
#   filter(!cyclicity_dormancy == "annual") %>%
#   filter(dorm_groups == "no_group")
# dorm_rest <- overview_dorm %>%
#   filter(!cyclicity_dormancy %in% c("annual")) %>%
#   filter(is.na(dorm_groups))
dorm_rest <- overview_dorm %>%
  filter(!cyclicity_dormancy %in% c("annual")) %>%
  filter(!dorm_groups %in% c("group1","group2"))


overview_turn <- read.csv("data/SI_table2_turnover.csv",
                          header = TRUE,
                          sep = ",",
                          stringsAsFactors = FALSE)
overview_turn$turn_groups <- ifelse(overview_turn$turn_groups == "", NA, overview_turn$turn_groups)
# # look for species with climate correlations that do not fit in a group
# overview_turn$turn_corrs <- ifelse(overview_turn$corr.precip != "" |
#                                      overview_turn$corr.sun != "" |
#                                      overview_turn$corr.tmax != "", "no_group", NA)
# overview_turn$turn_groups <- ifelse(is.na(overview_turn$turn_groups), overview_turn$turn_corrs, overview_turn$turn_groups)

# annual with different climate groups
turn_ann1 <- overview_turn %>%
  filter(cyclicity_turnover == "annual") %>%
  filter(turn_groups == "group1")
turn_ann2 <- overview_turn %>%
  filter(cyclicity_turnover == "annual") %>%
  filter(turn_groups == "group2")
# turn_ann_nogroup <- overview_turn %>%
#   filter(cyclicity_turnover == "annual") %>%
#   filter(turn_groups == "no_group")
# turn_ann_noclim <- overview_turn %>%
#   filter(cyclicity_turnover == "annual") %>%
#   filter(is.na(turn_groups))
turn_ann_nogroup <- overview_turn %>%
  filter(cyclicity_turnover == "annual") %>%
  filter(!turn_groups %in% c("group1","group2"))
# annual with different climate groups
turn_group1 <- overview_turn %>%
  filter(!cyclicity_turnover %in% c("annual")) %>%
  filter(turn_groups == "group1")
turn_group2 <- overview_turn %>%
  filter(!cyclicity_turnover %in% c("annual")) %>%
  filter(turn_groups == "group2")
# turn_nogroup <- overview_turn %>%
#   filter(!cyclicity_turnover == "annual") %>%
#   filter(turn_groups == "no_group")
# turn_rest <- overview_turn %>%
#   filter(!cyclicity_turnover %in% c("annual")) %>%
#   filter(is.na(turn_groups))
turn_rest <- overview_turn %>%
  filter(!cyclicity_turnover %in% c("annual")) %>%
  filter(!turn_groups %in% c("group1","group2"))

dorm_ann1 <- dorm_ann1$Species
dorm_ann2 <- dorm_ann2$Species
dorm_ann_nogroup <- dorm_ann_nogroup$Species
# dorm_ann_noclim <- dorm_ann_noclim$Species
dorm_group1 <- dorm_group1$Species
dorm_group2 <- dorm_group2$Species
# dorm_nogroup <- dorm_nogroup$Species
dorm_rest <- dorm_rest$Species

turn_ann1 <- turn_ann1$Species
turn_ann2 <- turn_ann2$Species
turn_ann_nogroup <- turn_ann_nogroup$Species
# turn_ann_noclim <- turn_ann_noclim$Species
turn_group1 <- turn_group1$Species
turn_group2 <- turn_group2$Species
# turn_nogroup <- turn_nogroup$Species
turn_rest <- turn_rest$Species


#----------------------------------------------------------------------
#----------------------------------------------------------------------
standlevel_ann_group1 <- standlevel_phen(data = data,
                                         species_list_dorm = dorm_ann1,
                                         species_list_turn = turn_ann1)
standlevel_ann_group1$cycl <- "annual"
standlevel_ann_group1$groups <- "group1"

standlevel_ann_group2 <- standlevel_phen(data = data,
                                         species_list_dorm = dorm_ann2,
                                         species_list_turn = turn_ann2)
standlevel_ann_group2$cycl <- "annual"
standlevel_ann_group2$groups <- "group2"

standlevel_ann_nogroup <- standlevel_phen(data = data,
                                          species_list_dorm = turn_ann_nogroup, # don't plot dormancy, flushing or senescence
                                          species_list_turn = turn_ann_nogroup)
standlevel_ann_nogroup$ss <- NA
standlevel_ann_nogroup$cycl <- "annual"
standlevel_ann_nogroup$groups <- "no_class"

# standlevel_ann_noclim <- standlevel_phen(data = data,
#                                           species_list_dorm = dorm_ann_noclim,
#                                           species_list_turn = turn_ann_noclim)
# standlevel_ann_noclim$cycl <- "annual"
# standlevel_ann_noclim$groups <- "no_clim"

standlevel_group1 <- standlevel_phen(data = data,
                                     species_list_dorm = dorm_group1,
                                     species_list_turn = turn_group1)
standlevel_group1$cycl <- "no"
standlevel_group1$groups <- "group1"

standlevel_group2 <- standlevel_phen(data = data,
                                     species_list_dorm = dorm_group2,
                                     species_list_turn = turn_group2)
standlevel_group2$cycl <- "no"
standlevel_group2$groups <- "group2"

# standlevel_nogroup <- standlevel_phen(data = data,
#                                    species_list_dorm = dorm_nogroup,
#                                    species_list_turn = dorm_nogroup) # for turnover only combretum lokele -> too few years for standlevel
# standlevel_nogroup$ss_turn <- NA
# standlevel_nogroup$cycl <- "no"
# standlevel_nogroup$groups <- "no_class"

standlevel_noclim <- standlevel_phen(data = data,
                                      species_list_dorm = dorm_rest,
                                      species_list_turn = turn_rest)
standlevel_noclim$cycl <- "no"
standlevel_noclim$groups <- "no_class" #no_clim

standlevel_subsets <- rbind(standlevel_ann_group1,
                            standlevel_ann_group2,
                            standlevel_ann_nogroup,
                            # standlevel_ann_noclim,
                            standlevel_group1,
                            standlevel_group2,
                            # standlevel_nogroup,
                            standlevel_noclim)

# #----------------------------------------------------------------------
# #----------------------------------------------------------------------
# #----------------------------------------------------------------------
# # check basics for the different groups
# #----------------------------------------------------------------------
# # percentage ba
# census_site$perc_ba <- census_site$basal_area_site / total_basal_area_site *100
# test <- census_site %>%
#   filter(species_full %in% turn_rest) %>%
#   arrange(desc(perc_ba)) %>%
#   mutate(cumsum = cumsum(perc_ba))
# test
# # area under curve
# AUC(standlevel_rest$week, standlevel_rest$ss_turn, method="step") #trapezoid
# AUC(standlevel_rest$week, standlevel_rest$ss, method="step")

#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# manuscript - figure 5
#----------------------------------------------------------------------
#----------------------------------------------------------------------

p1 <- ggplot(data = standlevel_subsets) +
  annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 1.01, alpha = .1) + # jan - febr
  annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 1.01, alpha = .1) + # jun - jul
  annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 1.01, alpha = .1) + # dec
  # turnover
  geom_smooth(aes(week, ss_turn, colour = "line2"), span = 0.2, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(aes(week, ss_turn),
             col="#dfc27d",
             shape = 1,
             stroke = 1.3) +
  # dormancy
  geom_smooth(aes(week, ss, colour = "line1"), span = 0.2, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(aes(week, ss),
             col="#8c510a",
             size=2) +
  scale_colour_manual(values = c("line1" = "#8c510a", "line2" = "#dfc27d"),
                      labels = c(" dormacy   "," turnover   ")) +
  scale_x_continuous(limits = c(1,49),
                     breaks = seq(1,48,4),
                     labels = month.abb) +
  scale_y_continuous(limits = c(0,1.01),
                     breaks = c(0,0.5,1)) +
  labs(y = "% of canopy in state of phenophase",
       x = "") +
  # theme_minimal() +
  theme(text = element_text(family = "Lato", color = "#22211d"),
        panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.major.y = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(size = 11),
        strip.background = element_rect(fill="grey80"), #, linetype="solid", color="grey60",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(vjust = 3),
        legend.position = c(0.1,0.95),
        legend.title = element_blank(),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.background = element_rect(fill = "transparent", colour = NA)
        # legend.text = element_text(size = 11)
        # plot.margin = unit(c(0,0,0,0.5),"cm")
  ) +
  facet_grid(groups ~ cycl)
p1

ann_text1 <- data.frame(week = 22,
                        ss_text = 0.9,
                        groups = c("group1","group1","group2","group2","no_class","no_class"),
                        cycl = c("annual","no","annual","no","annual","no"),
                        notes = c("D: 6 sp, 1.7% BA, 10.5% AUC\nT: 7 sp, 4.2% BA, 19.2% AUC",
                                  "D: 5 sp, 2.0% BA, 6.8% AUC\nT: 11 sp, 7.6% BA, 22.7% AUC",
                                  "D: 2 sp, 2.0% BA, 21.2% AUC\nT: 2 sp, 0.1% BA, 0.2% AUC",
                                  "D: 10 sp, 4.8% BA, 20.9% AUC\nT: 11 sp, 2.9% BA, 8.7% AUC",
                                  "D: 0 sp\nT: 5 sp, 1.7% BA, 7.9% AUC",
                                  "D: 106 sp, 80.7% BA, 40.5% AUC\nT: 93 sp, 74.0% BA, 41.4% AUC"))
ann_segm <- data.frame(x = c(19.5, 23),
                       y = c(0.35, 0.62),
                       xend = c(23, 25.5),
                       yend = c(0.62, 0.62),
                       groups = c("group1","group1"),
                       cycl = c("no","no"))
ann_segm2 <- data.frame(x = c(15, 13),
                       y = c(0.92, 0.92),
                       xend = c(18.5, 15),
                       yend = c(0.7, 0.92),
                       groups = c("no_class","no_class"),
                       cycl = c("no","no"))
ann_text2 <- data.frame(week = c(26.5, 4),
                        ss_text = c(0.55, 0.85),
                        groups = c("group1","no_class"),
                        cycl = c("no","no"),
                        notes = c("P. macrocarpus\n6.6% BA\nperiodically sub-annual",
                                  "S. zenkeri\n16.3% BA\nsingular event"))

p2 <- p1 +
  geom_text(data = ann_text1,
                     aes(week, ss_text, label = notes),
                     size = 3.2,
                     hjust = 0) +
  geom_segment(data = ann_segm,
               aes(x=x,y=y,yend=yend,xend=xend),
               inherit.aes=FALSE) +
  geom_segment(data = ann_segm2,
               aes(x=x,y=y,yend=yend,xend=xend),
               inherit.aes=FALSE) +
  geom_text(data = ann_text2,
            aes(week, ss_text, label = notes),
            size = 3.2,
            hjust = 0,
            # fontface = "italic",
            lineheight = .9)

p2

# use small hack to make 2 modis plots using facet
# by getting facet with same number of columns as p2
# we can force the same plot width via a gtable
VI_s2 <- VI_s
VI_s$hack <- "hack1"
VI_s2$hack <- "hack2"
VI_hack <- rbind(VI_s,VI_s2)

p_modis1 <- ggplot(VI_hack) +
  annotate("rect", xmin = 0, xmax = 61, ymin = 0.25, ymax = 0.77, alpha = .1) + # jan - febr
  annotate("rect", xmin = 152.5, xmax = 213.5, ymin = 0.25, ymax = 0.77, alpha = .1) + # jun - jul
  annotate("rect", xmin = 335.5, xmax = 365, ymin = 0.25, ymax = 0.77, alpha = .1) + # dec
  # geom_point(aes(doy, EVI, shape = site), col = "grey40") +
  geom_ribbon(data = VI_conf, aes(x = doy, ymin = EVImin ,ymax = EVImax), fill="grey60", alpha=0.5) +
  geom_point(data = VI_conf, aes(doy, EVImean), col = "grey30") +
  geom_smooth(aes(doy, EVI), span = 0.3, se = FALSE, col = "grey30", size = 1.2) +
  #geom_line(aes(doy, EVI + EVI_sd)) +
  #geom_line(aes(doy, EVI - EVI_sd)) +
  labs(#title = "MOD13Q1",
    #subtitle = "mean ...",
    x = "",
    y = "EVI") +
  scale_x_continuous(limits = c(0,365),
                     breaks = seq(0,365,30.5),
                     labels = month.abb) +
  # scale_y_continuous(limits = c(0.4,0.6),
  #                    breaks = c(0.4,0.5,0.6),
  #                    labels = scales::number_format(accuracy = 0.1)) +
  theme_minimal() +
  theme(text = element_text(family = "Lato", color = "#22211d"),
        panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_blank(),
        axis.line.x = element_blank(),
        # axis.text.x=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 1.5,size = 10), # vjust to center the label
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 3),
        legend.position = "none"
        # plot.margin=unit(c(0,0,0,0.5),"cm")
  ) +
  facet_grid(~ hack)



p_modis1 <- ggplot_gtable(ggplot_build(p_modis1))
p2 <- ggplot_gtable(ggplot_build(p2))
p_modis1$widths <-p2$widths


p_fig5 <- grid.arrange(p2, p_modis1, heights = c(3,1))

