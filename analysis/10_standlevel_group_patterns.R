#----- reset your R session. -----------------------------------------
rm(list=ls())
graphics.off()
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
# based on subsetted species lists
#----------------------------------------------------------------------
overview_dorm <- read.csv("data/SI_table2_dormancy.csv",
                          header = TRUE,
                          sep = ",",
                          stringsAsFactors = FALSE)
# overview_dorm$groups <- ifelse(overview_dorm$groups == "", NA, overview_dorm$groups)
overview_dorm$groups <- ifelse(overview_dorm$total_nr_events == 0, "no_dormancy", overview_dorm$groups)
# annual with different climate groups
dorm_ann1 <- overview_dorm %>%
  filter(cycle_category %in% c("annual","sub-annual")) %>%
  filter(groups == "group1")
dorm_ann2 <- overview_dorm %>%
  filter(cycle_category %in% c("annual","sub-annual")) %>%
  filter(groups == "group2")
dorm_ann_nogroup <- overview_dorm %>%
  filter(cycle_category %in% c("annual","sub-annual")) %>%
  filter(!groups %in% c("group1","group2","no_dormancy"))
# annual with different climate groups
dorm_group1 <- overview_dorm %>%
  filter(!cycle_category %in% c("annual","sub-annual")) %>%
  filter(groups == "group1")
dorm_group2 <- overview_dorm %>%
  filter(!cycle_category %in% c("annual","sub-annual")) %>%
  filter(groups == "group2")
dorm_rest <- overview_dorm %>%
  filter(!cycle_category %in% c("annual","sub-annual")) %>%
  filter(!groups %in% c("group1","group2","no_dormancy"))


overview_turn <- read.csv("data/SI_table3_turnover.csv",
                          header = TRUE,
                          sep = ",",
                          stringsAsFactors = FALSE)
# overview_turn$turn_groups <- ifelse(overview_turn$turn_groups == "", NA, overview_turn$turn_groups)
overview_turn$groups <- ifelse(overview_turn$total_nr_events == 0, "no_turnover", overview_turn$groups)

# annual with different climate groups
turn_ann1 <- overview_turn %>%
  filter(cycle_category %in% c("annual","sub-annual")) %>%
  filter(groups == "group1")
turn_ann2 <- overview_turn %>%
  filter(cycle_category %in% c("annual","sub-annual")) %>%
  filter(groups == "group2")
turn_ann_nogroup <- overview_turn %>%
  filter(cycle_category %in% c("annual","sub-annual")) %>%
  filter(!groups %in% c("group1","group2","no_turnover"))
# annual with different climate groups
turn_group1 <- overview_turn %>%
  filter(!cycle_category %in% c("annual","sub-annual")) %>%
  filter(groups == "group1")
turn_group2 <- overview_turn %>%
  filter(!cycle_category %in% c("annual","sub-annual")) %>%
  filter(groups == "group2")
turn_rest <- overview_turn %>%
  filter(!cycle_category %in% c("annual","sub-annual")) %>%
  filter(!groups %in% c("group1","group2","no_turnover"))

dorm_ann1 <- dorm_ann1$species_full
dorm_ann2 <- dorm_ann2$species_full
dorm_ann_nogroup <- dorm_ann_nogroup$species_full
dorm_group1 <- dorm_group1$species_full
dorm_group2 <- dorm_group2$species_full
dorm_rest <- dorm_rest$species_full

turn_ann1 <- turn_ann1$species_full
turn_ann2 <- turn_ann2$species_full
turn_ann_nogroup <- turn_ann_nogroup$species_full
turn_group1 <- turn_group1$species_full
turn_group2 <- turn_group2$species_full
turn_rest <- turn_rest$species_full

#----------------------------------------------------------------------
# annual, group 1
standlevel_ann_group1 <- standlevel_phen(data = data,
                                         census_site = census_site,
                                         total_basal_area_site = total_basal_area_site,
                                         species_list_dorm = dorm_ann1,
                                         species_list_turn = turn_ann1,
                                         minimum_siteyears = minimum_siteyears)
standlevel_ann_group1$cycl <- "annual"
standlevel_ann_group1$groups <- "group1"
# annual, group 2
standlevel_ann_group2 <- standlevel_phen(data = data,
                                         census_site = census_site,
                                         total_basal_area_site = total_basal_area_site,
                                         species_list_dorm = dorm_ann2,
                                         species_list_turn = turn_ann2,
                                         minimum_siteyears = minimum_siteyears)
standlevel_ann_group2$cycl <- "annual"
standlevel_ann_group2$groups <- "group2"
# annual, no_class
standlevel_ann_nogroup <- standlevel_phen(data = data,
                                          census_site = census_site,
                                          total_basal_area_site = total_basal_area_site,
                                          species_list_dorm = dorm_ann_nogroup,
                                          species_list_turn = turn_ann_nogroup,
                                          minimum_siteyears = minimum_siteyears)
standlevel_ann_nogroup$cycl <- "annual"
standlevel_ann_nogroup$groups <- "no_class"
# no cyc, group 1
standlevel_group1 <- standlevel_phen(data = data,
                                     census_site = census_site,
                                     total_basal_area_site = total_basal_area_site,
                                     species_list_dorm = dorm_group1,
                                     species_list_turn = turn_group1,
                                     minimum_siteyears = minimum_siteyears)
standlevel_group1$cycl <- "no"
standlevel_group1$groups <- "group1"
# no cyc, group 2
standlevel_group2 <- standlevel_phen(data = data,
                                     census_site = census_site,
                                     total_basal_area_site = total_basal_area_site,
                                     species_list_dorm = dorm_group2,
                                     species_list_turn = turn_group2,
                                     minimum_siteyears = minimum_siteyears)
standlevel_group2$cycl <- "no"
standlevel_group2$groups <- "group2"
# no cyc, no_class
standlevel_noclim <- standlevel_phen(data = data,
                                     census_site = census_site,
                                     total_basal_area_site = total_basal_area_site,
                                      species_list_dorm = dorm_rest,
                                      species_list_turn = turn_rest,
                                     minimum_siteyears = minimum_siteyears)
standlevel_noclim$cycl <- "no"
standlevel_noclim$groups <- "no_class"

standlevel_subsets <- rbind(standlevel_ann_group1,
                            standlevel_ann_group2,
                            standlevel_ann_nogroup,
                            standlevel_group1,
                            standlevel_group2,
                            standlevel_noclim)

#----------------------------------------------------------------------
# use function standlevel_phen_plotlevel
# weekly min/max for the different plots gives the range
# to get stand-level weighted means of phenological pattern
# of plot-level phenological patterns
#----------------------------------------------------------------------
# annual, group 1
standlevel_ann_group1_plot <- standlevel_phen_plotlevel(data = data,
                                                        census_plot = census_plot,
                                                        total_basal_area_plot = total_basal_area_plot,
                                                        species_list_dorm = dorm_ann1,
                                                        species_list_turn = turn_ann1,
                                                        minimum_siteyears = minimum_siteyears)

standlevel_ann_group1_range <- standlevel_ann_group1_plot %>%
  group_by(week) %>%
  dplyr::summarise(turn_min = min(ss_turn),
                   turn_max = max(ss_turn, na.rm=TRUE),
                   dorm_min = min(ss_dorm),
                   dorm_max = max(ss_dorm, na.rm=TRUE))
standlevel_ann_group1_range$cycl <- "annual"
standlevel_ann_group1_range$groups <- "group1"
# annual, group2
standlevel_ann_group2_plot <- standlevel_phen_plotlevel(data = data,
                                                        census_plot = census_plot,
                                                        total_basal_area_plot = total_basal_area_plot,
                                                        species_list_dorm = dorm_ann2,
                                                        species_list_turn = turn_ann2,
                                                        minimum_siteyears = minimum_siteyears)

standlevel_ann_group2_range <- standlevel_ann_group2_plot %>%
  group_by(week) %>%
  dplyr::summarise(turn_min = min(ss_turn),
                   turn_max = max(ss_turn, na.rm=TRUE),
                   dorm_min = min(ss_dorm),
                   dorm_max = max(ss_dorm, na.rm=TRUE))
standlevel_ann_group2_range$cycl <- "annual"
standlevel_ann_group2_range$groups <- "group2"
standlevel_ann_group2_range$turn_min <- ifelse(is.na(standlevel_ann_group2_range$turn_min), 0, standlevel_ann_group2_range$turn_min)
standlevel_ann_group2_range$dorm_min <- ifelse(is.na(standlevel_ann_group2_range$dorm_min), 0, standlevel_ann_group2_range$dorm_min)
# annual, no_class
standlevel_ann_nogroup_plot <- standlevel_phen_plotlevel(data = data,
                                                         census_plot = census_plot,
                                                         total_basal_area_plot = total_basal_area_plot,
                                                         species_list_dorm = dorm_ann_nogroup,
                                                         species_list_turn = turn_ann_nogroup,
                                                         minimum_siteyears = minimum_siteyears)

standlevel_ann_nogroup_range <- standlevel_ann_nogroup_plot %>%
  group_by(week) %>%
  dplyr::summarise(turn_min = min(ss_turn),
                   turn_max = max(ss_turn, na.rm=TRUE),
                   dorm_min = min(ss_dorm),
                   dorm_max = max(ss_dorm))
standlevel_ann_nogroup_range$cycl <- "annual"
standlevel_ann_nogroup_range$groups <- "no_class"
# no cyc, group 1
standlevel_group1_plot <- standlevel_phen_plotlevel(data = data,
                                                    census_plot = census_plot,
                                                    total_basal_area_plot = total_basal_area_plot,
                                                    species_list_dorm = dorm_group1,
                                                    species_list_turn = turn_group1,
                                                    minimum_siteyears = minimum_siteyears)

standlevel_group1_range <- standlevel_group1_plot %>%
  group_by(week) %>%
  dplyr::summarise(turn_min = min(ss_turn),
                   turn_max = max(ss_turn, na.rm=TRUE),
                   dorm_min = min(ss_dorm),
                   dorm_max = max(ss_dorm, na.rm=TRUE))
standlevel_group1_range$cycl <- "no"
standlevel_group1_range$groups <- "group1"
# no cyc, group 2
standlevel_group2_plot <- standlevel_phen_plotlevel(data = data,
                                                    census_plot = census_plot,
                                                    total_basal_area_plot = total_basal_area_plot,
                                                    species_list_dorm = dorm_group2,
                                                    species_list_turn = turn_group2,
                                                    minimum_siteyears = minimum_siteyears)

standlevel_group2_range <- standlevel_group2_plot %>%
  group_by(week) %>%
  dplyr::summarise(turn_min = min(ss_turn),
                   turn_max = max(ss_turn, na.rm=TRUE),
                   dorm_min = min(ss_dorm),
                   dorm_max = max(ss_dorm, na.rm=TRUE))
standlevel_group2_range$cycl <- "no"
standlevel_group2_range$groups <- "group2"
# no cyc, no_class
standlevel_noclim_plot <- standlevel_phen_plotlevel(data = data,
                                                    census_plot = census_plot,
                                                    total_basal_area_plot = total_basal_area_plot,
                                                    species_list_dorm = dorm_rest,
                                                    species_list_turn = turn_rest,
                                                    minimum_siteyears = minimum_siteyears)

standlevel_noclim_range <- standlevel_noclim_plot %>%
  group_by(week) %>%
  dplyr::summarise(turn_min = min(ss_turn),
                   turn_max = max(ss_turn, na.rm=TRUE),
                   dorm_min = min(ss_dorm),
                   dorm_max = max(ss_dorm, na.rm=TRUE))
standlevel_noclim_range$cycl <- "no"
standlevel_noclim_range$groups <- "no_class"

ranges_subsets <- rbind(standlevel_ann_group1_range,
                        standlevel_ann_group2_range,
                        standlevel_ann_nogroup_range,
                        standlevel_group1_range,
                        standlevel_group2_range,
                        standlevel_noclim_range)

#----------------------------------------------------------------------

# #----------------------------------------------------------------------
# #----------------------------------------------------------------------
# #----------------------------------------------------------------------
# # check basics for the different groups
# #----------------------------------------------------------------------
# # percentage ba
# census_site$perc_ba <- census_site$basal_area_site / total_basal_area_site *100
# test <- census_site %>%
#   filter(species_full %in% dorm_rest) %>% # dorm_ann1, dorm_ann2, dorm_group1, dorm_ann_nogroup, dorm_rest
#   arrange(desc(perc_ba)) %>%
#   mutate(cumsum = cumsum(perc_ba))
# test
# # area under curve
# AUC(standlevel_noclim$week, standlevel_noclim$ss_turn, method="step") # standlevel_ann_group1, standlevel_group1, standlevel_ann_nogroup, standlevel_noclim
# AUC(standlevel_noclim$week, standlevel_noclim$ss_dorm, method="step")

#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# manuscript - figure 5
#----------------------------------------------------------------------
#----------------------------------------------------------------------
cycl.names <- c("annual" = "(sub-)annual", "no" = "non-annual")
group.names <- c("group1" = "class 1\nin-phase: - prep; + sun / tmax", "group2" = "class 2\nlag: - sun / tmax","no_class" = "no class")

p1 <- ggplot(data = standlevel_subsets) +
  annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 1.93, alpha = .1) + # jan - febr #1.61
  annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 1.93, alpha = .1) + # jun - jul
  annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 1.93, alpha = .1) + # dec
  # turnover
  geom_ribbon(data = ranges_subsets,
              aes(x = week, ymin = turn_min, ymax = turn_max), fill="#018571", alpha=0.3) +
  geom_smooth(aes(week, ss_turn, colour = "line2"), span = 0.2, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(aes(week, ss_turn),
             col="#018571",
             shape = 1,
             stroke = 1.3) +
  # dormancy
  geom_ribbon(data = ranges_subsets,
              aes(x = week, ymin = dorm_min, ymax = dorm_max), fill="#8c510a", alpha=0.3) +
  geom_smooth(aes(week, ss_dorm, colour = "line1"), span = 0.2, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(aes(week, ss_dorm),
             col="#8c510a",
             size=2) +

  scale_colour_manual(values = c("line1" = "#8c510a", "line2" = "#018571"),
                      labels = c(" senescence   "," turnover   ")) +
  scale_x_continuous(limits = c(1,49),
                     breaks = seq(1,48,4),
                     labels = month.abb) +
  scale_y_continuous(limits = c(0,1.93),
                     breaks = c(0,0.5,1,1.5)) +
  labs(y = "% of canopy in phenophase",
       x = "") +
  # theme_minimal() +
  theme(text = element_text(family = "Lato", color = "#22211d"),
        panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.major.y = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(size = 10),
        strip.background = element_rect(fill="grey80"), #, linetype="solid", color="grey60",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(vjust = 3),
        legend.position = c(0.9,0.87),
        legend.title = element_blank(),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.background = element_rect(fill = "transparent", colour = NA)
        # legend.text = element_text(size = 11)
        # plot.margin = unit(c(0,0,0,0.5),"cm")
  ) +
  facet_grid(groups ~ cycl,
             labeller = labeller(cycl = cycl.names,
                                 groups = group.names))
p1


# annual = annual + sub-annual
ann_text1 <- data.frame(week = 21,
                        ss_text = 1.7,
                        groups = c("group1","group1","group2","group2","no_class","no_class"),
                        cycl = c("annual","no","annual","no","annual","no"),
                        notes = c("D: 8 sp, 1.9% BA, 11.0% TC\nT: 10 sp, 11.3% BA, 36.6% TC",
                                  "D: 3 sp, 1.8% BA, 6.3% TC\nT: 8 sp, 1.7% BA, 5.2% TC",
                                  "D: 7 sp, 3.3% BA, 24.6% TC\nT: 3 sp, 0.3% BA, 1.7% TC",
                                  "D: 5 sp, 3.5% BA, 17.5% TC\nT: 10 sp, 2.6% BA, 7.3% TC",
                                  "D: 4 sp, 0.25% BA, 2.0% TC\nT: 8 sp, 2.1% BA, 10.2% TC",
                                  "D: 54 sp, 57.7% BA, 38.6% TC\nT: 51 sp, 52.4% BA, 39.0% TC"))
ann_segm <- data.frame(x = c(19.5, 24),
                       y = c(0.35, 1.1),
                       xend = c(24, 26.5),
                       yend = c(1.1, 1.1),
                       groups = c("group1","group1"),
                       cycl = c("annual","annual"))
ann_segm2 <- data.frame(x = c(13, 11),
                       y = c(1.2, 1.2),
                       xend = c(18.5, 13),
                       yend = c(0.7, 1.2),
                       groups = c("no_class","no_class"),
                       cycl = c("no","no"))
ann_text2 <- data.frame(week = c(27.5, 1),
                        ss_text = c(0.98, 1.1),
                        groups = c("group1","no_class"),
                        cycl = c("annual","no"),
                        notes = c("P. macrocarpus\n6.6% BA\nperiodically sub-annual",
                                  "S. zenkeri\n16.3% BA\nsingular event"))

p2 <- p1 +
  geom_text(data = ann_text1,
                     aes(week, ss_text, label = notes),
                     size = 3.2,
                     hjust = 0,
            family = "Lato", color = "#22211d") +
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
            lineheight = .9,
            family = "Lato", color = "#22211d")

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

# head(standlevel_subsets)
#
# standlevel_subsets$Month <- ifelse(standlevel_subsets$week %in% c(1,2,3,4),1,
#                          ifelse(standlevel_subsets$week %in% c(5,6,7,8),2,
#                                 ifelse(standlevel_subsets$week %in% c(9,10,11,12),3,
#                                        ifelse(standlevel_subsets$week %in% c(13,14,15,16),4,
#                                               ifelse(standlevel_subsets$week %in% c(17,18,19,20),5,
#                                                      ifelse(standlevel_subsets$week %in% c(21,22,23,24),6,
#                                                             ifelse(standlevel_subsets$week %in% c(25,26,27,28),7,
#                                                                    ifelse(standlevel_subsets$week %in% c(29,30,31,32),8,
#                                                                           ifelse(standlevel_subsets$week %in% c(33,34,35,36),9,
#                                                                                  ifelse(standlevel_subsets$week %in% c(37,38,39,40),10,
#                                                                                         ifelse(standlevel_subsets$week %in% c(41,42,43,44),11,
#                                                                                                ifelse(standlevel_subsets$week %in% c(45,46,47,48),12,
#                                                                                                       NA))))))))))))
#
# standlevel_subsets_month <- standlevel_subsets %>%
#   group_by(Month, cycl, groups) %>%
#   summarise(dormancy_stand = mean(ss_dorm, na.rm = T),
#             turnover_stand = mean(ss_turn))
#
# climate.corr <- merge(climate, standlevel_subsets_month, by = c("Month"), all.x = TRUE)
#
# climate.corr <- climate.corr %>%
#   filter(cycl == "no",
#          groups == "group2")
#
#
# # # turnover
# # cor.test(climate.corr$insol_JR, climate.corr$turnover_stand, method = 'pearson')
# # cor.test(climate.corr$prec_JR, climate.corr$turnover_stand, method = 'pearson')
# # cor.test(climate.corr$tmax_JR, climate.corr$turnover_stand, method = 'pearson')
# #
# #
# # # dormancy
# # cor.test(climate.corr$insol_JR, climate.corr$dormancy_stand, method = 'pearson')
# # cor.test(climate.corr$prec_JR, climate.corr$dormancy_stand, method = 'pearson')
# # cor.test(climate.corr$tmax_JR, climate.corr$dormancy_stand, method = 'pearson')
#
# # # MODIS
# # VI_s$date <- as.Date(VI_s$doy, origin="2010-12-31")
# # VI_s$Month <- format(as.Date(VI_s$date), "%m")
# # modis <- VI_s %>%
# #   group_by(Month)%>%
# #   dplyr::summarise(EVIm = mean(EVI))
# # modis$Month <- as.numeric(modis$Month)
#
# climate.corr <- merge(climate.corr, modis, by = c("Month"), all.x = TRUE)
#
# # cor.test(climate.corr$EVIm, climate.corr$prec_all, method = 'pearson')
# # cor.test(climate.corr$EVIm, climate.corr$insol_all, method = 'pearson')
# # cor.test(climate.corr$EVIm, climate.corr$PAR_Ygb_Hauser, method = 'pearson')
# # cor.test(climate.corr$EVIm, climate.corr$tmax_JR, method = 'pearson')
#
# cor.test(climate.corr$EVIm, climate.corr$dormancy_stand, method = 'pearson')
# cor.test(climate.corr$EVIm, climate.corr$turnover_stand, method = 'pearson')




