#----- reset your R session. ------------------------------------------
rm(list=ls())
graphics.off()
#----- load required packages -----------------------------------------
library(zoo)
library(tidyverse)
library(circular)
library(ggplot2)
library(gridExtra)
library(scales)
library(showtext)
font_add_google(
  "Lato",
  regular.wt = 300,
  bold.wt = 700)
#----- source required functions --------------------------------------
source("R/event_length.R")
source("R/timeline_gap_fill.R")
#----------------------------------------------------------------------

# Suppress summarise info
# because `summarise()` ungrouping output (override with `.groups` argument)
# comes up a lot in the consol, which is annoying
options(dplyr.summarise.inform = FALSE)

#----------------------------------------------------------------------
#--------   Phenology data --------------------------------------------
#----------------------------------------------------------------------
data <- readRDS("data/jungle_rhythms_data_cleaned.rds")
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#-------- get the species list ----------------------------------------
#----------------------------------------------------------------------
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
#-- get full range timelines ------------------------------------------
#----------------------------------------------------------------------
species_list <- overview$species_full
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
data_timeline <- merge(data_timeline, overview, by = "species_full", all.x = TRUE)


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# calculate median data + CI of onset phenophases
# function event_length() gives start and end timing of each event
# this function only analyses 1 phenophase at a time
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# use function event_length to get the starting week of an event
transition_dates_turnover <- event_length(data = data_timeline,
                                          species_name = species_list,
                                          pheno = "leaf_turnover")
transition_dates_dormancy <- event_length(data = data_timeline,
                                          species_name = species_list,
                                          pheno = "leaf_dormancy")
transition_dates <- rbind(transition_dates_turnover, transition_dates_dormancy)


# weeks to degrees
transition_dates <- transition_dates %>%
  filter(!is.na(year_start)) %>%
  mutate(degree_start = (week_start-1) * 360/48) # -1 so that week 1 in janurari is degree 0

# summary of onset dates per species
onset <- transition_dates %>%
  group_by(species_full, phenophase) %>%
  dplyr::summarise(
    mean_degree = mean.circular(
      circular(degree_start, units = "degrees")),
    median_degree = median.circular(
      circular(degree_start, units = "degrees")),
    lower = mle.vonmises.bootstrap.ci(
      circular(degree_start, units = "degrees"))$mu.ci[1],
    upper = mle.vonmises.bootstrap.ci(
      circular(degree_start, units = "degrees"))$mu.ci[2],
    nr_events = length(week_start))


#--------------------------------------------------------------------------
#- rescaling and segmenting CIs for figure
#--------------------------------------------------------------------------
# merge with overview
# sd_intrasp will be used
onset <- merge(onset, overview, by = "species_full", all.x = TRUE)
# median onset in degrees are given from -180 to 180
# lower - upper if going over zero also from -180 to 180; if not, between 0 and 360
# -> rescale all between 0 and 360 degrees
onset$median_rescaled <- ifelse(onset$median_degree < 0, onset$median_degree +360, onset$median_degree)
onset$upper_rescaled <- ifelse (onset$upper < 0, onset$upper + 360, onset$upper)
onset$lower_rescaled <- ifelse (onset$lower < 0, onset$lower + 360, onset$lower)
# give year-round uncertainty to species with NA sd_intrasp
onset$lower_rescaled <- ifelse(onset$phenophase == "leaf_turnover" & is.na(onset$sd_intrasp_onset_leaf_turnover_weeks), 0, onset$lower_rescaled)
onset$upper_rescaled <- ifelse(onset$phenophase == "leaf_turnover" & is.na(onset$sd_intrasp_onset_leaf_turnover_weeks), 360, onset$upper_rescaled)
onset$lower_rescaled <- ifelse(onset$phenophase == "leaf_dormancy" & is.na(onset$sd_intrasp_onset_leaf_dormancy_weeks), 0, onset$lower_rescaled)
onset$upper_rescaled <- ifelse(onset$phenophase == "leaf_dormancy" & is.na(onset$sd_intrasp_onset_leaf_dormancy_weeks), 360, onset$upper_rescaled)
#------------------------------------------------------------------------
# to sort by intra-species variability
# set y value by sorted SD number
# needs to be done for evergreen/deciduous and turnover/dormancy seperately
onset$sd_intrasp_onset_leaf_dormancy_weeks <- ifelse(is.na(onset$sd_intrasp_onset_leaf_dormancy_weeks),52, onset$sd_intrasp_onset_leaf_dormancy_weeks)
onset$sd_intrasp_onset_leaf_turnover_weeks <- ifelse(is.na(onset$sd_intrasp_onset_leaf_turnover_weeks),52, onset$sd_intrasp_onset_leaf_turnover_weeks)

onset_evergreen_turn <- onset %>%
  filter(grepl("evergreen",deciduousness),
         phenophase == "leaf_turnover") %>%
  arrange(desc(sd_intrasp_onset_leaf_turnover_weeks)) %>%
  mutate(y_value = ((1:length(species_full)/length(species_full))), # divided by length(species_full) so all on scale 0-1
         deciduousness = "1evergreen")
onset_evergreen_dorm <- onset %>%
  filter(grepl("evergreen",deciduousness),
         phenophase == "leaf_dormancy") %>%
  arrange(desc(sd_intrasp_onset_leaf_dormancy_weeks)) %>%
  mutate(y_value = ((1:length(species_full)/length(species_full))),
         deciduousness = "1evergreen")
onset_deciduous_turn <- onset %>%
  filter(grepl("deciduous",deciduousness),
         phenophase == "leaf_turnover") %>%
  arrange(desc(sd_intrasp_onset_leaf_turnover_weeks)) %>%
  mutate(y_value = ((1:length(species_full)/length(species_full))),
         deciduousness = "deciduous")
onset_deciduous_dorm <- onset %>%
  filter(grepl("deciduous",deciduousness),
         phenophase == "leaf_dormancy") %>%
  arrange(desc(sd_intrasp_onset_leaf_dormancy_weeks)) %>%
  mutate(y_value = ((1:length(species_full)/length(species_full))),
         deciduousness = "deciduous")
onset <- rbind(onset_evergreen_turn, onset_evergreen_dorm, onset_deciduous_turn, onset_deciduous_dorm)
#------------------------------------------------------------------------
# to plot over the 0-360 line, segment CI when necessary: lower - 360 and 0 - upper
onset_CItwosegments <- onset %>%
  filter(median_rescaled < lower_rescaled | median_rescaled > upper_rescaled)
onset_CItwosegments$seg1_xend <- 360  # end of segment 1 --> lower:360
onset_CItwosegments$seg2_x <- 0       # beginning of segment 2 --> 0:upper
# for CI not going over zero, the segment can just be lower - upper
onset_CIonesegments <- onset %>%
  filter(median_rescaled >= lower_rescaled & median_rescaled <= upper_rescaled)


#----------------------------------------------------------------------
#----------------------------------------------------------------------
# manuscript - figure 2
#----------------------------------------------------------------------
#----------------------------------------------------------------------
dec_labels <- c("1evergreen" = "evergreen", "deciduous" = "deciduous")
pheno_labels <- c("leaf_turnover" = "turnover", "leaf_dormancy" = "dormancy")

p1 <- ggplot() +
  annotate("rect", xmin = 330, xmax = 360, ymin = -0.25, ymax = 1, alpha = .2) + #Dec
  annotate("rect", xmin = 0, xmax = 60, ymin = -0.25, ymax = 1, alpha = .2) + # jan - feb
  annotate("rect", xmin = 150, xmax = 210, ymin = -0.25, ymax = 1, alpha = .2) + # jun-jul

  annotate("text", x = 0, y = -0.1, label = "LD", col = "grey50") +
  annotate("text", x = 90, y = -0.1, label = "SW", col = "grey50") +
  annotate("text", x = 180, y = -0.1, label = "SD", col = "grey50") +
  annotate("text", x = 270, y = -0.1, label = "LW", col = "grey50") +

  geom_segment(data = onset_CIonesegments,
               aes(x = lower_rescaled,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = onset_CItwosegments,
               aes(x = lower_rescaled,
                   xend = seg1_xend,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = onset_CItwosegments,
               aes(x = seg2_x,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_point(data = onset,
             aes(x = median_rescaled, y = y_value, shape = phenophase)) +
  scale_shape_manual(values = c(19,1)) +

  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0,359,30),
                     labels = month.abb) +
  scale_y_continuous(limits = c(-0.25,1)) +  #staring y-axis at '-0.25' so that points are squeezed in the center
  coord_polar() +
  labs(x="",
       y="",
       title = "") +
  theme(text = element_text(family = "Lato", color = "#22211d"),
        panel.grid.major.x = element_line(colour = "grey75",
                                          size = 0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        plot.title = element_blank(),
        legend.position = "none",
        plot.margin=unit(c(0,0,0,0),"cm")) +
  facet_grid(phenophase ~ deciduousness,
             switch = "y",
             labeller = labeller(deciduousness = dec_labels, phenophase = pheno_labels))

p1
