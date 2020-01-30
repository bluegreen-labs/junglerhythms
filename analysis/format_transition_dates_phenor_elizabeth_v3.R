#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
# graphics.off()
#----- load required packages --------------------------------------------------#
library(zoo)
library(tidyverse)
library(circular)
library(ggplot2)
library(gridExtra)
library(scales)
#----- source required functions -----------------------------------------------#
source("R/event_length.R")
#-------------------------------------------------------------------------------#


#----------------------------------------------------------------------
#-------- input that you can change   ---------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------
minimum_events = 5    # species will only be included in this analysis
# when they have a minimum set of events
#----------------------------------------------------------------------



#----------------------------------------------------------------------
#--------   Phenology data - species correction Meise   ---------------
#----------------------------------------------------------------------
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
df <- df[which(df$value != 0),]
df$join_id <- paste0("R",df$image,"-",df$image_row)

metadata <- read.csv("data/phenology_archives_species_long_format_20190626.csv",
                     header = TRUE, sep = ",")
metadata$join_id <- paste(metadata$image,metadata$row, sep = "-")

# test merge the two tables based upon the unique join ID
data <- merge(df, metadata, by = c("join_id"), all.x = TRUE)
data$species_full <- paste(data$genus_Meise, data$species_Meise)

# remove column id.x and rename id.y to id (--> in id.y, empty ids are renamed to EK1, EK2, etc...)
data$id <- as.character(data$id.y)
data = data[,!(names(data) %in% c("id.x","id.y"))]

# remove rows with NA's in year -> individuals with 'no_data' in the archive
data <- data[!(is.na(data$year)),]


# #----------------------------------------------------------------------
# #-------- get the species list ----------------------------------------
# #----------------------------------------------------------------------
overview <- read.csv("data/species_meta_data.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)

# merge number of events to data, so you can filter on minimum number of events
overview <- overview[,(names(overview) %in% c("species_full",
                                              "deciduousness",
                                              "basal_area_site"))]

#--------------------------------
# merge with data to filter on deciduousness
data <- merge(data, overview, by = "species_full", all.x = TRUE)
#--------------------------------

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------- This is onset of event
#---------- turnover --------------------------------------
#---------- DECIDUOUS
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# use function event_length to get the starting week of an event
transition_dates_turnover_dec <- data %>%
  filter(phenophase == "leaf_turnover") %>%
  filter(grepl("deciduous",deciduousness)) %>%
  group_by(species_full, id) %>%
  do(event_length(.))

# weeks to degrees
transition_dates_turnover_dec <- transition_dates_turnover_dec %>%
  mutate(degree = (week_start-1) * 360/48) %>%# -1 so that week 1 in janurari is degree 0
  mutate(duration = (phenophase_length) * 360/48)


# summary of onset dates per species
onset_turnover_dec <- transition_dates_turnover_dec %>%
  group_by(species_full) %>%
  filter(length(week_start) >= minimum_events) %>%    # only species with a minimum of events
  summarise(
    mean_degree = mean.circular(
      circular(degree, units = "degrees")
    ),
    median_degree = median.circular(
      circular(degree, units = "degrees")
    ),
    lower = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees")
    )$mu.ci[1],
    upper = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees")
    )$mu.ci[2],
    nr_events = length(week_start),
    duration = mean(duration)
  )

# timing flushing added to match columns dormancy
onset_turnover_dec$flushing_degree <- NA

# merge with census and remove rows (species) not included in the Yangambi mixed forest census
onset_turnover_dec <- merge(onset_turnover_dec, overview, by = "species_full", all.x = TRUE)
onset_turnover_dec <- onset_turnover_dec[!(is.na(onset_turnover_dec$basal_area_site)),]

# sort by BA
# set y value by sorted BA number
onset_turnover_dec <- onset_turnover_dec %>%
  arrange(basal_area_site) %>%
  mutate(y_value = (1:length(species_full)))

# rescaling between 0 and 360 degrees
onset_turnover_dec$mean_rescaled <- ifelse(onset_turnover_dec$mean_degree < 0, onset_turnover_dec$mean_degree +360, onset_turnover_dec$mean_degree)
onset_turnover_dec$median_rescaled <- ifelse(onset_turnover_dec$median_degree < 0, onset_turnover_dec$median_degree +360, onset_turnover_dec$median_degree)
onset_turnover_dec$upper_rescaled <- ifelse (onset_turnover_dec$upper < 0, onset_turnover_dec$upper + 360, onset_turnover_dec$upper)
onset_turnover_dec$lower_rescaled <- ifelse (onset_turnover_dec$lower < 0, onset_turnover_dec$lower + 360, onset_turnover_dec$lower)
onset_turnover_dec$flushing_rescaled <- NA

# in order to plot the CI if CI goes over zero, two segments needed: lower - 360 and 0 - upper
onset_turnover_dec_CItwosegments <- onset_turnover_dec %>%
  filter(mean_rescaled < lower_rescaled | mean_rescaled > upper_rescaled)
onset_turnover_dec_CItwosegments$seg1_xend <- 360  # end of segment 1 --> lower:360
onset_turnover_dec_CItwosegments$seg2_x <- 0       # beginning of segment 2 --> 0:upper
# for CI not going over zero, the segment can just be lower - upper
onset_turnover_dec_CIonesegments <- onset_turnover_dec %>%
  filter(mean_rescaled > lower_rescaled & mean_rescaled < upper_rescaled)
#------------------------------------------------------------------------



#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------- This is onset of event
#---------- DORMANCY --------------------------------------
#---------- DECIDUOUS
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# use function event_length to get the starting week of an event
transition_dates_dormancy_dec <- data %>%
  filter(phenophase == "leaf_dormancy") %>%
  filter(grepl("deciduous",deciduousness)) %>%
  group_by(species_full, id) %>%
  do(event_length(.))

# weeks to degrees
transition_dates_dormancy_dec <- transition_dates_dormancy_dec %>%
  mutate(degree = (week_start-1) * 360/48) %>%# -1 so that week 1 in janurari is degree 0
  mutate(duration = (phenophase_length) * 360/48)

# summary of onset dates per species
onset_dormancy_dec <- transition_dates_dormancy_dec %>%
  group_by(species_full) %>%
  filter(length(week_start) >= minimum_events) %>%    # only species with a minimum of events
  summarise(
    mean_degree = mean.circular(
      circular(degree, units = "degrees")
    ),
    median_degree = median.circular(
      circular(degree, units = "degrees")
    ),
    lower = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees")
    )$mu.ci[1],
    upper = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees")
    )$mu.ci[2],
    nr_events = length(week_start),
    duration = mean(duration)
  )
# timing flushing = degree + duration + (7.5 = 1 week (1*360/48))
onset_dormancy_dec$flushing_degree <- onset_dormancy_dec$median_degree + onset_dormancy_dec$duration + 7.5

# merge with census and remove rows (species) not included in the Yangambi mixed forest census
onset_dormancy_dec <- merge(onset_dormancy_dec, overview, by = "species_full", all.x = TRUE)
onset_dormancy_dec <- onset_dormancy_dec[!(is.na(onset_dormancy_dec$basal_area_site)),]

# sort by BA
# set y value by sorted BA number
onset_dormancy_dec <- onset_dormancy_dec %>%
  arrange(basal_area_site) %>%
  mutate(y_value = (1:length(species_full) + max(onset_turnover_dec$y_value) + 10))

# rescaling between 0 and 360 degrees
onset_dormancy_dec$mean_rescaled <- ifelse(onset_dormancy_dec$mean_degree < 0, onset_dormancy_dec$mean_degree +360, onset_dormancy_dec$mean_degree)
onset_dormancy_dec$median_rescaled <- ifelse(onset_dormancy_dec$median_degree < 0, onset_dormancy_dec$median_degree +360, onset_dormancy_dec$median_degree)
onset_dormancy_dec$upper_rescaled <- ifelse (onset_dormancy_dec$upper < 0, onset_dormancy_dec$upper + 360, onset_dormancy_dec$upper)
onset_dormancy_dec$lower_rescaled <- ifelse (onset_dormancy_dec$lower < 0, onset_dormancy_dec$lower + 360, onset_dormancy_dec$lower)
onset_dormancy_dec$flushing_rescaled <- ifelse(onset_dormancy_dec$flushing_degree < 0, onset_dormancy_dec$flushing_degree +360, onset_dormancy_dec$flushing_degree)

# in order to plot the CI if CI goes over zero, two segments needed: lower - 360 and 0 - upper
onset_dormancy_dec_CItwosegments <- onset_dormancy_dec %>%
  filter(mean_rescaled < lower_rescaled | mean_rescaled > upper_rescaled)
onset_dormancy_dec_CItwosegments$seg1_xend <- 360  # end of segment 1 --> lower:360
onset_dormancy_dec_CItwosegments$seg2_x <- 0       # beginning of segment 2 --> 0:upper
# for CI not going over zero, the segment can just be lower - upper
onset_dormancy_dec_CIonesegments <- onset_dormancy_dec %>%
  filter(mean_rescaled > lower_rescaled & mean_rescaled < upper_rescaled)
#------------------------------------------------------------------------


# merge tunover and dormancy for deciduous species
onset_turnover_dec$phase <- "turnover"
onset_dormancy_dec$phase <- "dormancy"
onset_dec <- rbind(onset_turnover_dec, onset_dormancy_dec)
onset_turnover_dec_CIonesegments$phase <- "turnover"
onset_dormancy_dec_CIonesegments$phase <- "dormancy"
onset_dec_CIonesegments <- rbind(onset_turnover_dec_CIonesegments, onset_dormancy_dec_CIonesegments)
onset_turnover_dec_CItwosegments$phase <- "turnover"
onset_dormancy_dec_CItwosegments$phase <- "dormancy"
onset_dec_CItwosegments <- rbind(onset_turnover_dec_CItwosegments, onset_dormancy_dec_CItwosegments)

#
onset_dec_allphases <- onset_dec[,(names(onset_dec) %in% c("species_full",
                                                           "y_value",
                                                         "median_rescaled",
                                                         "phase"))]
onset_dec_flushing <- onset_dec[,(names(onset_dec) %in% c("species_full",
                                                          "y_value",
                                                "flushing_rescaled"
                                                ))]
colnames(onset_dec_flushing)[3] <- "median_rescaled"
onset_dec_flushing <- onset_dec_flushing[!(is.na(onset_dec_flushing$median_rescaled)),]
onset_dec_flushing$phase <- "flushing"

onset_dec_allphases <- rbind(onset_dec_allphases, onset_dec_flushing)

#------------------------------------------------------------------------
# FIGURE DECIDUOUS
#------------------------------------------------------------------------
sp_clust <- read.csv("data/cluster_test.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)
# merge with data to filter on deciduousness
onset_dec_allphases <- merge(onset_dec_allphases, sp_clust, by = "species_full", all.x = TRUE)
onset_dec_allphases$cluster_id <- ifelse(onset_dec_allphases$phase %in% 'dormancy', onset_dec_allphases$cluster_id, "0")
onset_dec_allphases$cluster_id <- as.character(onset_dec_allphases$cluster_id)

#--------------------------------
p_onset_dec <- ggplot(data = onset_dec) +
  annotate("rect", xmin = 330, xmax = 360, ymin = 0, ymax = max(onset_dec$y_value) , alpha = .2) + #Dec
  annotate("rect", xmin = 0, xmax = 60, ymin = 0, ymax = max(onset_dec$y_value) , alpha = .2) + # jan - feb
  annotate("rect", xmin = 150, xmax = 210, ymin = 0, ymax = max(onset_dec$y_value) , alpha = .2) + # jun-jul
  geom_segment(data = onset_dec_CIonesegments,
               aes(x = lower_rescaled,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = onset_dec_CItwosegments,
               aes(x = lower_rescaled,
                   xend = seg1_xend,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = onset_dec_CItwosegments,
               aes(x = seg2_x,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_point(data = onset_dec_allphases,
             aes(x = median_rescaled, y = y_value, shape = phase, colour = cluster_id)) + #
  # geom_point(aes(x = mean_rescaled, y = y_value, col = phase),
  #            shape = 4) +
  scale_shape_manual(values = c(19,4,1),
                     name = "onset of:",
                     labels = c("dormancy", "flushing","turnover")) +
  scale_colour_manual(values = c("black","red","green","blue","purple","darkblue","orange")) +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0,359,30),
                     labels = month.abb) +
  scale_y_continuous(limits = c(-20,max(onset_dec$y_value))) +  #staring y-axis at '-20' so that points are squeezed in the center
  # minor_breaks = seq(1, 48, 1), breaks = seq(1, 48, 1)) +
  coord_polar() +
  labs(x="",
       y="",
       title = paste("Deciduous")) +
  theme(panel.grid.major.x = element_line(colour = "grey75",
                                          size = 0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),# element_line(colour = "grey90", size = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11),
        strip.text = element_text(face = "italic", size = 13),
        strip.background = element_rect(fill="white"),
        legend.position = "bottom",
        # legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.key = element_rect(fill = "white"),
        plot.margin=unit(c(0,0,0,0),"cm")
  )
p_onset_dec
#------------------------------------------------------------------------

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------- This is onset of event
#---------- turnover --------------------------------------
#---------- EVERGREEN
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# use function event_length to get the starting week of an event
transition_dates_turnover_ever <- data %>%
  filter(phenophase == "leaf_turnover") %>%
  filter(grepl("evergreen",deciduousness)) %>%
  group_by(species_full, id) %>%
  do(event_length(.))

# weeks to degrees
transition_dates_turnover_ever <- transition_dates_turnover_ever %>%
  mutate(degree = (week_start-1) * 360/48) %>%# -1 so that week 1 in janurari is degree 0
  mutate(duration = (phenophase_length) * 360/48)

# summary of onset dates per species
onset_turnover_ever <- transition_dates_turnover_ever %>%
  group_by(species_full) %>%
  filter(length(week_start) >= minimum_events) %>%    # only species with a minimum of events
  summarise(
    mean_degree = mean.circular(
      circular(degree, units = "degrees")
    ),
    median_degree = median.circular(
      circular(degree, units = "degrees")
    ),
    lower = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees")
    )$mu.ci[1],
    upper = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees")
    )$mu.ci[2],
    nr_events = length(week_start),
    duration = mean(duration)
  )

# timing flushing to link to column names dormancy
onset_turnover_ever$flushing_degree <- NA

# merge with census and remove rows (species) not included in the Yangambi mixed forest census
onset_turnover_ever <- merge(onset_turnover_ever, overview, by = "species_full", all.x = TRUE)
onset_turnover_ever <- onset_turnover_ever[!(is.na(onset_turnover_ever$basal_area_site)),]

# sort by BA
# set y value by sorted BA number
onset_turnover_ever <- onset_turnover_ever %>%
  arrange(basal_area_site) %>%
  mutate(y_value = (1:length(species_full)))

# rescaling between 0 and 360 degrees
onset_turnover_ever$mean_rescaled <- ifelse(onset_turnover_ever$mean_degree < 0, onset_turnover_ever$mean_degree +360, onset_turnover_ever$mean_degree)
onset_turnover_ever$median_rescaled <- ifelse(onset_turnover_ever$median_degree < 0, onset_turnover_ever$median_degree +360, onset_turnover_ever$median_degree)
onset_turnover_ever$upper_rescaled <- ifelse (onset_turnover_ever$upper < 0, onset_turnover_ever$upper + 360, onset_turnover_ever$upper)
onset_turnover_ever$lower_rescaled <- ifelse (onset_turnover_ever$lower < 0, onset_turnover_ever$lower + 360, onset_turnover_ever$lower)
onset_turnover_ever$flushing_rescaled <- NA

# in order to plot the CI if CI goes over zero, two segments needed: lower - 360 and 0 - upper
onset_turnover_ever_CItwosegments <- onset_turnover_ever %>%
  filter(mean_rescaled < lower_rescaled | mean_rescaled > upper_rescaled)
onset_turnover_ever_CItwosegments$seg1_xend <- 360  # end of segment 1 --> lower:360
onset_turnover_ever_CItwosegments$seg2_x <- 0       # beginning of segment 2 --> 0:upper
# for CI not going over zero, the segment can just be lower - upper
onset_turnover_ever_CIonesegments <- onset_turnover_ever %>%
  filter(mean_rescaled > lower_rescaled & mean_rescaled < upper_rescaled)
#------------------------------------------------------------------------



#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------- This is onset of event
#---------- DORMANCY --------------------------------------
#---------- EVERGREEN
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# use function event_length to get the starting week of an event
transition_dates_dormancy_ever <- data %>%
  filter(phenophase == "leaf_dormancy") %>%
  filter(grepl("evergreen",deciduousness)) %>%
  group_by(species_full, id) %>%
  do(event_length(.))

# weeks to degrees
transition_dates_dormancy_ever <- transition_dates_dormancy_ever %>%
  mutate(degree = (week_start-1) * 360/48) %>%# -1 so that week 1 in janurari is degree 0
  mutate(duration = (phenophase_length) * 360/48)

# summary of onset dates per species
onset_dormancy_ever <- transition_dates_dormancy_ever %>%
  group_by(species_full) %>%
  filter(length(week_start) >= minimum_events) %>%    # only species with a minimum of events
  summarise(
    mean_degree = mean.circular(
      circular(degree, units = "degrees")
    ),
    median_degree = median.circular(
      circular(degree, units = "degrees")
    ),
    lower = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees")
    )$mu.ci[1],
    upper = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees")
    )$mu.ci[2],
    nr_events = length(week_start),
    duration = mean(duration)
  )
# timing flushing = degree + duration + (7.5 = 1 week (1*360/48))
onset_dormancy_ever$flushing_degree <- onset_dormancy_ever$median_degree + onset_dormancy_ever$duration + 7.5

# merge with census and remove rows (species) not included in the Yangambi mixed forest census
onset_dormancy_ever <- merge(onset_dormancy_ever, overview, by = "species_full", all.x = TRUE)
onset_dormancy_ever <- onset_dormancy_ever[!(is.na(onset_dormancy_ever$basal_area_site)),]

# sort by BA
# set y value by sorted BA number
onset_dormancy_ever <- onset_dormancy_ever %>%
  arrange(basal_area_site) %>%
  mutate(y_value = (1:length(species_full) + max(onset_turnover_ever$y_value) + 8))

# rescaling between 0 and 360 degrees
onset_dormancy_ever$mean_rescaled <- ifelse(onset_dormancy_ever$mean_degree < 0, onset_dormancy_ever$mean_degree +360, onset_dormancy_ever$mean_degree)
onset_dormancy_ever$median_rescaled <- ifelse(onset_dormancy_ever$median_degree < 0, onset_dormancy_ever$median_degree +360, onset_dormancy_ever$median_degree)
onset_dormancy_ever$upper_rescaled <- ifelse (onset_dormancy_ever$upper < 0, onset_dormancy_ever$upper + 360, onset_dormancy_ever$upper)
onset_dormancy_ever$lower_rescaled <- ifelse (onset_dormancy_ever$lower < 0, onset_dormancy_ever$lower + 360, onset_dormancy_ever$lower)
onset_dormancy_ever$flushing_rescaled <- ifelse(onset_dormancy_ever$flushing_degree < 0, onset_dormancy_ever$flushing_degree +360, onset_dormancy_ever$flushing_degree)

# in order to plot the CI if CI goes over zero, two segments needed: lower - 360 and 0 - upper
onset_dormancy_ever_CItwosegments <- onset_dormancy_ever %>%
  filter(mean_rescaled < lower_rescaled | mean_rescaled > upper_rescaled)
onset_dormancy_ever_CItwosegments$seg1_xend <- 360  # end of segment 1 --> lower:360
onset_dormancy_ever_CItwosegments$seg2_x <- 0       # beginning of segment 2 --> 0:upper
# for CI not going over zero, the segment can just be lower - upper
onset_dormancy_ever_CIonesegments <- onset_dormancy_ever %>%
  filter(mean_rescaled > lower_rescaled & mean_rescaled < upper_rescaled)
#------------------------------------------------------------------------


# merge tunover and dormancy for deciduous species
onset_turnover_ever$phase <- "turnover"
onset_dormancy_ever$phase <- "dormancy"
onset_ever <- rbind(onset_turnover_ever, onset_dormancy_ever)
onset_turnover_ever_CIonesegments$phase <- "turnover"
onset_dormancy_ever_CIonesegments$phase <- "dormancy"
onset_ever_CIonesegments <- rbind(onset_turnover_ever_CIonesegments, onset_dormancy_ever_CIonesegments)
onset_turnover_ever_CItwosegments$phase <- "turnover"
onset_dormancy_ever_CItwosegments$phase <- "dormancy"
onset_ever_CItwosegments <- rbind(onset_turnover_ever_CItwosegments, onset_dormancy_ever_CItwosegments)

#
onset_ever_allphases <- onset_ever[,(names(onset_ever) %in% c("y_value",
                                                           "median_rescaled",
                                                           "phase"))]
onset_ever_flushing <- onset_ever[,(names(onset_ever) %in% c("y_value",
                                                          "flushing_rescaled"
))]
colnames(onset_ever_flushing)[2] <- "median_rescaled"
onset_ever_flushing <- onset_ever_flushing[!(is.na(onset_ever_flushing$median_rescaled)),]
onset_ever_flushing$phase <- "flushing"

onset_ever_allphases <- rbind(onset_ever_allphases, onset_ever_flushing)

#------------------------------------------------------------------------
# FIGURE EVERGREEN
#------------------------------------------------------------------------
p_onset_ever <- ggplot(data = onset_ever) +
  annotate("rect", xmin = 330, xmax = 360, ymin = 0, ymax = max(onset_ever$y_value) , alpha = .2) + #Dec
  annotate("rect", xmin = 0, xmax = 60, ymin = 0, ymax = max(onset_ever$y_value) , alpha = .2) + # jan - feb
  annotate("rect", xmin = 150, xmax = 210, ymin = 0, ymax = max(onset_ever$y_value) , alpha = .2) + # jun-jul
  geom_segment(data = onset_ever_CIonesegments,
               aes(x = lower_rescaled,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = onset_ever_CItwosegments,
               aes(x = lower_rescaled,
                   xend = seg1_xend,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = onset_ever_CItwosegments,
               aes(x = seg2_x,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_point(data = onset_ever_allphases,
             aes(x = median_rescaled, y = y_value, shape = phase)) +
  scale_shape_manual(values = c(19,4,1),
                     name = "onset of",
                     labels = c("dormancy", "flushing","turnover")) +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0,359,30),
                     labels = month.abb) +
  scale_y_continuous(limits = c(-20,max(onset_ever$y_value))) +  #staring y-axis at '-20' so that points are squeezed in the center
  # minor_breaks = seq(1, 48, 1), breaks = seq(1, 48, 1)) +
  coord_polar() +
  labs(x="",
       y="",
       title = paste("Evergreen")) +
  theme(panel.grid.major.x = element_line(colour = "grey75",
                                          size = 0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),# element_line(colour = "grey90", size = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11),
        strip.text = element_text(face = "italic", size = 13),
        strip.background = element_rect(fill="white"),
        legend.position = "none",
        plot.margin=unit(c(0.2,0,0,0),"cm")
  )
#------------------------------------------------------------------------

p_onset <- grid.arrange(p_onset_ever, p_onset_dec, heights = c(0.9,1))


pdf("~/Desktop/figure2_onset.pdf",4.5,10)
plot(p_onset)
dev.off()


