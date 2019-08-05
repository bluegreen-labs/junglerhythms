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
data = data[,!(names(data) %in% "id.x")]
data <- data %>%
  rename("id" = id.y)
data$id <- as.character(data$id)


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
                                              "mating_system",
                                              "diaspora",
                                              "basal_area_site"))]
overview$mating_system <- ifelse(overview$mating_system %in% "andromonoecious", "monoecious", overview$mating_system)
overview$mating_system <- ifelse(overview$mating_system %in% "monoecious or hermaphrodite", "hermaphrodite", overview$mating_system)
overview$mating_system <- ifelse(overview$mating_system %in% "monoecious or dioecious", "unknown", overview$mating_system)

#--------------------------------
# merge with data to filter on mating_system of diaspora
data <- merge(data, overview, by = "species_full", all.x = TRUE)
#--------------------------------

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------- This is onset of event
#---------- flowering --------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# use function event_length to get the starting week of an event
transition_dates_flower <- data %>%
  filter(phenophase == "fruit") %>% #
  # filter(grepl("evergreen",deciduousness)) %>%
  group_by(species_full, id) %>%
  do(event_length(.))

# weeks to degrees
transition_dates_flower <- transition_dates_flower %>%
  mutate(degree = (week_start-1) * 360/48)# -1 so that week 1 in janurari is degree 0


# summary of onset dates per species
onset_flower <- transition_dates_flower %>%
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
    nr_events = length(week_start)
  )

# merge with census and remove rows (species) not included in the Yangambi mixed forest census
onset_flower <- merge(onset_flower, overview, by = "species_full", all.x = TRUE)
onset_flower <- onset_flower[!(is.na(onset_flower$basal_area_site)),]

# sort by BA
# set y value by sorted BA number
onset_flower <- onset_flower %>%
  arrange(basal_area_site) %>%
  mutate(y_value = (1:length(species_full)))

# rescaling between 0 and 360 degrees
onset_flower$mean_rescaled <- ifelse(onset_flower$mean_degree < 0, onset_flower$mean_degree +360, onset_flower$mean_degree)
onset_flower$median_rescaled <- ifelse(onset_flower$median_degree < 0, onset_flower$median_degree +360, onset_flower$median_degree)
onset_flower$upper_rescaled <- ifelse (onset_flower$upper < 0, onset_flower$upper + 360, onset_flower$upper)
onset_flower$lower_rescaled <- ifelse (onset_flower$lower < 0, onset_flower$lower + 360, onset_flower$lower)

# in order to plot the CI if CI goes over zero, two segments needed: lower - 360 and 0 - upper
onset_flower_CItwosegments <- onset_flower %>%
  filter(mean_rescaled < lower_rescaled | mean_rescaled > upper_rescaled)
onset_flower_CItwosegments$seg1_xend <- 360  # end of segment 1 --> lower:360
onset_flower_CItwosegments$seg2_x <- 0       # beginning of segment 2 --> 0:upper
# for CI not going over zero, the segment can just be lower - upper
onset_flower_CIonesegments <- onset_flower %>%
  filter(mean_rescaled > lower_rescaled & mean_rescaled < upper_rescaled)
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# FIGURE Flowering
#------------------------------------------------------------------------
p_onset_flower <- ggplot(data = onset_flower) +
  annotate("rect", xmin = 330, xmax = 360, ymin = 0, ymax = max(onset_flower$y_value) , alpha = .2) + #Dec
  annotate("rect", xmin = 0, xmax = 60, ymin = 0, ymax = max(onset_flower$y_value) , alpha = .2) + # jan - feb
  annotate("rect", xmin = 150, xmax = 210, ymin = 0, ymax = max(onset_flower$y_value) , alpha = .2) + # jun-jul
  geom_segment(data = onset_flower_CIonesegments,
               aes(x = lower_rescaled,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = onset_flower_CItwosegments,
               aes(x = lower_rescaled,
                   xend = seg1_xend,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = onset_flower_CItwosegments,
               aes(x = seg2_x,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_point(data = onset_flower,
             aes(x = median_rescaled, y = y_value, color = mating_system)) +
  scale_color_manual(values = c("red","blue","purple","black")) + #, name = "onset of" ,labels = c("dormancy", "flushing","turnover")
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0,359,30),
                     labels = month.abb) +
  scale_y_continuous(limits = c(-20,max(onset_flower$y_value))) +  #staring y-axis at '-20' so that points are squeezed in the center
  # minor_breaks = seq(1, 48, 1), breaks = seq(1, 48, 1)) +
  coord_polar() +
  labs(x="",
       y="",
       title = paste("Onset flowering")) +
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
        plot.margin=unit(c(0.2,0,0,0),"cm")
  )
# p_onset_flower
#------------------------------------------------------------------------
p_onset_flower_2 <- ggplot(data = onset_flower) +
  annotate("rect", xmin = 330, xmax = 360, ymin = 0, ymax = max(onset_flower$y_value) , alpha = .2) + #Dec
  annotate("rect", xmin = 0, xmax = 60, ymin = 0, ymax = max(onset_flower$y_value) , alpha = .2) + # jan - feb
  annotate("rect", xmin = 150, xmax = 210, ymin = 0, ymax = max(onset_flower$y_value) , alpha = .2) + # jun-jul
  geom_segment(data = onset_flower_CIonesegments,
               aes(x = lower_rescaled,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = onset_flower_CItwosegments,
               aes(x = lower_rescaled,
                   xend = seg1_xend,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_segment(data = onset_flower_CItwosegments,
               aes(x = seg2_x,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey60") +
  geom_point(data = onset_flower,
             aes(x = median_rescaled, y = y_value, color = diaspora)) +
  scale_color_manual(values = c("red","blue","black","purple")) + #, name = "onset of" ,labels = c("dormancy", "flushing","turnover")
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0,359,30),
                     labels = month.abb) +
  scale_y_continuous(limits = c(-20,max(onset_flower$y_value))) +  #staring y-axis at '-20' so that points are squeezed in the center
  # minor_breaks = seq(1, 48, 1), breaks = seq(1, 48, 1)) +
  coord_polar() +
  labs(x="",
       y="",
       title = paste("Onset flowering")) +
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
        plot.margin=unit(c(0.2,0,0,0),"cm")
  )
# p_onset_flower_2

p_onset <- grid.arrange(p_onset_flower, p_onset_flower_2, heights = c(1,1))


pdf("~/Desktop/flowering.pdf",6,10)
plot(p_onset)
dev.off()
