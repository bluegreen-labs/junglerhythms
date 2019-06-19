#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
# graphics.off()
#----- load required packages --------------------------------------------------#
library(zoo)
library(tidyverse)
library(circular)
library(scales)
#----- source required functions -----------------------------------------------#
source("R/event_length.R")
#-------------------------------------------------------------------------------#


#----------------------------------------------------------------------
#-------- input that you can change   ---------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------
phenophase_selected = "leaf_turnover"
title_name = "(b) Canopy turnover"
minimum_events = 5    # species will only be included in this analysis
# when they have a minimum set of events
#----------------------------------------------------------------------



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
data <- merge(df, metadata, by = c("join_id"), all.x = TRUE)
data$species_full <- paste(data$genus_Meise, data$species_Meise)

# remove column id.x and rename id.y to id (--> in id.y, empty ids are renamed to EK1, EK2, etc...)
data = data[,!(names(data) %in% "id.x")]
data <- data %>%
  rename("id" = id.y)
data$id <- as.character(data$id)

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
  filter(grepl("MIX",Plot)) %>%
  group_by(species_full) %>%
  dplyr::summarise(basal_area_site = sum(basal_area))
#----------------------------------------------------------------------




#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------- This is onset of event    --------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# use function event_length to get the starting week of an event
transition_dates <- data %>%
  filter(phenophase == phenophase_selected) %>%
  group_by(species_full, id, phenophase) %>%
  do(event_length(.))

# weeks to degrees
transition_dates <- transition_dates %>%
  mutate(degree = (week_start-1) * 360/48) # -1 so that week 1 in janurari is degree 0

# summary of onset dates per species
onset <- transition_dates %>%
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
onset <- merge(onset, census_site, by = "species_full", all.x = TRUE)
onset <- onset[!(is.na(onset$basal_area_site)),]

# sort by BA
# set y value by sorted BA number
onset <- onset %>%
  arrange(basal_area_site) %>%
  mutate(y_value = (1:length(species_full)))

# rescaling between 0 and 360 degrees
onset$mean_rescaled <- ifelse(onset$mean_degree < 0, onset$mean_degree +360, onset$mean_degree)
onset$median_rescaled <- ifelse(onset$median_degree < 0, onset$median_degree +360, onset$median_degree)
onset$upper_rescaled <- ifelse (onset$upper < 0, onset$upper + 360, onset$upper)
onset$lower_rescaled <- ifelse (onset$lower < 0, onset$lower + 360, onset$lower)

# in order to plot the CI if CI goes over zero, two segments needed: lower - 360 and 0 - upper
onset_CItwosegments <- onset %>%
  filter(mean_rescaled < lower_rescaled | mean_rescaled > upper_rescaled)
onset_CItwosegments$seg1_xend <- 360  # end of segment 1 --> lower:360
onset_CItwosegments$seg2_x <- 0       # beginning of segment 2 --> 0:upper
# for CI not going over zero, the segment can just be lower - upper
onset_CIonesegments <- onset %>%
  filter(mean_rescaled > lower_rescaled & mean_rescaled < upper_rescaled)

# for the species we have in the onset df, get the y_value
# link to the original transition_dates, so that we have all datapoints of a species, linked to y_value
overview <- onset[,(names(onset) %in% c("species_full","y_value"))]
overview <- inner_join(transition_dates, overview, by = c("species_full"))
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#----- calculate pairwise distances
#----- bootstrapped 95% confidence interval around mean (cross)
#------------------------------------------------------------------------
onset$median_rad <- rad(onset$median_degree)
a <- as.matrix(dist(onset$median_rad), labels = TRUE)
rownames(a) <- onset[['species_full']]
colnames(a) <- rownames(a)
b <- deg(a)
c <- ifelse(b > 180, 360-b, b)
d <- as.data.frame(rowMeans(c))
d$species_full <- rownames(d)
colnames(d)[1] <- "mean_distance_onset_weeks"


# # Ward Hierarchical Clustering
# # d <- dist(mydata, method = "euclidean") # distance matrix
# a <- dist(onset$median_rad, method = "euclidean")
# b <- deg(a)
# # c <- ifelse(b > 180, 360-b, b)
# fit <- hclust(b, method="ward")
# plot(fit) # display dendogram
# groups <- cutree(fit, k=4) # cut tree into 5 clusters
# # draw dendogram with red borders around the 5 clusters
# rect.hclust(fit, k=4, border="red")
#
# groups.cluster <- as.data.frame(groups)
# species_dates = onset[,(names(onset) %in% c("species_full",
#                                        "median_degree"))]
# species_groups <- cbind(species_dates, groups.cluster)
#
# plot(species_groups$groups, species_groups$median_degree)

#------------------------------------------------------------------------
#----- circular plot of median onset dates of event
#----- bootstrapped 95% confidence interval around mean (cross)
#------------------------------------------------------------------------
p_onset <- ggplot(data = onset) +
  geom_segment(data = onset_CIonesegments,
               aes(x = lower_rescaled,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey45") +
  geom_segment(data = onset_CItwosegments,
               aes(x = lower_rescaled,
                   xend = seg1_xend,
                   y = y_value,
                   yend = y_value),
               color = "grey45") +
  geom_segment(data = onset_CItwosegments,
               aes(x = seg2_x,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value),
               color = "grey45") +
  # geom_point(data = overview,
  #            aes(x = degree, y = y_value),
  #            color = "grey45",
  #            shape = 4) +
  geom_point(aes(x = median_rescaled, y = y_value)) +
  geom_point(aes(x = mean_rescaled, y = y_value),
             shape = 4) +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0,359,30),
                     labels = month.abb) +
  scale_y_continuous(limits = c(-20,length(onset$median_rescaled))) +  #staring y-axis at '-20' so that points are squeezed in the center
                     # minor_breaks = seq(1, 48, 1), breaks = seq(1, 48, 1)) +
  annotate("rect", xmin = 330, xmax = 360, ymin = 0, ymax = length(onset$median_rescaled) , alpha = .2) + #Dec
  annotate("rect", xmin = 0, xmax = 60, ymin = 0, ymax = length(onset$median_rescaled) , alpha = .2) + # jan - feb
  annotate("rect", xmin = 150, xmax = 210, ymin = 0, ymax = length(onset$median_rescaled) , alpha = .2) + # jun-jul
  coord_polar() +
  labs(x="",
       y="",
       title = title_name) +
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
        # legend.position = "right",
        plot.margin=unit(c(1,0,0,0),"cm")
  )
#------------------------------------------------------------------------

# pdf("~/Desktop/leaf_turnover_onset.pdf",7,7)
# plot(p_onset)
# dev.off()



# #--------------------------------------------------------------------------
# #--------------------------------------------------------------------------
# #--------- This is mean event time, not onset of event    -----------------
# #--------------------------------------------------------------------------
# #--------------------------------------------------------------------------
# #--------------------------------------------------------------------------
# test <- data %>%
#   filter(phenophase == "leaf_turnover",
#          value == 1) %>%
#   # filter(species_full == "Prioria balsamifera") %>%
#   mutate(degree = week * 360/48)
#
#
# bla <- test %>%
#   group_by(species_full) %>%
#   filter(length(value) > 10) %>%     # we only look at species with more than 10 events
#   summarise(
#     mean_degree = mean.circular(
#       circular(degree, units = "degrees")
#     ),
#     lower = mle.vonmises.bootstrap.ci(
#       circular(degree, units = "degrees")
#       )$mu.ci[1],
#     upper = mle.vonmises.bootstrap.ci(
#       circular(degree, units = "degrees")
#     )$mu.ci[2]
#     )
#
# # merge with census and remove rows (species) not included in the Yangambi mixed forest census
# bla <- merge(bla, census, by = "species_full", all.x = TRUE)
# bla <- bla[!(is.na(bla$BAperc)),]
#
# # sort by BA
# # number by sorted BA
# # set y value by sorted BA number
# bla <- bla %>%
#   arrange(BA) %>%
#   mutate(y_value = 1:length(species_full))
#
# # rescaling between 0 and 360 degrees
# bla$mean_rescaled <- ifelse(bla$mean_degree < 0, bla$mean_degree +360, bla$mean_degree)
# bla$upper_rescaled <- ifelse (bla$upper < 0, bla$upper + 360, bla$upper)
# bla$upper_rescaled <- circular(bla$upper_rescaled, units = "degrees")
# bla$lower_rescaled <- ifelse (bla$lower < 0, bla$lower + 360, bla$lower)
# bla$lower_rescaled <- circular(bla$lower_rescaled, units = "degrees")
#
# bla_subset <- bla %>%
#   filter(mean_rescaled < lower_rescaled | mean_rescaled > upper_rescaled)
# bla_subset$seg1_x <- 0
# bla_subset$seg2_x <- 360
#
# bla_ok <- bla %>%
#   filter(mean_rescaled > lower_rescaled & mean_rescaled < upper_rescaled)
#
#
# # set x end points by upper and lower limits
# # plot mean as point
#
# ggplot(data = bla_ok) +
#   geom_point(aes(x = mean_rescaled, y = y_value)) +
#   geom_segment(aes(x = lower_rescaled,
#                    xend = upper_rescaled,
#                    y = y_value,
#                    yend = y_value)) +
#   geom_point(data=bla_subset, aes(x = mean_rescaled, y = y_value)) +
#   geom_segment(data=bla_subset, aes(x = upper_rescaled,
#                    xend = seg1_x,
#                    y = y_value,
#                    yend = y_value)) +
#   geom_segment(data=bla_subset, aes(x = lower_rescaled,
#                                     xend = seg2_x,
#                                     y = y_value,
#                                     yend = y_value)) +
#   scale_x_continuous(limits = c(0,360),
#                      breaks = seq(0,359,30),
#                      labels = month.abb) +
#   scale_y_continuous(limits = c(-20,length(bla$mean_rescaled))) + #staring at '-10' so that points are squeezed in the center
#
#   annotate("rect", xmin = 330, xmax = 360, ymin = 0, ymax = length(bla$mean_rescaled), alpha = .2) + #Dec
#   annotate("rect", xmin = 0, xmax = 60, ymin = 0, ymax = length(bla$mean_rescaled), alpha = .2) + # jan - feb
#   annotate("rect", xmin = 150, xmax = 210, ymin = 0, ymax = length(bla$mean_rescaled), alpha = .2) + # jun-jul
#   coord_polar() +
#   labs(x="",
#        y="",
#        title = "Leaf turnover") +
#   theme(panel.grid.major.x = element_line(colour = "grey75",
#                                           size = 0.3),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.background = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.text.x = element_text(size = 11),
#         strip.text = element_text(face = "italic", size = 13),
#         strip.background = element_rect(fill="white"),
#         # legend.position = "right",
#         plot.margin=unit(c(1,0,0,0),"cm")
#   )
# #--------------------------------------------------------------------------
# #--------------------------------------------------------------------------
# #--------------------------------------------------------------------------
