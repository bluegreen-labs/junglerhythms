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
data$species_full <- as.character(paste(data$genus_Meise, data$species_Meise))

# remove column id.x and rename id.y to id (--> in id.y, empty ids are renamed to EK1, EK2, etc...)
data$id <- as.character(data$id.y)
data = data[,!(names(data) %in% c("id.x","id.y"))]


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
#------------ Leaf turnover  -------------------------------------------
#-----------------------------------------------------------------------
data_LT <- data %>%
  filter(phenophase == "leaf_turnover") %>%
  # filter(layer == "canopy") %>%
  # filter(Ecology == "shade") %>%
  na.omit()

data_LT <- data_LT %>%
  group_by(species_full, week) %>%
  summarise(mean_week = mean(value, na.rm = TRUE),
            count_week = sum(value),
            total_week = length(value))

test <- data_LT %>%
  filter(species_full == "Combretum lokele")

# filter out species with too few total siteyears to have a meaningfull average year
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

#-----------------------------------------------------------------------
#------------ Leaf dormancy  -------------------------------------------
#-----------------------------------------------------------------------
data_LD <- data %>%
  filter(phenophase == "leaf_dormancy") %>%
  # filter(layer == "canopy") %>%
  # filter(Ecology == "shade") %>%
  na.omit()

data_LD <- data_LD %>%
  group_by(species_full, week) %>%
  summarise(mean_week = mean(value, na.rm = TRUE),
            count_week = sum(value),
            total_week = length(value))

# filter out species with too few total siteyears to have a meaningfull average year
data_LD <- data_LD %>%
  filter(total_week >= minimum_siteyears)

# data_LD <- data_LD %>%
#   filter(mean_week >= minimum_event_frequency)
data_LD$mean_week <- ifelse(data_LD$mean_week < minimum_event_frequency, 0, data_LD$mean_week)

# if(rescale_event){
#   data_LD <- data_LD %>%
#     group_by(species_full) %>%
#     mutate(mean_week = rescale(count_week, c(minimum_event_frequency,1)))
# }

data_LD_site <- inner_join(data_LD, census_site, by = c("species_full"))
data_LD_plot <- inner_join(data_LD, census_plot, by = c("species_full"))

final_LD_site <- data_LD_site %>%
  group_by(week) %>%
  summarise(ss = sum(mean_week*basal_area_site, na.rm = TRUE)/total_basal_area_site *100) #,
# ss2 = wtd.mean(mean_week, basal_area_site, na.rm = TRUE)*100,
# weighted_sd = sqrt(wtd.var(mean_week, basal_area_site, na.rm = TRUE))*100,
# weighted_CI = ss2 + 1.96*weighted_sd)

final_LD_plot <- data_LD_plot %>%
  group_by(Plot, week) %>%
  summarise(ss = sum(mean_week*basal_area_plot, na.rm = TRUE))
final_LD_plot <- inner_join(final_LD_plot, total_basal_area_plot, by = c("Plot"))
final_LD_plot <- final_LD_plot %>%
  mutate(ss = ss/total_basal_area_plot*100)
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
#------------ Leaf flushing  -------------------------------------------
#-----------------------------------------------------------------------
data_flush <- data %>%
  filter(phenophase == "leaf_dormancy") %>%
  filter(basal_area_site > 0)
data_flush$flushing_date <- paste(data_flush$year, data_flush$week, sep = "-")

data_flush_ones <- data_flush[which(data_flush$value != 0),]
flushing_timing <- data_flush_ones %>%
  group_by(species_full, id) %>%
  do(event_length(.))

flushing_timing$flushing_date <- paste(flushing_timing$year_end, flushing_timing$week_end +1, sep = "-")

flushing_timing <- flushing_timing[,(names(flushing_timing) %in% c("species_full",
                                                                   "id",
                                                                   "flushing_date"))]
flushing_timing$flushing_value <- 1

data_flush <- merge(data_flush, flushing_timing, by = c("species_full","id","flushing_date"), all.x = TRUE)
data_flush$flushing_value <- ifelse(is.na(data_flush$flushing_value),"0", data_flush$flushing_value)
data_flush$flushing_value <- ifelse(is.na(data_flush$value),NA, data_flush$flushing_value)
data_flush$flushing_value <- as.numeric(data_flush$flushing_value)

data_flush <- data_flush %>%
  group_by(species_full, week) %>%
  summarise(mean_week = mean(flushing_value, na.rm = TRUE),
            count_week = sum(flushing_value),
            total_week = length(flushing_value))

# filter out species with too few total siteyears to have a meaningfull average year
data_flush <- data_flush %>%
  filter(total_week >= minimum_siteyears)

data_flush$mean_week <- ifelse(data_flush$mean_week < minimum_event_frequency, 0, data_flush$mean_week)

data_flush_site <- inner_join(data_flush, census_site, by = c("species_full"))
data_flush_plot <- inner_join(data_flush, census_plot, by = c("species_full"))


final_flush_site <- data_flush_site %>%
  group_by(week) %>%
  summarise(ss = sum(mean_week*basal_area_site, na.rm = TRUE)/total_basal_area_site *100) #,

final_flush_plot <- data_flush_plot %>%
  group_by(Plot, week) %>%
  summarise(ss = sum(mean_week*basal_area_plot, na.rm = TRUE))
final_flush_plot <- inner_join(final_flush_plot, total_basal_area_plot, by = c("Plot"))
final_flush_plot <- final_flush_plot %>%
  mutate(ss = ss/total_basal_area_plot*100)
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# contribution of species to different peaks
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# ------  turnover -----------------------------------------------------
# # first peak
# data_peak1 <- data_LT_site %>%
#   filter(week %in% c(3:11))
# final_peak1 <- data_peak1 %>%
#   group_by(species_full) %>%
#   summarise(mean_freq = mean(mean_week),
#             basal_area = mean(basal_area_site),
#             ss = mean_freq * basal_area) %>%
#   arrange(desc(ss))
# total1 <- sum(final_peak1$ss)
# final_peak1$ss_perc <- final_peak1$ss/total1*100
# head(final_peak1)
# # second peak
# data_peak2 <- data_LT_site %>%
#   filter(week %in% c(17:20))
# final_peak2 <- data_peak2 %>%
#   group_by(species_full) %>%
#   summarise(mean_freq = mean(mean_week),
#             basal_area = mean(basal_area_site),
#             ss = mean_freq * basal_area) %>%
#   arrange(desc(ss))
# total2 <- sum(final_peak2$ss)
# final_peak2$ss_perc <- final_peak2$ss/total2*100
# head(final_peak2)
# # third peak, small incline in sept
# data_peak3 <- data_LT_site %>%
#   filter(week %in% c(34:36))
# final_peak3 <- data_peak3 %>%
#   group_by(species_full) %>%
#   summarise(mean_freq = mean(mean_week),
#             basal_area = mean(basal_area_site),
#             ss = mean_freq * basal_area) %>%
#   arrange(desc(ss))
# total3 <- sum(final_peak3$ss)
# final_peak3$ss_perc <- final_peak3$ss/total3*100
# head(final_peak3)
# # ------  dormancy -----------------------------------------------------
# # first peak
# data_LD_peak1 <- data_LD_site %>%
#   filter(week %in% c(5:7)) #%>% # c(1:8,41:48)
# final_LD_peak1 <- data_LD_peak1 %>%
#   group_by(species_full) %>%
#   summarise(mean_freq = mean(mean_week),
#             basal_area = mean(basal_area_site),
#             ss = mean_freq * basal_area) %>%
#   arrange(desc(ss))
# total1 <- sum(final_LD_peak1$ss)
# final_LD_peak1$ss_perc <- final_LD_peak1$ss/total1*100
# head(final_LD_peak1)
# # second peak
# data_LD_peak2 <- data_LD_site %>%
#   filter(week %in% c(9)) #%>% # c(1:8,41:48)
# final_LD_peak2 <- data_LD_peak2 %>%
#   group_by(species_full) %>%
#   summarise(mean_freq = mean(mean_week),
#             basal_area = mean(basal_area_site),
#             ss = mean_freq * basal_area) %>%
#   arrange(desc(ss))
# total2 <- sum(final_LD_peak2$ss)
# final_LD_peak2$ss_perc <- final_LD_peak2$ss/total2*100
# head(final_LD_peak2)
# # third peak, small incline in sept
# data_LD_peak3 <- data_LD_site %>%
#   filter(week %in% c(21:22)) #%>% # c(1:8,41:48)
# final_LD_peak3 <- data_LD_peak3 %>%
#   group_by(species_full) %>%
#   summarise(mean_freq = mean(mean_week),
#             basal_area = mean(basal_area_site),
#             ss = mean_freq * basal_area) %>%
#   arrange(desc(ss))
# total3 <- sum(final_LD_peak3$ss)
# final_LD_peak3$ss_perc <- final_LD_peak3$ss/total3*100
# head(final_LD_peak3)
# # forth peak, small incline in sept
# data_LD_peak4 <- data_LD_site %>%
#   filter(week %in% c(35:38)) #%>% # c(1:8,41:48)
# final_LD_peak4 <- data_LD_peak4 %>%
#   group_by(species_full) %>%
#   summarise(mean_freq = mean(mean_week),
#             basal_area = mean(basal_area_site),
#             ss = mean_freq * basal_area) %>%
#   arrange(desc(ss))
# total4 <- sum(final_LD_peak4$ss)
# final_LD_peak4$ss_perc <- final_LD_peak4$ss/total4*100
# head(final_LD_peak4)
# # fifth peak, small incline in sept
# data_LD_peak5 <- data_LD_site %>%
#   filter(week %in% c(45:47)) #%>% # c(1:8,41:48)
# final_LD_peak5 <- data_LD_peak5 %>%
#   group_by(species_full) %>%
#   summarise(mean_freq = mean(mean_week),
#             basal_area = mean(basal_area_site),
#             ss = mean_freq * basal_area) %>%
#   arrange(desc(ss))
# total5 <- sum(final_LD_peak5$ss)
# final_LD_peak5$ss_perc <- final_LD_peak5$ss/total5*100
# head(final_LD_peak5)
#-----------------------------------------------------------------------
# sp_driven <- data_LT_site %>%
#   filter(mean_week > 0)
# sp_driven$wtd <- sp_driven$mean_week * sp_driven$basal_area_site
# sp_driven$month <- ifelse(sp_driven$week %in% c(1:4),"Jan",
#                           ifelse(sp_driven$week %in% c(5:8),"Feb",
#                           ifelse(sp_driven$week %in% c(9:12),"Mar",
#                           ifelse(sp_driven$week %in% c(11:16),"Apr",
#                           ifelse(sp_driven$week %in% c(17:20),"May",
#                           ifelse(sp_driven$week %in% c(21:24),"Jun",
#                           ifelse(sp_driven$week %in% c(25:28),"Jul",
#                           ifelse(sp_driven$week %in% c(29:32),"Aug",
#                           ifelse(sp_driven$week %in% c(33:36),"Sep",
#                           ifelse(sp_driven$week %in% c(37:40),"Oct",
#                           ifelse(sp_driven$week %in% c(41:44),"Nov",
#                           ifelse(sp_driven$week %in% c(45:48),"Dec",
#                           NA))))))))))))
# # write to file
# write.table(sp_driven, "data/species_contribution_standlevel_turnover.csv",
#             quote = FALSE,
#             col.names = TRUE,
#             row.names = FALSE,
#             sep = ",")
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#------------ plots -------  -------------------------------------------
#-----------------------------------------------------------------------
combined <- final_LT_site
colnames(combined)[2] <- "ss_turn"
combined <- merge(combined, final_LD_site, by = c("week"), all.x = TRUE)
combined$ss_senescence <- combined$ss_turn + combined$ss

colnames(final_flush_site)[2] <- "ss_flush"
combined <- merge(combined, final_flush_site, by = c("week"), all.x = TRUE)


p_combined <- ggplot() +
  annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 2.2, alpha = .1) + # jan - febr
  annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 2.2, alpha = .1) + # jun - jul
  annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 2.2, alpha = .1) + # dec
  # senescence
  geom_smooth(data = combined,
              aes(week, ss_senescence, colour = "line4"), span = 0.2, se = FALSE, size = 1.2, linetype = "twodash",
              show.legend = TRUE) +
  # turnover
  geom_smooth(data = final_LT_site,
              aes(week, ss, colour = "line2"), span = 0.2, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(data = final_LT_site,
             aes(week, ss),
             col="#d8b365",
             shape = 1,
             stroke = 1.3
             ) +
  # dormancy
  geom_smooth(data = final_LD_site,
              aes(week, ss, colour = "line1"), span = 0.2, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(data = final_LD_site,
             aes(week, ss),
             col="#8c510a",
             size=2) +
  # flushing
  geom_smooth(data = combined,
              aes(week, ss_flush, colour = "line3"), span = 0.2, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(data = combined,
             aes(week, ss_flush),
             col="#018571",
             size=2) +

  scale_colour_manual(values = c("line1" = "#8c510a", "line2" = "#d8b365", "line3" = "#018571", "line4" = "grey30"),
                    labels = c(" dormacy   "," turnover   ", " flushing   "," senescence   ")) +
  scale_x_continuous(limits = c(1,49),
                     breaks = seq(1,48,4),
                     labels = month.abb) +
  labs(y = "freq. canopy phenophase",
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
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        plot.margin = unit(c(0,0,0,0.5),"cm")
  )
# p_combined

p_turnover <- ggplot() +
  geom_point(data = final_LT_plot,
             aes(week, ss, shape = Plot),
             col="grey40") +
  geom_smooth(data = final_LT_site,
              aes(week, ss), span = 0.2, se = FALSE, col = "black", size = 1.2) +
  # geom_line(data = final_LT_site,
  #           aes(week, ss),
  #           col="black",
  #           size=1.2) +
  geom_point(data = final_LT_site,
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
  labs(y = "freq. canopy turnover",
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
# print(p_turnover)


p_dormancy <- ggplot() +
  geom_point(data = final_LD_plot,
             aes(week, ss, shape = Plot),
             col="grey40") +
  geom_smooth(data = final_LD_site,
              aes(week, ss), span = 0.2, se = FALSE, col = "black", size = 1.2) +
  # geom_line(data = final_LD_site,
  #           aes(week, ss),
  #           col="black",
  #           size=1.2) +
  geom_point(data = final_LD_site,
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
  labs(y = "freq. canopy dormancy",
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
# print(p_dormancy)


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
p_modis <- ggplot_gtable(ggplot_build(p_modis))
p_turnover <- ggplot_gtable(ggplot_build(p_turnover))
p_precip <- ggplot_gtable(ggplot_build(p_precip))
p_sun <- ggplot_gtable(ggplot_build(p_sun))
p_tmax <- ggplot_gtable(ggplot_build(p_tmax))
p_combined <- ggplot_gtable(ggplot_build(p_combined))

p_modis$widths <-p_turnover$widths
p_precip$widths <-p_turnover$widths
p_sun$widths <-p_turnover$widths
p_tmax$widths <-p_turnover$widths
p_combined$widths <-p_turnover$widths

# p_all <- grid.arrange(p_combined, p_tmax, p_sun, p_precip, heights = c(5,1,1,4))
p_all <- grid.arrange(p_combined, p_modis, p_tmax, p_sun, p_precip, heights = c(5,4,1,1,4))

# supplementary information
p_all_phases <- grid.arrange(p_turnover, p_dormancy, heights = c(1,1.2))
p_all_modis <- grid.arrange(p_modis, heights = c(1))


# p_all <- grid.arrange(p_sun, p_precip, p_modis, p_dormancy, p_turnover, heights = c(1,2,3,3,3.4)) #


pdf("~/Desktop/figure4_standlevel.pdf",4,8.5) # 5,10)
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

combined$Month <- ifelse(combined$week %in% c(1,2,3,4),1,
                              ifelse(combined$week %in% c(5,6,7,8),2,
                                     ifelse(combined$week %in% c(9,10,11,12),3,
                                            ifelse(combined$week %in% c(13,14,15,16),4,
                                                   ifelse(combined$week %in% c(17,18,19,20),5,
                                                          ifelse(combined$week %in% c(21,22,23,24),6,
                                                                 ifelse(combined$week %in% c(25,26,27,28),7,
                                                                        ifelse(combined$week %in% c(29,30,31,32),8,
                                                                               ifelse(combined$week %in% c(33,34,35,36),9,
                                                                                      ifelse(combined$week %in% c(37,38,39,40),10,
                                                                                             ifelse(combined$week %in% c(41,42,43,44),11,
                                                                                                    ifelse(combined$week %in% c(45,46,47,48),12,
                                                                                                           NA))))))))))))

combined_month <- combined %>%
  group_by(Month) %>%
  summarise(dormancy_stand = mean(ss),
            turnover_stand = mean(ss_turn),
            senescence_stand = mean(ss_senescence),
            flushing_stand = mean(ss_flush))

climate.corr <- merge(climate, combined_month, by = c("Month"), all.x = TRUE)



# turnover
cor.test(climate.corr$insol_JR, climate.corr$turnover_stand, method = 'pearson')
cor.test(climate.corr$prec_JR, climate.corr$turnover_stand, method = 'pearson')
cor.test(climate.corr$tmax_JR, climate.corr$turnover_stand, method = 'pearson')


# dormancy
cor.test(climate.corr$insol_JR, climate.corr$dormancy_stand, method = 'pearson')
cor.test(climate.corr$prec_JR, climate.corr$dormancy_stand, method = 'pearson')
cor.test(climate.corr$tmax_JR, climate.corr$dormancy_stand, method = 'pearson')

# flushin
cor.test(climate.corr$insol_JR, climate.corr$flushing_stand, method = 'pearson')
cor.test(climate.corr$prec_JR, climate.corr$flushing_stand, method = 'pearson')
cor.test(climate.corr$tmax_JR, climate.corr$flushing_stand, method = 'pearson')

# senescence
cor.test(climate.corr$insol_JR, climate.corr$senescence_stand, method = 'pearson')
cor.test(climate.corr$prec_JR, climate.corr$senescence_stand, method = 'pearson')
cor.test(climate.corr$tmax_JR, climate.corr$senescence_stand, method = 'pearson')

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
cor.test(climate.corr$EVIm, climate.corr$tmax_JR, method = 'pearson')

cor.test(climate.corr$EVIm, climate.corr$dormancy_stand, method = 'pearson')
cor.test(climate.corr$EVIm, climate.corr$turnover_stand, method = 'pearson')
cor.test(climate.corr$EVIm, climate.corr$senescence_stand, method = 'pearson')
cor.test(climate.corr$EVIm, climate.corr$flushing_stand, method = 'pearson')
