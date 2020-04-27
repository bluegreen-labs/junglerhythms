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
source("R/timeline_gap_fill.R")
source("R/standlevel_phen.R")
#-------------------------------------------------------------------------------#


#-------------------------------------------------------------------------------#
# SETTINGS                                                                      #
# Change these setting according to the output you want.                        #
#-------------------------------------------------------------------------------#
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
minimum_siteyears = 10
minimum_event_frequency = 0 # range c(0,1)
#-------------------------------------------------------------------------------#


#----------------------------------------------------------------------
#--------   Phenology data - species correction Meise   ---------------
#----------------------------------------------------------------------
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
df$join_id <- paste0("R",df$image,"-",df$image_row)

metadata <- read.csv("data/phenology_archives_species_long_format_20200324.csv",
                     header = TRUE, sep = ",")
metadata$join_id <- paste(metadata$image,metadata$row, sep = "-")
# merge the two tables based upon the unique join ID
data <- merge(df, metadata, by = c("join_id"), all.x = TRUE)
data$species_full <- paste(data$genus_Meise, data$species_Meise)

# remove column id.x and rename id.y to id (--> in id.y, empty ids are renamed to EK1, EK2, etc...)
data$id <- as.character(data$id.y)
data = data[,!(names(data) %in% c("id.x","id.y"))]

# remove rows with NA's in year -> individuals with 'no_data' in the archive
data <- data[!(is.na(data$year)),]

# sum events for each id, each year, across phenophases
# years with zero observations across phenophases are possibly not observed
empty_years <- data %>%
  group_by(species_full,join_id,year) %>%
  dplyr::summarise(check_empty_years = sum(value))
data <- merge(data, empty_years, by = c("join_id","species_full","year"), all.x = TRUE)
data <- data %>%
  filter(check_empty_years > 0)
#----------------------------------------------------------------------
# only select parameters you need, more clear structure to work with
data <- data %>%
  select(species_full,
         id,
         phenophase,
         year,
         week,
         value)
#----------------------------------------------------------------------
rm(df,metadata, empty_years)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#-------- load climate data -------------------------------------------
#----------------------------------------------------------------------
climate <- read.csv("data/ClimData_monthly_avg.csv",
                    header = TRUE,
                    sep = ",")
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
  dplyr::rename("species_full" = Species)
census$species_full <- ifelse(census$species_full == "Unknown", NA, census$species_full)

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
#-- for this species list at ID level: --------------------------------
#-- get full range timelines with 2-year gap filled -------------------
#----------------------------------------------------------------------
overview_dorm <- read.csv("data/SI_table2_dormancy.csv",
                          header = TRUE,
                          sep = ",",
                          stringsAsFactors = FALSE)
dorm_sp1 <- overview_dorm %>%
  filter(cyclicity_dormancy %in% c("annual"))
dorm_rest <- overview_dorm %>%
  filter(!cyclicity_dormancy == "annual")

dorm_sp1 <- dorm_sp1$Species
dorm_rest <- dorm_rest$Species

overview_turn <- read.csv("data/SI_table2_turnover.csv",
                          header = TRUE,
                          sep = ",",
                          stringsAsFactors = FALSE)
turn_sp1 <- overview_turn %>%
  filter(cyclicity_turnover %in% c("annual"))
turn_rest <- overview_turn %>%
  filter(!cyclicity_turnover == "annual")

turn_sp1 <- turn_sp1$Species
turn_rest <- turn_rest$Species

all_species_list <- overview_dorm$Species


#----------------------------------------------------------------------
#----------------------------------------------------------------------
standlevel_ann <- standlevel_phen(data = data,
                                  species_list_dorm = dorm_sp1,
                                  species_list_turn = turn_sp1)
standlevel_rest <- standlevel_phen(data = data,
                                 species_list_dorm = dorm_rest,
                                 species_list_turn = turn_rest)

standlevel_full <- standlevel_phen(data = data,
                                   species_list_dorm = all_species_list,
                                   species_list_turn = all_species_list)
#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# left side panel
#----------------------------------------------------------------------
#----------------------------------------------------------------------
p_combined_all <- ggplot() +
  annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 2.2, alpha = .1) + # jan - febr
  annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 2.2, alpha = .1) + # jun - jul
  annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 2.2, alpha = .1) + # dec
  # senescence
  geom_smooth(data = standlevel_full,
              aes(week, ss_senescence, colour = "line4"), span = 0.2, se = FALSE, size = 1.2, linetype = "twodash",
              show.legend = TRUE) +
  # turnover
  geom_smooth(data = standlevel_full,
              aes(week, ss_turn, colour = "line2"), span = 0.2, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(data = standlevel_full,
             aes(week, ss_turn),
             col="#d8b365",
             shape = 1,
             stroke = 1.3) +
  # dormancy
  geom_smooth(data = standlevel_full,
              aes(week, ss, colour = "line1"), span = 0.2, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(data = standlevel_full,
             aes(week, ss),
             col="#8c510a",
             size=2) +
  # flushing
  geom_smooth(data = standlevel_full,
              aes(week, ss_flush, colour = "line3"), span = 0.2, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(data = standlevel_full,
             aes(week, ss_flush),
             col="#018571",
             size=2) +

  scale_colour_manual(values = c("line1" = "#8c510a", "line2" = "#d8b365", "line3" = "#018571", "line4" = "grey30"),
                      labels = c(" dormacy   "," turnover   ", " flushing   "," senescence   ")) +
  scale_x_continuous(limits = c(1,49),
                     breaks = seq(1,48,4),
                     labels = month.abb) +
  scale_y_continuous(limits = c(0,2.2)) +
  labs(y = "% of canopy in state of phenophase",
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
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        plot.margin = unit(c(0,0,0,0.5),"cm")
  )

p_precip <- ggplot(climate) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = 0, ymax = 330, alpha = .1) + # jan - febr
  annotate("rect", xmin = 5.5, xmax = 7.5, ymin = 0, ymax = 330, alpha = .1) + # jun - jul
  annotate("rect", xmin = 11.5, xmax = 12.5, ymin = 0, ymax = 330, alpha = .1) + # dec

  annotate("text", x = 1.5, y = 310, label = "LD", col = "grey50") +
  annotate("text", x = 4, y = 310, label = "SW", col = "grey50") +
  annotate("text", x = 6.5, y = 310, label = "SD", col = "grey50") +
  annotate("text", x = 9.5, y = 310, label = "LW", col = "grey50") +

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

#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# right side panel
#----------------------------------------------------------------------
#----------------------------------------------------------------------
p_combined_ann <- ggplot() +
  annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 1.6, alpha = .1) + # jan - febr
  annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 1.6, alpha = .1) + # jun - jul
  annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 1.6, alpha = .1) + # dec

  annotate("text", x = 35, y = 1.5, label = "species with annual cycle", col = "grey10") +

  # senescence
  geom_smooth(data = standlevel_ann,
              aes(week, ss_senescence, colour = "line4"), span = 0.2, se = FALSE, size = 1.2, linetype = "twodash",
              show.legend = TRUE) +
  # turnover
  geom_smooth(data = standlevel_ann,
              aes(week, ss_turn, colour = "line2"), span = 0.2, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(data = standlevel_ann,
             aes(week, ss_turn),
             col="#d8b365",
             shape = 1,
             stroke = 1.3) +
  # dormancy
  geom_smooth(data = standlevel_ann,
              aes(week, ss, colour = "line1"), span = 0.2, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(data = standlevel_ann,
             aes(week, ss),
             col="#8c510a",
             size=2) +
  # flushing
  geom_smooth(data = standlevel_ann,
              aes(week, ss_flush, colour = "line3"), span = 0.2, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(data = standlevel_ann,
             aes(week, ss_flush),
             col="#018571",
             size=2) +

  scale_colour_manual(values = c("line1" = "#8c510a", "line2" = "#d8b365", "line3" = "#018571", "line4" = "grey30"),
                      labels = c(" dormacy   "," turnover   ", " flushing   "," senescence   ")) +
  scale_x_continuous(limits = c(1,49),
                     breaks = seq(1,48,4),
                     labels = month.abb) +
  scale_y_continuous(limits = c(0,1.6),
                     breaks = c(0,0.5,1,1.5)) +
  labs(y = "",
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
        axis.title.y = element_text(vjust = 3),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        plot.margin = unit(c(0,0,0,0.5),"cm")
  )
#----------------------------------------------------------------------
p_combined_rest <- ggplot() +
  annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 1.6, alpha = .1) + # jan - febr
  annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 1.6, alpha = .1) + # jun - jul
  annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 1.6, alpha = .1) + # dec

  annotate("text", x = 35, y = 1.5, label = "species with no annual cycle", col = "grey10") +

  # senescence
  geom_smooth(data = standlevel_rest,
              aes(week, ss_senescence, colour = "line4"), span = 0.2, se = FALSE, size = 1.2, linetype = "twodash",
              show.legend = TRUE) +
  # turnover
  geom_smooth(data = standlevel_rest,
              aes(week, ss_turn, colour = "line2"), span = 0.2, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(data = standlevel_rest,
             aes(week, ss_turn),
             col="#d8b365",
             shape = 1,
             stroke = 1.3) +
  # dormancy
  geom_smooth(data = standlevel_rest,
              aes(week, ss, colour = "line1"), span = 0.2, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(data = standlevel_rest,
             aes(week, ss),
             col="#8c510a",
             size=2) +
  # flushing
  geom_smooth(data = standlevel_rest,
              aes(week, ss_flush, colour = "line3"), span = 0.2, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(data = standlevel_rest,
             aes(week, ss_flush),
             col="#018571",
             size=2) +

  scale_colour_manual(values = c("line1" = "#8c510a", "line2" = "#d8b365", "line3" = "#018571", "line4" = "grey30"),
                      labels = c(" dormacy   "," turnover   ", " flushing   "," senescence   ")) +
  scale_x_continuous(limits = c(1,49),
                     breaks = seq(1,48,4),
                     labels = month.abb) +
  scale_y_continuous(limits = c(0,1.6),
                     breaks = c(0,0.5,1,1.5)) +
  labs(y = "                                              % of canopy in state of phenophase",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 1.5, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 3),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        plot.margin = unit(c(0,0,0,0.5),"cm")
  )

p_species_ann <- ggplot(standlevel_ann, aes(week, ss_senescence)) +
  annotate("segment", x = 4, xend = 4, y = 0, yend = 1) +
  annotate("text", x = 4.5, y = 1.5, label = "species1", col = "grey10") +
  annotate("text", x = 4.5, y = 1.4, label = "species2", col = "grey10") +
  annotate("text", x = 4.5, y = 1.3, label = "species3", col = "grey10") +
  annotate("text", x = 4.5, y = 1.2, label = "species4", col = "grey10") +
  scale_x_continuous(limits = c(1,49)) +
  scale_y_continuous(limits = c(0,2)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

p_species_rest <- ggplot(standlevel_rest, aes(week, ss_senescence)) +
  annotate("segment", x = 9, xend = 9, y = 0, yend = 1) +
  annotate("text", x = 9.5, y = 1.5, label = "species1", col = "grey10") +
  annotate("text", x = 9.5, y = 1.4, label = "species2", col = "grey10") +
  annotate("text", x = 9.5, y = 1.3, label = "species3", col = "grey10") +
  annotate("text", x = 9.5, y = 1.2, label = "species4", col = "grey10") +


  annotate("segment", x = 17, xend = 17, y = 0, yend = 1) +
  annotate("text", x = 17.5, y = 1.5, label = "species1", col = "grey10") +
  annotate("text", x = 17.5, y = 1.4, label = "species2", col = "grey10") +
  annotate("text", x = 17.5, y = 1.3, label = "species3", col = "grey10") +
  annotate("text", x = 17.5, y = 1.2, label = "species4", col = "grey10") +

  annotate("segment", x = 34, xend = 34, y = 0, yend = 1) +
  annotate("text", x = 34.5, y = 1.5, label = "species1", col = "grey10") +
  annotate("text", x = 34.5, y = 1.4, label = "species2", col = "grey10") +
  annotate("text", x = 34.5, y = 1.3, label = "species3", col = "grey10") +
  annotate("text", x = 34.5, y = 1.2, label = "species4", col = "grey10") +

  scale_x_continuous(limits = c(1,49)) +
  scale_y_continuous(limits = c(0,2)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

p_empty <- ggplot() +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())



#-----------------------------------------------------------------------
# combine plots
#-----------------------------------------------------------------------
p_modis <- ggplot_gtable(ggplot_build(p_modis))
p_precip <- ggplot_gtable(ggplot_build(p_precip))
p_sun <- ggplot_gtable(ggplot_build(p_sun))
p_tmax <- ggplot_gtable(ggplot_build(p_tmax))
p_combined_all <- ggplot_gtable(ggplot_build(p_combined_all))
p_combined_ann <- ggplot_gtable(ggplot_build(p_combined_ann))
p_combined_rest <- ggplot_gtable(ggplot_build(p_combined_rest))
p_species_ann <- ggplot_gtable(ggplot_build(p_species_ann))
p_species_rest <- ggplot_gtable(ggplot_build(p_species_rest))


p_modis$widths <-p_combined_all$widths
p_precip$widths <-p_combined_all$widths
p_sun$widths <-p_combined_all$widths
p_tmax$widths <-p_combined_all$widths
p_combined_ann$widths <-p_combined_all$widths
p_combined_rest$widths <-p_combined_all$widths
p_species_ann$widths <-p_combined_all$widths
p_species_rest$widths <-p_combined_all$widths

p_empty <- ggplot() + theme_void()

p_all <- grid.arrange(p_combined_all, p_modis, p_tmax, p_sun, p_precip, heights = c(4,4,1,1,4))
# p_sub <- grid.arrange(p_combined_ann,p_species_ann,p_species_rest, p_combined_rest,heights = c(4,4,4,5))
p_sub <- grid.arrange(p_combined_ann, p_combined_rest, p_empty, heights = c(4,5,5))
p_overview <- grid.arrange(p_all, p_sub, ncol=2)


pdf("~/Desktop/test4.pdf",9,10)
plot(p_overview)
dev.off()

