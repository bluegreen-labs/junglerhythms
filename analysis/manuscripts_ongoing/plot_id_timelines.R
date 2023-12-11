rm(list=ls())
#----- load required libraries -------------------------------------------------#
library(tidyverse)
library(plyr)
library(ggplot2)
library(gtools)
library(grid)
library(gridExtra)
library(tseries)
#----- source required functions -----------------------------------------------#
source("R/timeline_gap_fill.R")
source("R/consec_years.R")
#-------------------------------------------------------------------------------#


#----------------------------------------------------------------------
#--------   Phenology data - species correction Meise   ---------------
#----------------------------------------------------------------------
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
df$join_id <- paste0("R",df$image,"-",df$image_row)
metadata <- read.csv("data/phenology_archives_species_long_format_20200324.csv",
                     header = TRUE, sep = ",")
metadata$join_id <- paste(metadata$image,metadata$row, sep = "-")

# test merge the two tables based upon the unique join ID
data <- merge(df, metadata, by = c("join_id"), all.x = TRUE)
data$species_full <- paste(data$genus_Meise, data$species_Meise)

# remove column id.x and rename id.y to id (--> in id.y, empty ids are renamed to EK1, EK2, etc...)
data$id <- as.character(data$id.y)
data = data[,!(names(data) %in% c("id.x","id.y"))]
# data$id <- as.character(data$id)
# data$species_full <- as.character(data$species_full)

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

#----------------------------------------------------------------------------------------------------------------------
#--- get selected species and clean time series -----------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
overview <- read.csv("data/species_meta_data_phase2.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)

species_list <- overview$species_full
#----------------------------------------------------------------------------------------------------------------------
overview_dorm <- read.csv("data/SI_table2_dormancy.csv",
                          header = TRUE,
                          sep = ",",
                          stringsAsFactors = FALSE)

overview_dorm$dorm_corrs <- ifelse(overview_dorm$corr.precip != "" |
                                     overview_dorm$corr.sun != "" |
                                     overview_dorm$corr.tmax != "", "yes", NA)
dorm_sp1 <- overview_dorm %>%
  filter(cyclicity_dormancy %in% c("annual"))
dorm_sp2 <- overview_dorm %>%
  filter(!cyclicity_dormancy == "annual") %>%
  filter(dorm_corrs == "yes")

dorm_rest <- overview_dorm %>%
  filter(!cyclicity_dormancy == "annual") %>%
  filter(is.na(dorm_corrs))

dorm_sp1 <- dorm_sp1$Species
dorm_sp2 <- dorm_sp2$Species
dorm_rest <- dorm_rest$Species

#--- leaf turnover ----------------------------------------------------------------------------------------------------
# for selected species and phenophase: get extended timelines at ID level with 2 year-gaps filled with zero
# timelines_id <- missing_year_gaps(data = data,
#                                        species_name = species_list,
#                                        pheno = "flowers",
#                                        gapfill_missingyears = 0)

timelines_id <- missing_year_gaps(data = data,
                                  species_name = dorm_sp1,
                                  pheno = "leaf_dormancy",
                                  gapfill_missingyears = 0)


#------------------------------------------------------------------------
# linear leaf turnover
#------------------------------------------------------------------------


# species_name <- c("Erythrophleum suaveolens", "Panda oleosa")
# species_name <- species_list
#
# for (j in 1:length(species_name)){
#
#   # subset on species
#   data_subset <- timelines_id %>%
#     filter(species_full == species_name[j])
#
#   plot_name <- paste("~/Desktop/ids/flowers/",species_name[j], "_flowers",".png",sep = "")
#
#
#
#   p_lin <- ggplot(data_subset, aes(x = date,
#                                    y = value)) +
#     geom_line() +
#     theme_minimal() +
#     labs(title = species_name[j],
#          y = "Obs",
#          x = "Year") +
#     scale_x_date(date_breaks = "1 years",
#                  date_labels = "%Y",
#                  limits = as.Date(c('1937-01-01','1956-12-31')),
#                  expand = c(0, 0)) +
#     scale_y_continuous(limits = c(0,1),
#                        breaks = c(0,0.5,1)) +
#     theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
#           panel.grid.minor.x = element_blank(),
#           panel.grid.minor.y = element_blank(),
#           panel.background = element_blank(),
#           plot.background = element_rect(fill = 'white', colour = 'white'),
#           plot.title = element_text(),
#           strip.text = element_text(hjust = 0, size = 13, face = "italic"),
#           axis.line.x = element_blank(),
#           axis.title.x = element_blank(),
#           axis.text.x = element_text(angle = 90, hjust = 1),
#           axis.ticks.x = element_blank(),
#           legend.position = "none",
#           plot.margin = unit(c(1,1,1,1),"cm")
#     ) +
#     facet_wrap( ~ id, ncol = 1)
#
#   png(plot_name, width = 925, height = 700)
#   plot(p_lin)
#   dev.off()
#
# }

species_name <- dorm_sp1

for (j in 1:length(species_name)){

  # subset on species
  data_subset <- timelines_id %>%
    filter(species_full == species_name[j])

  # plot_name <- paste("~/Desktop/ids/flowers/",species_name[j], "_flowers",".png",sep = "")



  p_lin <- ggplot(data_subset, aes(x = date,
                                   y = value)) +
    geom_line() +
    theme_minimal() +
    labs(title = species_name[j],
         y = "Obs",
         x = "Year") +
    scale_x_date(date_breaks = "1 years",
                 date_labels = "%Y",
                 limits = as.Date(c('1937-01-01','1956-12-31')),
                 expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,1),
                       breaks = c(0,0.5,1)) +
    theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = 'white', colour = 'white'),
          plot.title = element_text(),
          strip.text = element_text(hjust = 0, size = 13, face = "italic"),
          axis.line.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(1,1,1,1),"cm")
    ) +
    facet_wrap( ~ id, ncol = 1)

  # png(plot_name, width = 925, height = 700)
  # plot(p_lin)
  # dev.off()

  plot(p_lin)

}

data_subset <- data %>%
  filter(species_full %in% dorm_sp1) %>%
  filter(!id == "101")

timelines_id <- missing_year_gaps(data = data_subset,
                                  species_name = dorm_sp1,
                                  pheno = "leaf_dormancy",
                                  gapfill_missingyears = 0)

head(timelines_id)

timelines_id <- timelines_id %>%
  dplyr::filter(year > 1938 & year < 1953)

timeline_compiled <- timelines_id %>%
  group_by(species_full, date, phenophase) %>%
  dplyr::summarise(mean_value = mean(value, na.rm=TRUE))

p_lin <- ggplot(timeline_compiled, aes(x = date,
                                 y = mean_value)) +
  geom_line() +
  theme_minimal() +
  labs(y = "Obs",
       x = "Year") +
  scale_x_date(date_breaks = "1 years",
               date_labels = "%Y",
               limits = as.Date(c('1939-01-01','1952-12-31')),
               expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0,0.5,1)) +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(),
        strip.text = element_text(hjust = 0, size = 13, face = "italic"),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm")
  ) +
  facet_wrap( ~ species_full, ncol = 1)

write.table(timeline_compiled, "data/timeseries_dormancy_annual_sp.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")
