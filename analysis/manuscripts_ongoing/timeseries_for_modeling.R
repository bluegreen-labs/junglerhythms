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
species_list <- c("Anthonotha macrophylla",
                  "Entandrophragma utile",
                  "Erythrophleum suaveolens",
                  "Fernandoa adolfi-friderici",
                  "Milicia excelsa",
                  "Monodora angolensis",
                  "Pericopsis elata",
                  "Ricinodendron heudelotii"  )

#--- timelines; no gap filling (zero missing years filled) ------------------------------------------------------------
timelines_id <- missing_year_gaps(data = data,
                                  species_name = species_list,
                                  pheno = "leaf_dormancy",
                                  gapfill_missingyears = 0)


# check what the timeseries of each individual of a species looks like
for (j in 1:length(species_list)){

  data_subset <- timelines_id %>%
    filter(species_full == species_list[j])


  p_lin <- ggplot(data_subset, aes(x = date,
                                   y = value)) +
    geom_line() +
    theme_minimal() +
    labs(title = species_list[j],
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

  plot(p_lin)

}

# Fernandoa adolfi-friderici, ID 101 is strange outlier. removed for this analysis
timelines_id <- timelines_id %>%
  filter(!id == "101")

# most timeseries few data after 1952 and before 1939
timelines_id <- timelines_id %>%
  dplyr::filter(year > 1938 & year < 1953)

# mean values for each date at species level
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
