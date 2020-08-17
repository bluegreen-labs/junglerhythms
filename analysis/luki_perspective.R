#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
# graphics.off()
#----- load required libraries -------------------------------------------------#
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(gridExtra)
#----- source required functions -----------------------------------------------#
source("R/luki_circular_linear.R")
#-------------------------------------------------------------------------------#

# #----------------------------------------------------------------------
# #--------   Phenology data - Luki -------------------------------------
# #----------------------------------------------------------------------
# df <- readRDS("data/luki_weekly_annotations.rds")
# metadata <- read.csv("data/luki_species_corrections.csv",
#                      header = TRUE, sep = ",")
# df <- df %>%
#   rename(luki_original = species_full)
# df <- merge(df, metadata, by = c("luki_original"), all.x = TRUE)
#
# species_list <- as.data.frame(unique(df$luki_original))
# colnames(species_list)[1] <- "luki_original"
# species_list <- merge(species_list, metadata, by = c("luki_original"), all.x = TRUE)


#-----------------------------------------

car <- df %>%
  filter(species_full %in% "Entandrophragma utile")#"Carapa procera")

# group by week and take the mean value
car <- car %>%
  group_by(species_full, week, phenophase) %>%
  dplyr::summarise(percent_value = mean(value, na.rm=TRUE) *100)

# # add a date column
# car$date <- as.Date(paste(car$year,
#                           round((car$week*10.167)-9),sep="-"), "%Y-%j")

species_list <- unique(df$species_full)[1:5]

species_list <- paste(species_list,collapse = "|")

plot_1 <- luki_circular_linear_plot(df,
                                  species_name = species_list,
                                  leg_gradient = c(0,0.3,1),
                                  title_name = "luki")
