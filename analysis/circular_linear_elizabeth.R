#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
# graphics.off()
#----- load required libraries -------------------------------------------------#
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(gridExtra)
#----- source required functions -----------------------------------------------#
# source("analysis/circular_plot_leaf_turnover_elizabeth.r")
# source("analysis/linear_plot_senescence_elizabeth.R")
source("analysis/circular_linear_function_elizabeth.R")
#-------------------------------------------------------------------------------#

#----------------------------------------------------------------------
#--------   Phenology data - species correction Meise   ---------------
#----------------------------------------------------------------------
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
# df <- df[which(df$value != 0),]
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
#----------------------------------------------------------------------

# specieslist based on selection made in advance
sp_evergreen <- c("Scorodophloeus zenkeri",
                  "Strombosia pustulata",
                  "Synsepalum subcordatum",
                  "Dacryodes osika",
                  "Quassia silvestris")
sp_deciduous <- c("Petersianthus macrocarpus",
                  "Irvingia grandifolia",
                  "Erythrophleum suaveolens",
                  "Antrocaryon nannanii",
                  "Pericopsis elata")

query_evergreen <- paste(sp_evergreen,collapse = "|")
query_deciduous <- paste(sp_deciduous,collapse = "|")

p_deciduous <- circular_linear_plot(data,
                                    species_name = query_deciduous,
                                    # leg_pos = c(1,0.1),
                                    leg_gradient = c(0,0.3,1),
                                    title_name = "(b) Deciduous")
pdf("~/Desktop/figure1b_deciduous.pdf",8.85,12.1)
plot(p_deciduous)
dev.off()

p_evergreen <- circular_linear_plot(data,
                                    species_name = query_evergreen,
                                    title_name = "(a) Evergreen")
pdf("~/Desktop/figure1a_evergreen.pdf",8.85,12.1)
plot(p_evergreen)
dev.off()

#---------------------------------------------------------------------
#---------------------------------------------------------------------
# getting the max frequency in the color scales of figure 1 evergreen v deciduous
#---------------------------------------------------------------------
#---------------------------------------------------------------------
# data_subset_circ <- data %>%
#   filter(species_full %in% "Irvingia grandifolia")
#
# # group by week and take the mean value
# data_subset_circ <- data_subset_circ %>%
#   group_by(species_full, week, phenophase) %>%
#   dplyr::summarise(mean_value = mean(value, na.rm=TRUE))#,
#
# data_subset_circ <- data_subset_circ %>%
#   filter(phenophase %in% "leaf_turnover")
#
# max(data_subset_circ$mean_value)


