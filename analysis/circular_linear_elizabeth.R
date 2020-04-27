#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
# graphics.off()
#----- load required libraries -------------------------------------------------#
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(gridExtra)
#----- source required functions -----------------------------------------------#
source("R/circular_linear_function_elizabeth.R")
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
rm(df,metadata, empty_years)
#----------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------
#--- get selected species and clean time series -----------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
species_list <- c("Scorodophloeus zenkeri",
                  "Strombosia pustulata",
                  "Synsepalum subcordatum",
                  "Dacryodes osika",
                  "Quassia silvestris",
                  "Petersianthus macrocarpus",
                  "Vitex ferruginea",
                  "Erythrophleum suaveolens",
                  "Pterocarpus soyauxii",
                  "Pericopsis elata")

# for selected species and phenophase: get extended timelines at ID level with 2 year-gaps filled with zero
timelines_id_dorm <- two_year_gaps(data = data,
                              species_name = species_list,
                              pheno = "leaf_dormancy")
timelines_id_turn <- two_year_gaps(data = data,
                                   species_name = species_list,
                                   pheno = "leaf_turnover")
timelines_id <- rbind(timelines_id_dorm, timelines_id_turn)
#----------------------------------------------------------------------
rm(timelines_id_dorm, timelines_id_turn)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# specieslist based on selection made in advance
sp_evergreen <- c("Scorodophloeus zenkeri",
                  "Strombosia pustulata",
                  "Synsepalum subcordatum",
                  "Dacryodes osika",
                  "Quassia silvestris")
sp_deciduous <- c("Petersianthus macrocarpus",
                  "Vitex ferruginea",
                  "Erythrophleum suaveolens",
                  "Pterocarpus soyauxii",
                  "Pericopsis elata")

query_evergreen <- paste(sp_evergreen, collapse = "|")
query_deciduous <- paste(sp_deciduous, collapse = "|")

p_deciduous <- circular_linear_plot(timelines_id,
                                    species_name = query_deciduous,
                                    # leg_pos = c(1,0.1),
                                    leg_gradient = c(0,0.3,1),
                                    title_name = "(b) Deciduous")
pdf("~/Desktop/figure1b_deciduous.pdf",8.85,12.1)
plot(p_deciduous)
dev.off()


p_evergreen <- circular_linear_plot(timelines_id,
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


