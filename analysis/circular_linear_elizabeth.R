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
#--------   Phenology data --------------------------------------------
#----------------------------------------------------------------------
data <- readRDS("data/jungle_rhythms_data_cleaned.rds")
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
timelines_id_dorm <- missing_year_gaps(data = data,
                                       species_name = species_list,
                                       pheno = "leaf_dormancy",
                                       gapfill_missingyears = 0)
timelines_id_turn <- missing_year_gaps(data = data,
                                       species_name = species_list,
                                       pheno = "leaf_turnover",
                                       gapfill_missingyears = 0)
timelines_id <- rbind(timelines_id_dorm, timelines_id_turn)
#----------------------------------------------------------------------
rm(timelines_id_dorm, timelines_id_turn)
#----------------------------------------------------------------------

# #----------------------------------------------------------------------
# # specieslist based on selection made in advance
# sp_evergreen <- c("Scorodophloeus zenkeri",
#                   "Strombosia pustulata",
#                   "Synsepalum subcordatum",
#                   "Dacryodes osika",
#                   "Quassia silvestris")
# sp_deciduous <- c("Petersianthus macrocarpus",
#                   "Vitex ferruginea",
#                   "Erythrophleum suaveolens",
#                   "Pterocarpus soyauxii",
#                   "Pericopsis elata")
#
# query_evergreen <- paste(sp_evergreen, collapse = "|")
# query_deciduous <- paste(sp_deciduous, collapse = "|")
#
# p_deciduous <- circular_linear_plot(timelines_id,
#                                     species_name = query_deciduous,
#                                     # leg_pos = c(1,0.1),
#                                     leg_gradient = c(0,0.3,1),
#                                     title_name = "(b) Deciduous")
# # pdf("~/Desktop/figure1b_deciduous.pdf",8.85,12.1)
# # plot(p_deciduous)
# # dev.off()
#
#
# p_evergreen <- circular_linear_plot(timelines_id,
#                                     species_name = query_evergreen,
#                                     title_name = "(a) Evergreen")
# pdf("~/Desktop/figure1a_evergreen.pdf",8.85,12.1)
# plot(p_evergreen)
# dev.off()



sp_all <- c("Dacryodes osika",
            "Quassia silvestris",
            "Scorodophloeus zenkeri",
            "Strombosia pustulata",
            "Synsepalum subcordatum",

            "Erythrophleum suaveolens",
            "Pericopsis elata",
            "Petersianthus macrocarpus",
            "Pterocarpus soyauxii",
            "Vitex ferruginea")

dummy <- c("a","b","c","d","e",
           "f","g","h","i","j")
sp_title <- as.data.frame(cbind(sp_all,  dummy))
sp_title <- sp_title %>%
  dplyr::rename(species_full = sp_all)
# facet labels
sp_title$dummy <- factor(sp_title$dummy,
                         levels = c("a","b","c","d","e",
                                    "f","g","h","i","j"),
                         labels = c("D. osika",
                                    "Q. silvestris",
                                    "S. zenkeri",
                                    "S. pustulata",
                                    "S. subcordatum",

                                    "E. suaveolens",
                                    "P. elata",
                                    "P. macrocarpus",
                                    "P. soyauxii",
                                    "V. ferruginea"))

timelines_id <- merge(timelines_id, sp_title, by = "species_full", all.x = TRUE)
sp_all <- paste(sp_all, collapse = "|")
p_all <- circular_plot(timelines_id,
                       species_name = sp_all,
                       title_name = "(a) evergreen                  (b) deciduous")




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


