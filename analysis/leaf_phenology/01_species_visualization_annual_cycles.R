#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
graphics.off()
#----- load required libraries -------------------------------------------------#
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(gridExtra)
#----- source required functions -----------------------------------------------#
source("R/circular_linear_plots.R")
#-------------------------------------------------------------------------------#

#----------------------------------------------------------------------
#--------   Phenology data --------------------------------------------
#----------------------------------------------------------------------
data <- readRDS("data/jungle_rhythms_data_manuscript_leaf_repro.rds")
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#--- get selected species ---------------------------------------------
#----------------------------------------------------------------------
species_list <- c("Dacryodes osika", # evergreens
                  "Quassia silvestris",
                  "Scorodophloeus zenkeri",
                  "Strombosia pustulata",
                  "Synsepalum subcordatum",
                  "Erythrophleum suaveolens", # deciduous
                  "Pericopsis elata",
                  "Petersianthus macrocarpus",
                  "Pterocarpus soyauxii",
                  "Vitex ferruginea")

# dummies to change alphabetical species order for the figure
dummy <- c("a","b","c","d","e",
           "f","g","h","i","j")
sp_title <- as.data.frame(cbind(species_list,  dummy))
# link facet labels to dymmy
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
sp_title <- sp_title %>%
  dplyr::rename(species_full = species_list)
data <- merge(data, sp_title, by = "species_full", all.x = TRUE)

#----------------------------------------------------------------------
# render figure for selected species
#----------------------------------------------------------------------
species_list <- paste(species_list, collapse = "|")
p_all <- circular_plot(data = data,
                       species_name = species_list,
                       title_name = "(a) evergreen                  (b) deciduous")

ggsave("manuscript/leaf_phenology/figures/fig1.png", p_all, device = "png", width = 6, height = 12)



