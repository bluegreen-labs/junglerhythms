#----- reset your R session. ------------------------------------------
rm(list=ls())
#----- load required libraries ----------------------------------------
library(tidyverse)
library(plyr)
library(ggplot2)
library(gtools)
library(grid)
library(gridExtra)
#----- source required functions --------------------------------------
source("R/consec_years.R")
source("R/fourier_function_species_level.R")
#----------------------------------------------------------------------

# Suppress summarise info
# because `summarise()` ungrouping output (override with `.groups` argument)
# comes up a lot in the consol, which is annoying
options(dplyr.summarise.inform = FALSE)

#----------------------------------------------------------------------
#----------------------------------------------------------------------
# specify phenophase
# output file will have this in filename
pheno_investigated = "leaf_dormancy" #leaf_turnover  #leaf_dormancy
filename = paste("data/fourier_",pheno_investigated,".csv",sep="")
#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#--------   Phenology data --------------------------------------------
#----------------------------------------------------------------------
data <- readRDS("data/jungle_rhythms_data_manuscript_leaf_repro.rds")
species_list <- unique(data$species_full)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#--- for fourier, no gaps in timelines allowed ------------------------
#--- get longest consecutive timeline; at species level ---------------
#----------------------------------------------------------------------
timelines_sp_consec <- consecutive_timeline_sp(data = data,
                                               species_name = species_list,
                                               pheno = pheno_investigated)

#---------------------------------------------------------------------
#--- fourier at the species level ------------------------------------
#---------------------------------------------------------------------
fourier_sp <- peak_detection(data = timelines_sp_consec,
                            species_name = species_list,
                            pheno = pheno_investigated,
                            perio_plot = FALSE)

fourier_sp <- fourier_sp %>%
  mutate(cycle_dom_months = cycle_dom/48*12) %>%
  dplyr::select(species_full,
                phenophase,
                cycle_dom_months,
                sig_95,
                cycle_category)

fourier_sp$cycle_category <- ifelse(fourier_sp$sig_95 == FALSE, NA,
                               ifelse(fourier_sp$cycle_category == "No cyclicity", NA,
                               ifelse(fourier_sp$cycle_dom_months < 11, "sub-annual",
                                      ifelse(fourier_sp$cycle_dom_months > 13, "supra-annual","annual"))))
fourier_sp$cycle_dom_months <- ifelse(is.na(fourier_sp$cycle_category), NA, fourier_sp$cycle_dom_months)

#--------------------------------------------------------------------
#--------------------------------------------------------------------
# write to file
#--------------------------------------------------------------------
#--------------------------------------------------------------------
write.table(fourier_sp, filename,
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")
