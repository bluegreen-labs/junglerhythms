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
source("R/timeline_gap_fill.R")
source("R/consec_years.R")
source("R/fourier_function_species_level.R")
source("R/climate_coherence_function.R")
#----------------------------------------------------------------------

# Suppress summarise info
# because `summarise()` ungrouping output (override with `.groups` argument)
# comes up a lot in the consol, which is annoying
options(dplyr.summarise.inform = FALSE)

#----------------------------------------------------------------------
#----------------------------------------------------------------------
# specify phenophase
# out file will have this in filename
pheno_investigated = "leaf_turnover" #leaf_turnover  #leaf_dormancy
filename = paste("data/fourier_coherence_",pheno_investigated,".csv",sep="")
#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#--------   Phenology data --------------------------------------------
#----------------------------------------------------------------------
data <- readRDS("data/jungle_rhythms_data_manuscript_leaf_repro.rds")
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#--- get selected species and clean time series -----------------------
#----------------------------------------------------------------------
species_list <- unique(data$species_full)

# for selected species and phenophase: get extended timelines at ID level
timelines_id <- missing_year_gaps(data = data,
                                  species_name = species_list,
                                  pheno = pheno_investigated,
                                  gapfill_missingyears = 0)
# for fourier, no gaps in timelines allowed
# get longest consecutive timeline; at species level
timelines_sp_consec <- consecutive_timeline_sp(data = timelines_id,
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
                # sig_99,
                cycle_category)
                # spec_dom,
                # spec_zero,
                # spec_ratio)

fourier_sp$cycle_category <- ifelse(fourier_sp$sig_95 == FALSE, NA,
                               ifelse(fourier_sp$cycle_category == "No cyclicity", NA,
                               ifelse(fourier_sp$cycle_dom_months < 11, "sub-annual",
                                      ifelse(fourier_sp$cycle_dom_months > 13, "supra-annual","annual"))))
fourier_sp$cycle_dom_months <- ifelse(is.na(fourier_sp$cycle_category), NA, fourier_sp$cycle_dom_months)

# #---------------------------------------------------------------------
# #--- COHERENCE TEST with climate data at the species level -----------
# #---------------------------------------------------------------------
# # read climate data
# climate <- read.csv("data/yangambi_km5_monthly_kasongo.csv")
# climate$date <- paste(climate$year,climate$month,"15",sep = "-")
# climate$date <- as.Date(climate$date, "%Y-%m-%d")
# climate$date_monthly <- format(as.Date(climate$date), "%Y-%m")
#
# climate.avg <- read.csv("data/ClimData_monthly_avg.csv",header = TRUE,sep = ",")
# climate.avg = climate.avg[,(names(climate.avg) %in% c("Month","insol_JR","tmax_JR"))]
# colnames(climate.avg)[1] <- "month"
# climate <- merge(climate, climate.avg, by = "month", all.x = TRUE)
# climate <- climate %>%
#   dplyr::arrange(date)
#
# # coherence analysis
# coherence <- climate_coherence(data = timelines_sp_consec,
#                                climate = climate,
#                                species_name = species_list,
#                                pheno = pheno_investigated)
# # precip
# coherence$phase_diff_precip <- ifelse(coherence$signif == TRUE, coherence$phase_diff_precip, NA)
# coherence$coh_precip <- ifelse(coherence$signif == TRUE, coherence$coh_precip, NA)
# coherence$phase_diff_precip_neg <- ifelse(coherence$signif == TRUE, coherence$phase_diff_precip_neg, NA)
# coherence$coh_precip_neg <- ifelse(coherence$signif == TRUE, coherence$coh_precip_neg, NA)
# # sun
# coherence$phase_diff_sun <- ifelse(coherence$signif == TRUE, coherence$phase_diff_sun, NA)
# coherence$coh_sun <- ifelse(coherence$signif == TRUE, coherence$coh_sun, NA)
# coherence$phase_diff_sun_neg <- ifelse(coherence$signif == TRUE, coherence$phase_diff_sun_neg, NA)
# coherence$coh_sun_neg <- ifelse(coherence$signif == TRUE, coherence$coh_sun_neg, NA)
# # temp
# coherence$phase_diff_temp <- ifelse(coherence$signif == TRUE, coherence$phase_diff_temp, NA)
# coherence$coh_temp <- ifelse(coherence$signif == TRUE, coherence$coh_temp, NA)
# coherence$phase_diff_temp_neg <- ifelse(coherence$signif == TRUE, coherence$phase_diff_temp_neg, NA)
# coherence$coh_temp_neg <- ifelse(coherence$signif == TRUE, coherence$coh_temp_neg, NA)
#
#
# #---------------------------------------------------------------------
# #--- merge coherence and fourier ------------------------- -----------
# #---------------------------------------------------------------------
# fourier_stats <- merge(fourier_sp, coherence, by = c("species_full","phenophase"), all.x = TRUE)
fourier_stats <- fourier_sp
#--------------------------------------------------------------------
#--------------------------------------------------------------------
# write to file
#--------------------------------------------------------------------
#--------------------------------------------------------------------
write.table(fourier_stats, filename,
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")
