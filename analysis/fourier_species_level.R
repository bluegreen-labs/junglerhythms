#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
#----- load required libraries -------------------------------------------------#
# library(TSA)
library(tidyverse)
library(plyr)
library(ggplot2)
library(gtools)
library(grid)
library(gridExtra)
# library(circular)
# library(fields)
# library(DescTools)
# library(truncnorm)
#----- source required functions -----------------------------------------------#
source("R/timeline_gap_fill.R")
source("R/consec_years.R")
source("R/fourier_function_species_level.R")
source("R/climate_coherence_function.R")
#-------------------------------------------------------------------------------#

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
rm(df, metadata, empty_years)
#----------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------
#--- get selected species and clean time series -----------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
overview <- read.csv("data/species_meta_data_phase2.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)
species_list <- overview$species_full

# for selected species and phenophase: get extended timelines at ID level with 2 year-gaps filled with zero
timelines_id <- two_year_gaps(data = data,
                            species_name = species_list,
                            pheno = pheno_investigated)
# for fourier, no gaps in timelines allowed
# get longest consecutive timeline; at species level
timelines_sp_consec <- consecutive_timeline_sp(data = timelines_id,
                             species_name = species_list,
                             pheno = pheno_investigated)

#----------------------------------------------------------------------------------------------------------------------
#--- fourier at the species level -------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
fourier_sp <- peak_detection(data = timelines_sp_consec,
                            species_name = species_list,
                            pheno = pheno_investigated,
                            perio_plot = FALSE)

fourier_sp <- fourier_sp %>%
  # filter(sig_cycle == "TRUE") %>%
  dplyr::select(species_full,
                phenophase,
                cycle_dom,
                sig_95,
                sig_99,
                cycle_category,
                spec_dom,
                spec_zero,
                spec_ratio) %>%
  mutate(cycle_dom_months = cycle_dom/48*12)
fourier_sp$cyclicity <- ifelse(fourier_sp$sig_95 == FALSE, NA,
                               ifelse(fourier_sp$cycle_category == "No cyclicity", NA,
                               ifelse(fourier_sp$cycle_dom_months < 11, "sub-annual",
                                      ifelse(fourier_sp$cycle_dom_months > 13, "supra-annual","annual"))))
fourier_sp$cycle_dom_months <- ifelse(is.na(fourier_sp$cyclicity), NA, fourier_sp$cycle_dom_months)


#----------------------------------------------------------------------------------------------------------------------
#--- COHERENCE TEST with climate data at the species level ------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
# read climate data
climate <- read.csv("data/yangambi_km5_monthly_kasongo.csv")
climate$date <- paste(climate$year,climate$month,"15",sep = "-")
climate$date <- as.Date(climate$date, "%Y-%m-%d")
climate$date_monthly <- format(as.Date(climate$date), "%Y-%m")

climate.avg <- read.csv("data/ClimData_monthly_avg.csv",header = TRUE,sep = ",")
climate.avg = climate.avg[,(names(climate.avg) %in% c("Month","insol_JR","tmax_JR"))]
colnames(climate.avg)[1] <- "month"
climate <- merge(climate, climate.avg, by = "month", all.x = TRUE)
climate <- climate %>%
  dplyr::arrange(date)

# coherence analysis
coherence <- climate_coherence(data = timelines_sp_consec,
                               climate = climate,
                               species_name = species_list,
                               pheno = pheno_investigated)
# precip
coherence$phase_diff_precip <- ifelse(coherence$signif == TRUE, coherence$phase_diff_precip, NA)
coherence$coh_precip <- ifelse(coherence$signif == TRUE, coherence$coh_precip, NA)
coherence$phase_diff_precip_neg <- ifelse(coherence$signif == TRUE, coherence$phase_diff_precip_neg, NA)
coherence$coh_precip_neg <- ifelse(coherence$signif == TRUE, coherence$coh_precip_neg, NA)
# sun
coherence$phase_diff_sun <- ifelse(coherence$signif == TRUE, coherence$phase_diff_sun, NA)
coherence$coh_sun <- ifelse(coherence$signif == TRUE, coherence$coh_sun, NA)
coherence$phase_diff_sun_neg <- ifelse(coherence$signif == TRUE, coherence$phase_diff_sun_neg, NA)
coherence$coh_sun_neg <- ifelse(coherence$signif == TRUE, coherence$coh_sun_neg, NA)
# temp
coherence$phase_diff_temp <- ifelse(coherence$signif == TRUE, coherence$phase_diff_temp, NA)
coherence$coh_temp <- ifelse(coherence$signif == TRUE, coherence$coh_temp, NA)
coherence$phase_diff_temp_neg <- ifelse(coherence$signif == TRUE, coherence$phase_diff_temp_neg, NA)
coherence$coh_temp_neg <- ifelse(coherence$signif == TRUE, coherence$coh_temp_neg, NA)


#----------------------------------------------------------------------------------------------------------------------
#--- merge coherence and fourier ------------------------- ------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
fourier_stats <- merge(fourier_sp, coherence, by = c("species_full","phenophase"), all.x = TRUE)

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


# #----- reset your R session. ---------------------------------------------------#
# rm(list=ls())
#
# #----------------------------------------------------------------------------------------------------------------------
# #----------------------------------------------------------------------------------------------------------------------
# four <- read.csv("data/fourier_coherence_leaf_turnover.csv",
#                  header = TRUE,
#                  sep = ",",
#                  stringsAsFactors = FALSE)
# corr <- read.csv("data/timeseries_correlations_phase2.csv",
#                  header = TRUE,
#                  sep = ",",
#                  stringsAsFactors = FALSE)
# df <- merge(four, corr, by = "species_full", all.x = TRUE)
#
# df <- df %>%
#   select(species_full,
#          cycle_dom_months,
#          phase_diff_precip,
#          phase_diff_precip_neg,
#          coh_precip,
#          corr_leaf_turnover_precip,
#          corr_leaf_turnover_precip_timing)
# df$cycle_dom_months <- round(df$cycle_dom_months, digits = 2)
# df$phase_diff_precip <- round(df$phase_diff_precip, digits = 2)
# df$phase_diff_precip_neg <- round(df$phase_diff_precip_neg, digits = 2)
# df$coh_precip <- round(df$coh_precip, digits = 2)
# df$corr_leaf_turnover_precip <- round(df$corr_leaf_turnover_precip, digits = 2)
#
#
#
# write.table(df, "data/test.csv",
#             quote = FALSE,
#             col.names = TRUE,
#             row.names = FALSE,
#             sep = ",")
