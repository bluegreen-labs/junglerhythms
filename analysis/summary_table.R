#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
# graphics.off()
#----- load required packages --------------------------------------------------#
library(tidyverse)

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# add family name !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# read in summary data (file made in summary_statistics_elizabeth.r)
#-----------------------------------------------------------------------
overview <- read.csv("data/species_meta_data.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)
# remove some columns you don't want in the overview table
overview = overview[,!(names(overview) %in% c("start_year",
                                              "end_year",
                                              "site_years_with_leaf_dormancy",
                                              "site_years_with_leaf_turnover",
                                              "basal_area_site",
                                              "abundance",
                                              "dbh_max_mm_census",
                                              "mating_system",
                                              "diaspora",
                                              "height_m_literature"))]
# round up some numbers
overview$percentage_site_years_with_leaf_dormancy <- round(overview$ratio_site_years_with_leaf_dormancy, digits=3) *100
overview$percentage_site_years_with_leaf_turnover <- round(overview$ratio_site_years_with_leaf_turnover, digits=3) *100
overview = overview[,!(names(overview) %in% c("ratio_site_years_with_leaf_dormancy",
                                              "ratio_site_years_with_leaf_turnover"))]
overview$basal_area_percentage <- round(overview$basal_area_percentage, digits=2)

#----------------------------------
# combine means with SD and (n = )
# and then remove the original columns to get an overview of what's left
overview$duration_dormancy_weeks <- paste (round(overview$mean_duration_leaf_dormancy_weeks, digits=1),
                                           "±",
                                           round(overview$sd_duration_leaf_dormancy_weeks, digits=1),
                                           "(n =",
                                           overview$total_nr_events_leaf_dormancy,
                                           ")")
overview = overview[,!(names(overview) %in% c("mean_duration_leaf_dormancy_weeks",
                                              "sd_duration_leaf_dormancy_weeks",
                                              "total_nr_events_leaf_dormancy"))]
#----------------------------------
overview$duration_turnover_weeks <- paste (round(overview$mean_duration_leaf_turnover_weeks, digits=1),
                                           "±",
                                           round(overview$sd_duration_leaf_turnover_weeks, digits=1),
                                           "(n =",
                                           overview$total_nr_events_leaf_turnover,
                                           ")")
overview = overview[,!(names(overview) %in% c("mean_duration_leaf_turnover_weeks",
                                              "sd_duration_leaf_turnover_weeks",
                                              "total_nr_events_leaf_turnover"))]
#----------------------------------
#----------------------------------
# overview$median_onset_dormancy_weeks <- paste (round(overview$median_intrasp_onset_leaf_dormancy_weeks, digits=1),
#                                            "(",
#                                            round(overview$mean_intrasp_onset_leaf_dormancy_weeks, digits=1),
#                                            "±",
#                                            round(overview$sd_intrasp_onset_leaf_dormancy_weeks, digits=1),
#                                            ")")
overview$median_onset_dormancy_weeks <- round(overview$median_intrasp_onset_leaf_dormancy_weeks, digits=1)
overview$intra_species_synchr_dormancy_weeks <- round(overview$sd_intrasp_onset_leaf_dormancy_weeks, digits=1)
overview = overview[,!(names(overview) %in% c("median_intrasp_onset_leaf_dormancy_weeks",
                                              "mean_intrasp_onset_leaf_dormancy_weeks",
                                              "sd_intrasp_onset_leaf_dormancy_weeks"))]
#----------------------------------
# overview$median_onset_turnover_weeks <- paste(round(overview$median_intrasp_onset_leaf_turnover_weeks, digits=1),
#                                            "(",
#                                            round(overview$mean_intrasp_onset_leaf_turnover_weeks, digits=1),
#                                            "±",
#                                            round(overview$sd_intrasp_onset_leaf_turnover_weeks, digits=1),
#                                            ")")
overview$median_onset_turnover_weeks <- round(overview$median_intrasp_onset_leaf_turnover_weeks, digits=1)
overview$intra_species_synchr_turnover_weeks <- round(overview$sd_intrasp_onset_leaf_turnover_weeks, digits=1)
overview = overview[,!(names(overview) %in% c("median_intrasp_onset_leaf_turnover_weeks",
                                              "mean_intrasp_onset_leaf_turnover_weeks",
                                              "sd_intrasp_onset_leaf_turnover_weeks"))]
#----------------------------------
#----------------------------------
# overview$mean_distance_onset_dormancy_weeks <- round(overview$mean_distance_onset_leaf_dormancy_weeks, digits=1)
# overview$minfreq_distance_onset_dormancy_weeks <- round(overview$mean_distance_onset_minfreq_leaf_dormancy_weeks, digits=1)
overview = overview[,!(names(overview) %in% c("mean_distance_onset_leaf_dormancy_weeks"))]
overview = overview[,!(names(overview) %in% c("mean_distance_onset_minfreq_leaf_dormancy_weeks"))]
#----------------------------------
# overview$mean_distance_onset_turnover_weeks <- round(overview$mean_distance_onset_leaf_turnover_weeks, digits=1)
# overview$minfreq_distance_onset_turnover_weeks <- round(overview$mean_distance_onset_minfreq_leaf_turnover_weeks, digits=1)
overview = overview[,!(names(overview) %in% c("mean_distance_onset_leaf_turnover_weeks"))]
overview = overview[,!(names(overview) %in% c("mean_distance_onset_minfreq_leaf_turnover_weeks"))]
#----------------------------------
#----------------------------------
overview$intra_annual_indiv_synchr_dormancy_weeks <- paste(round(overview$mean_synchrony_individuals_onset_leaf_dormancy_weeks, digits=1),
                                                           "±",
                                                           round(overview$sd_synchrony_individuals_onset_leaf_dormancy_weeks, digits=1))
overview = overview[,!(names(overview) %in% c("mean_synchrony_individuals_onset_leaf_dormancy_weeks",
                                              "sd_synchrony_individuals_onset_leaf_dormancy_weeks",
                                              "mean_nr_events_within_individuals_leaf_dormancy"))]
#----------------------------------
overview$intra_annual_indiv_synchr_turnover_weeks <- paste(round(overview$mean_synchrony_individuals_onset_leaf_turnover_weeks, digits=1),
                                                           "±",
                                                           round(overview$sd_synchrony_individuals_onset_leaf_turnover_weeks, digits=1))
overview = overview[,!(names(overview) %in% c("mean_synchrony_individuals_onset_leaf_turnover_weeks",
                                              "sd_synchrony_individuals_onset_leaf_turnover_weeks",
                                              "mean_nr_events_within_individuals_leaf_turnover"))]
#----------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------
# write to file
#--------------------------------------------------------------------
#--------------------------------------------------------------------
write.table(overview, "data/table_supplementary_information.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")
#--------------------------------------------------------------------
species_table_manuscript <- c("Scorodophloeus zenkeri",
                  "Panda oleosa",
                  "Staudtia kamerunensis",
                  "Synsepalum subcordatum",
                  "Prioria oxyphylla",
                  "Petersianthus macrocarpus",
                  "Irvingia grandifolia",
                  "Erythrophleum suaveolens",
                  "Autranella congolensis",
                  "Pericopsis elata")

table_manuscript <- overview %>%
  filter(species_full %in% species_table_manuscript)

write.table(table_manuscript, "data/table_manuscript.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")
#--------------------------------------------------------------------

