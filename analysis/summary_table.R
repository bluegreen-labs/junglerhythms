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

# add full name and family
metadata <- read.csv("data/phenology_archives_species_long_format_20190619.csv",
                     header = TRUE, sep = ",")
metadata$species_full <- paste(metadata$genus_Meise, metadata$species_Meise)
metadata$full_name <- paste(metadata$genus_Meise, metadata$species_Meise, metadata$author_Meise)
metadata = metadata[,(names(metadata) %in% c("species_full",
                                              "full_name",
                                              "family_Meise"))]
metadata <- metadata[!duplicated(metadata), ]
metadata$species_full <- ifelse(metadata$full_name %in% "Albizia adianthifolia (J.F.Gmel.) C.A.Sm.","NA",metadata$species_full)
overview <- merge(overview, metadata, by = "species_full", all.x = TRUE)


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

#
overview$ecology <- ifelse(is.na(overview$ecology), "",overview$ecology)
overview$dbh_max_cm_literature <- ifelse(is.na(overview$dbh_max_cm_literature), "",overview$dbh_max_cm_literature)

#----------------------------------
# combine means with SD and (n = )
# and then remove the original columns to get an overview of what's left
overview$duration_dormancy_weeks <- paste (round(overview$mean_duration_leaf_dormancy_weeks, digits=1),
                                           "±",
                                           round(overview$sd_duration_leaf_dormancy_weeks, digits=1),
                                           "(n =",
                                           overview$total_nr_events_leaf_dormancy,
                                           ")")
overview$duration_dormancy_weeks <- ifelse(overview$duration_dormancy_weeks %in% "NA ± NA (n = NA )", "",overview$duration_dormancy_weeks)
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
overview$duration_turnover_weeks <- ifelse(overview$duration_turnover_weeks %in% "NA ± NA (n = NA )", "",overview$duration_turnover_weeks)
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
overview$median_onset_dormancy_weeks <- ifelse(is.na(overview$median_onset_dormancy_weeks), "",overview$median_onset_dormancy_weeks)
overview$intra_species_synchr_dormancy_weeks <- round(overview$sd_intrasp_onset_leaf_dormancy_weeks, digits=1)
overview$intra_species_synchr_dormancy_weeks <- ifelse(is.na(overview$intra_species_synchr_dormancy_weeks), "",overview$intra_species_synchr_dormancy_weeks)
overview$intra_species_synchr_dormancy_weeks <- ifelse(overview$intra_species_synchr_dormancy_weeks %in% "0", "",overview$intra_species_synchr_dormancy_weeks)

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
overview$median_onset_turnover_weeks <- ifelse(is.na(overview$median_onset_turnover_weeks), "",overview$median_onset_turnover_weeks)
overview$intra_species_synchr_turnover_weeks <- round(overview$sd_intrasp_onset_leaf_turnover_weeks, digits=1)
overview$intra_species_synchr_turnover_weeks <- ifelse(is.na(overview$intra_species_synchr_turnover_weeks), "",overview$intra_species_synchr_turnover_weeks)
overview$intra_species_synchr_turnover_weeks <- ifelse(overview$intra_species_synchr_turnover_weeks %in% "0", "",overview$intra_species_synchr_turnover_weeks)
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
overview$intra_annual_indiv_synchr_dormancy_weeks <- ifelse(overview$intra_annual_indiv_synchr_dormancy_weeks %in% "NA ± NA", "",overview$intra_annual_indiv_synchr_dormancy_weeks)
overview$intra_annual_indiv_synchr_dormancy_weeks <- ifelse(overview$intra_annual_indiv_synchr_dormancy_weeks %in% "0 ± NA", "",overview$intra_annual_indiv_synchr_dormancy_weeks)
overview$intra_annual_indiv_synchr_dormancy_weeks <- ifelse(overview$intra_annual_indiv_synchr_dormancy_weeks %in% "0 ± 0", "",overview$intra_annual_indiv_synchr_dormancy_weeks)
overview = overview[,!(names(overview) %in% c("mean_synchrony_individuals_onset_leaf_dormancy_weeks",
                                              "sd_synchrony_individuals_onset_leaf_dormancy_weeks",
                                              "mean_nr_events_within_individuals_leaf_dormancy"))]
#----------------------------------
overview$intra_annual_indiv_synchr_turnover_weeks <- paste(round(overview$mean_synchrony_individuals_onset_leaf_turnover_weeks, digits=1),
                                                           "±",
                                                           round(overview$sd_synchrony_individuals_onset_leaf_turnover_weeks, digits=1))
overview$intra_annual_indiv_synchr_turnover_weeks <- ifelse(overview$intra_annual_indiv_synchr_turnover_weeks %in% "NA ± NA", "",overview$intra_annual_indiv_synchr_turnover_weeks)
overview$intra_annual_indiv_synchr_turnover_weeks <- ifelse(overview$intra_annual_indiv_synchr_turnover_weeks %in% "0 ± NA", "",overview$intra_annual_indiv_synchr_turnover_weeks)
overview$intra_annual_indiv_synchr_turnover_weeks <- ifelse(overview$intra_annual_indiv_synchr_turnover_weeks %in% "0 ± 0", "",overview$intra_annual_indiv_synchr_turnover_weeks)
overview = overview[,!(names(overview) %in% c("mean_synchrony_individuals_onset_leaf_turnover_weeks",
                                              "sd_synchrony_individuals_onset_leaf_turnover_weeks",
                                              "mean_nr_events_within_individuals_leaf_turnover"))]
#----------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------
# now include the timeseries correlations with precipitation
#--------------------------------------------------------------------
#--------------------------------------------------------------------

tseries <- read.csv("data/timeseries_correlations.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)

# round to 3 digits and add timing to corr
tseries$corr_dormancy_precip <- round(tseries$corr_dormancy_precip, digits=3)
tseries$corr_dormancy_precip <- ifelse(is.na(tseries$corr_dormancy_precip), "",tseries$corr_dormancy_precip)

tseries$ccf_dormancy_precip_lead <- paste(round(tseries$ccf_dormancy_precip_lead, digits=3),
                                          " (t",tseries$ccf_dormancy_precip_lead_timing,")", sep="")
tseries$ccf_dormancy_precip_lead <- ifelse(tseries$ccf_dormancy_precip_lead %in% "NA (tNA)", "",tseries$ccf_dormancy_precip_lead)

tseries$ccf_dormancy_precip_lag <- paste(round(tseries$ccf_dormancy_precip_lag, digits=3),
                                         " (t+",tseries$ccf_dormancy_precip_lag_timing,")", sep="")
tseries$ccf_dormancy_precip_lag <- ifelse(tseries$ccf_dormancy_precip_lag %in% "NA (t+NA)", "",tseries$ccf_dormancy_precip_lag)

tseries$ccf_dormancy_anom_lead <- paste(round(tseries$ccf_dormancy_anom_lead, digits=3),
                                        " (t",tseries$ccf_dormancy_anom_lead_timing,")", sep="")
tseries$ccf_dormancy_anom_lead <- ifelse(tseries$ccf_dormancy_anom_lead %in% "NA (tNA)", "",tseries$ccf_dormancy_anom_lead)


tseries$corr_turnover_precip <- round(tseries$corr_turnover_precip, digits=3)
tseries$corr_turnover_precip <- ifelse(is.na(tseries$corr_turnover_precip), "",tseries$corr_turnover_precip)

tseries$ccf_turnover_precip_lead <- paste(round(tseries$ccf_turnover_precip_lead, digits=3),
                                          " (t",tseries$ccf_turnover_precip_lead_timing,")", sep="")
tseries$ccf_turnover_precip_lead <- ifelse(tseries$ccf_turnover_precip_lead %in% "NA (tNA)", "",tseries$ccf_turnover_precip_lead)

tseries$ccf_turnover_precip_lag <- paste(round(tseries$ccf_turnover_precip_lag, digits=3),
                                         " (t+",tseries$ccf_turnover_precip_lag_timing,")", sep="")
tseries$ccf_turnover_precip_lag <- ifelse(tseries$ccf_turnover_precip_lag %in% "NA (t+NA)", "",tseries$ccf_turnover_precip_lag)

tseries$ccf_turnover_anom_lead <- paste(round(tseries$ccf_turnover_anom_lead, digits=3),
                                        " (t",tseries$ccf_turnover_anom_lead_timing,")", sep="")
tseries$ccf_turnover_anom_lead <- ifelse(tseries$ccf_turnover_anom_lead %in% "NA (tNA)", "",tseries$ccf_turnover_anom_lead)

tseries = tseries[,!(names(tseries) %in% c("ccf_dormancy_precip_lead_timing",
                                           "ccf_dormancy_precip_lag_timing",
                                           "ccf_dormancy_anom_lead_timing",
                                           "ccf_turnover_precip_lead_timing",
                                           "ccf_turnover_precip_lag_timing",
                                           "ccf_turnover_anom_lead_timing"
                                           ))]

overview <- merge(overview, tseries, by = "species_full", all.x = TRUE)
#--------------------------------------------------------------------
#--------------------------------------------------------------------
# set phenological label for species with NA in deciduousness
#--------------------------------------------------------------------
#--------------------------------------------------------------------
# star to add in table description that this is based on the actual data
# indicate specifically for genus level that this could potentially not be true for all species within the genus

# clear
overview$deciduousness <- ifelse(overview$species_full %in% "Pericopsis elata", "deciduous*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Diospyros sp.", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Beilschmiedia sp.", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Macaranga sp.", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Trichilia welwitschii", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Chytranthus sp.", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Copaifera mildbraedii", "deciduous*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Tridesmostemon omphalocarpoides", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Omphalocarpum lecomteanum", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Fernandoa adolfi-friderici", "deciduous*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Cola sp.", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Tessmannia sp.", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Irvingia sp.", "deciduous*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Tabernaemontana crassa", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Trichilia tessmannii", "deciduous*",overview$deciduousness)
# not so sure, limited data
overview$deciduousness <- ifelse(overview$species_full %in% "Trichilia sp.", "deciduous* (?)",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Trichilia gilletii", "evergreen* (?)",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Drypetes sp.", "evergreen* (?)",overview$deciduousness)
# not so sure, unclear phenological data
overview$deciduousness <- ifelse(overview$species_full %in% "Radlkofera calodendron", "evergreen* (?)",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Gilletiodendron mildbraedii", "evergreen* (?)",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Homalium sp.", "deciduous* (?)",overview$deciduousness)

## two stars, in literature found as evergreen or (sometimes) deciduous
## selected a class based on the actual data
overview$deciduousness <- ifelse(overview$species_full %in% "Celtis mildbraedii", "evergreen**",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Combretum lokele", "deciduous**",overview$deciduousness)

overview$deciduousness <- ifelse(overview$species_full %in% "Homalium africanum", "evergreen**",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Quassia silvestris", "evergreen**",overview$deciduousness)

# not so sure, unclear phenological data
overview$deciduousness <- ifelse(overview$species_full %in% "Homalium longistylum", "evergreen** (?)",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Irvingia gabonensis", "deciduous** (?)",overview$deciduousness)

# to sort easily for making the eventual tables -. uniform label deciduous or evergreen
overview$dec_label <- ifelse(overview$deciduousness %in% c("deciduous","deciduous*","deciduous**","deciduous* (?)","deciduous** (?)"), "dec",
                             ifelse(overview$deciduousness %in% c("evergreen","evergreen*","evergreen**","evergreen* (?)","evergreen** (?)"), "ever", "NA"))


#--------------------------------------------------------------------
#--------------------------------------------------------------------
# write to file
#--------------------------------------------------------------------
#--------------------------------------------------------------------
table1 <- overview[,c(2:9,28)]
table1 <- table1[c("dec_label",
                   "full_name",
                   "family_Meise",
                   "deciduousness",
                   "ecology",
                   "dbh_max_cm_literature",
                   "basal_area_percentage",
                   "nr_indiv",
                   "site_years"
                   )]
colnames(table1) <- c("dec_label",
                      "Species",
                      "Family",
                      "Leaf phenology",
                      "Ecology",
                      "DBH max",
                      "% BA Yangambi",
                      "Nr ind in JR",
                      "Total site-years in JR"
                      )

write.table(table1, "data/SI_table1.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")

#--------------------------------------------------------------------
#--------------------------------------------------------------------
table2 <- overview[,c(1:3,6,10:28)]
table2 <- table2[c("dec_label",
                   "deciduousness",
                   "species_full",
                   "nr_indiv",
                   "site_years",
                   "percentage_site_years_with_leaf_dormancy",
                   "duration_dormancy_weeks",
                   "median_onset_dormancy_weeks",
                   "intra_species_synchr_dormancy_weeks",
                   "intra_annual_indiv_synchr_dormancy_weeks",
                   "corr_dormancy_precip",
                   "ccf_dormancy_precip_lead",
                   "ccf_dormancy_precip_lag",
                   "ccf_dormancy_anom_lead",
                   "percentage_site_years_with_leaf_turnover",
                   "duration_turnover_weeks",
                   "median_onset_turnover_weeks",
                   "intra_species_synchr_turnover_weeks",
                   "intra_annual_indiv_synchr_turnover_weeks",
                   "corr_turnover_precip",
                   "ccf_turnover_precip_lead",
                   "ccf_turnover_precip_lag",
                   "ccf_turnover_anom_lead"
                   )]

write.table(table2, "data/SI_table2.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")
#--------------------------------------------------------------------
table2_dormancy <- table2[,c(1:5,6:14)]

colnames(table2_dormancy) <- c("dec_label",
                      "Leaf phenology",
                      "Species",
                      "Nr ind in JR",
                      "Total site-years in JR",
                      "% site-years with events",
                      "Length event (weeks)",
                      "Timing onset event (WOY)",
                      "DSI intra-species",
                      "DSI intra-annual",
                      "corr precip",
                      "CCF precip leading event",
                      "CCF precip lag after event",
                      "CCF anomaly in precip leading event")

write.table(table2_dormancy, "data/SI_table2_dormancy.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")
#--------------------------------------------------------------------
table2_turnover <- table2[,c(1:5,15:23)]

colnames(table2_turnover) <- c("dec_label",
                               "Leaf phenology",
                               "Species",
                               "Nr ind in JR",
                               "Total site-years in JR",
                               "% site-years with events",
                               "Length event (weeks)",
                               "Timing onset event (WOY)",
                               "DSI intra-species",
                               "DSI intra-annual",
                               "corr precip",
                               "CCF precip leading event",
                               "CCF precip lag after event",
                               "CCF anomaly in precip leading event")

write.table(table2_turnover, "data/SI_table2_turnover.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")
#--------------------------------------------------------------------
#--------------------------------------------------------------------
# table for in the manuscript
#--------------------------------------------------------------------
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

table1_ms <- table_manuscript[,c(2:9,28)]

table1_ms <- table1_ms[c("dec_label",
                   "full_name",
                   "family_Meise",
                   "deciduousness",
                   "ecology",
                   "dbh_max_cm_literature",
                   "basal_area_percentage",
                   "nr_indiv",
                   "site_years"
                   )]
colnames(table1_ms) <- c("dec_label",
                      "Species",
                      "Family",
                      "Leaf phenology",
                      "Ecology",
                      "DBH max",
                      "% BA Yangambi",
                      "Nr ind in JR",
                      "Total site-years in JR"
                      )

write.table(table1_ms, "data/table1_manuscript.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")
#--------------------------------------------------------------------
table2_ms <- table_manuscript[,c(1:3,6,10:28)]
table2_ms <- table2_ms[c("dec_label",
                         "deciduousness",
                   "species_full",
                   "nr_indiv",
                   "site_years",
                   "percentage_site_years_with_leaf_dormancy",
                   "duration_dormancy_weeks",
                   "median_onset_dormancy_weeks",
                   "intra_species_synchr_dormancy_weeks",
                   "intra_annual_indiv_synchr_dormancy_weeks",
                   "corr_dormancy_precip",
                   "ccf_dormancy_precip_lead",
                   "ccf_dormancy_precip_lag",
                   "ccf_dormancy_anom_lead",
                   "percentage_site_years_with_leaf_turnover",
                   "duration_turnover_weeks",
                   "median_onset_turnover_weeks",
                   "intra_species_synchr_turnover_weeks",
                   "intra_annual_indiv_synchr_turnover_weeks",
                   "corr_turnover_precip",
                   "ccf_turnover_precip_lead",
                   "ccf_turnover_precip_lag",
                   "ccf_turnover_anom_lead"
                   )]

table2_ms_dormancy <- table2_ms[,c(1:5,6:14)]
colnames(table2_ms_dormancy) <- c("dec_label",
                               "Leaf phenology",
                               "Species",
                               "Nr ind in JR",
                               "Total site-years in JR",
                               "% site-years with events",
                               "Length event (weeks)",
                               "Timing onset event (WOY)",
                               "DSI intra-species",
                               "DSI intra-annual",
                               "corr precip",
                               "CCF precip leading event",
                               "CCF precip lag after event",
                               "CCF anomaly in precip leading event")
table2_ms_turnover <- table2_ms[,c(1:5,15:23)]
colnames(table2_ms_turnover) <- c("dec_label",
                                  "Leaf phenology",
                                  "Species",
                                  "Nr ind in JR",
                                  "Total site-years in JR",
                                  "% site-years with events",
                                  "Length event (weeks)",
                                  "Timing onset event (WOY)",
                                  "DSI intra-species",
                                  "DSI intra-annual",
                                  "corr precip",
                                  "CCF precip leading event",
                                  "CCF precip lag after event",
                                  "CCF anomaly in precip leading event")

write.table(table2_ms_dormancy, "data/table2_manuscript_dormancy.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")

write.table(table2_ms_turnover, "data/table2_manuscript_turnover.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")

