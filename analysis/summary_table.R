#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
# graphics.off()
#----- load required packages --------------------------------------------------#
library(tidyverse)


#-----------------------------------------------------------------------
# read in summary data (file made in summary_statistics_elizabeth.r)
#-----------------------------------------------------------------------
overview <- read.csv("data/species_meta_data.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)

# add full name and family
metadata <- read.csv("data/phenology_archives_species_long_format_20190626.csv",
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

# nr of events, if NA then zero
overview$total_nr_events_leaf_dormancy <- ifelse(is.na(overview$total_nr_events_leaf_dormancy),0,overview$total_nr_events_leaf_dormancy)
overview$total_nr_events_leaf_turnover <- ifelse(is.na(overview$total_nr_events_leaf_turnover),0,overview$total_nr_events_leaf_turnover)
#----------------------------------
# combine means with SD
# and then remove the original columns to get an overview of what's left
overview$duration_dormancy_weeks <- paste (round(overview$mean_duration_leaf_dormancy_weeks, digits=1),
                                           "±",
                                           round(overview$sd_duration_leaf_dormancy_weeks, digits=1))
overview$duration_dormancy_weeks <- ifelse(overview$duration_dormancy_weeks %in% "NA ± NA", "",overview$duration_dormancy_weeks)
overview = overview[,!(names(overview) %in% c("mean_duration_leaf_dormancy_weeks",
                                              "sd_duration_leaf_dormancy_weeks"))]
#----------------------------------
overview$duration_turnover_weeks <- paste (round(overview$mean_duration_leaf_turnover_weeks, digits=1),
                                           "±",
                                           round(overview$sd_duration_leaf_turnover_weeks, digits=1))
overview$duration_turnover_weeks <- ifelse(overview$duration_turnover_weeks %in% "NA ± NA", "",overview$duration_turnover_weeks)
overview = overview[,!(names(overview) %in% c("mean_duration_leaf_turnover_weeks",
                                              "sd_duration_leaf_turnover_weeks"))]
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

# dormancy x precip
tseries$corr_dormancy_precip <- paste(round(tseries$corr_dormancy_precip, digits=3),
                                          " (t",tseries$corr_dormancy_precip_timing,")", sep="")
tseries$corr_dormancy_precip <- ifelse(tseries$corr_dormancy_precip %in% "NA (tNA)", "",tseries$corr_dormancy_precip)
# turnover x precip
tseries$corr_turnover_precip <- paste(round(tseries$corr_turnover_precip, digits=3),
                                          " (t",tseries$corr_turnover_precip_timing,")", sep="")
tseries$corr_turnover_precip <- ifelse(tseries$corr_turnover_precip %in% "NA (tNA)", "",tseries$corr_turnover_precip)

# dormancy x sun
tseries$corr_dormancy_insol_JR <- paste(round(tseries$corr_dormancy_insol_JR, digits=3),
                                      " (t",tseries$corr_dormancy_insol_JR_timing,")", sep="")
tseries$corr_dormancy_insol_JR <- ifelse(tseries$corr_dormancy_insol_JR %in% "NA (tNA)", "",tseries$corr_dormancy_insol_JR)
# turnover x sun
tseries$corr_turnover_insol_JR <- paste(round(tseries$corr_turnover_insol_JR, digits=3),
                                      " (t",tseries$corr_turnover_insol_JR_timing,")", sep="")
tseries$corr_turnover_insol_JR <- ifelse(tseries$corr_turnover_insol_JR %in% "NA (tNA)", "",tseries$corr_turnover_insol_JR)

# dormancy x tmax
tseries$corr_dormancy_tmax_JR <- paste(round(tseries$corr_dormancy_tmax_JR, digits=3),
                                        " (t",tseries$corr_dormancy_tmax_JR_timing,")", sep="")
tseries$corr_dormancy_tmax_JR <- ifelse(tseries$corr_dormancy_tmax_JR %in% "NA (tNA)", "",tseries$corr_dormancy_tmax_JR)
# turnover x tmax
tseries$corr_turnover_tmax_JR <- paste(round(tseries$corr_turnover_tmax_JR, digits=3),
                                        " (t",tseries$corr_turnover_tmax_JR_timing,")", sep="")
tseries$corr_turnover_tmax_JR <- ifelse(tseries$corr_turnover_tmax_JR %in% "NA (tNA)", "",tseries$corr_turnover_tmax_JR)

tseries = tseries[,!(names(tseries) %in% c("corr_dormancy_precip_timing",
                                           "corr_turnover_precip_timing",
                                           "corr_dormancy_insol_JR_timing",
                                           "corr_turnover_insol_JR_timing",
                                           "corr_dormancy_tmax_JR_timing",
                                           "corr_turnover_tmax_JR_timing"
                                           ))]

overview <- merge(overview, tseries, by = "species_full", all.x = TRUE)


#--------------------------------------------------------------------
#--------------------------------------------------------------------
# set phenological label for species with NA in deciduousness
#--------------------------------------------------------------------
#--------------------------------------------------------------------
# star to add in table description that this is based on the actual data
# indicate specifically for genus level that this could potentially not be true for all species within the genus

tapply(overview$species_full, overview$deciduousness, length)

# clear
overview$deciduousness <- ifelse(overview$species_full %in% "Pericopsis elata", "deciduous*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Trichilia welwitschii", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Copaifera mildbraedii", "deciduous*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Tridesmostemon omphalocarpoides", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Omphalocarpum lecomteanum", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Fernandoa adolfi-friderici", "deciduous*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Tabernaemontana crassa", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Trichilia tessmannii", "deciduous*",overview$deciduousness)
# not so sure, limited data
overview$deciduousness <- ifelse(overview$species_full %in% "Trichilia gilletii", "evergreen* (?)",overview$deciduousness)
# not so sure, unclear phenological data
overview$deciduousness <- ifelse(overview$species_full %in% "Radlkofera calodendron", "evergreen* (?)",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Gilletiodendron mildbraedii", "evergreen* (?)",overview$deciduousness)

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

tapply(overview$species_full, overview$dec_label, length)

# test <- overview %>%
#   # filter(grepl("evergreen",deciduousness)) %>%
#   # filter(total_nr_events_leaf_turnover == 0 & total_nr_events_leaf_dormancy == 0)
#   # filter(total_nr_events_leaf_dormancy < 5)
#   filter(total_nr_events_leaf_turnover < 5 & total_nr_events_leaf_dormancy < 5)


#--------------------------------------------------------------------
#--------------------------------------------------------------------
# write to file
#--------------------------------------------------------------------
#--------------------------------------------------------------------
table1 <- overview[,(names(overview) %in% c("dec_label",
                                            "full_name",
                                            "family_Meise",
                                            "deciduousness",
                                            "ecology",
                                            "dbh_max_cm_literature",
                                            "basal_area_percentage",
                                            "nr_indiv",
                                            "site_years"))]
# correct column order
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
# correct column names
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
table2 <- overview[,(names(overview) %in% c("dec_label",
                                            "deciduousness",
                                            "species_full",
                                            "nr_indiv",
                                            "site_years",
                                            "percentage_site_years_with_leaf_dormancy",
                                            "total_nr_events_leaf_dormancy",
                                            "duration_dormancy_weeks",
                                            "median_onset_dormancy_weeks",
                                            "intra_species_synchr_dormancy_weeks",
                                            "intra_annual_indiv_synchr_dormancy_weeks",
                                            "corr_dormancy_precip",
                                            "corr_dormancy_insol_JR",
                                            "corr_dormancy_tmax_JR",
                                            "percentage_site_years_with_leaf_turnover",
                                            "total_nr_events_leaf_turnover",
                                            "duration_turnover_weeks",
                                            "median_onset_turnover_weeks",
                                            "intra_species_synchr_turnover_weeks",
                                            "intra_annual_indiv_synchr_turnover_weeks",
                                            "corr_turnover_precip",
                                            "corr_turnover_insol_JR",
                                            "corr_turnover_tmax_JR"))]
table2 <- table2[c("dec_label",
                   "deciduousness",
                   "species_full",
                   "nr_indiv",
                   "site_years",
                   "percentage_site_years_with_leaf_dormancy",
                   "total_nr_events_leaf_dormancy",
                   "duration_dormancy_weeks",
                   "median_onset_dormancy_weeks",
                   "intra_species_synchr_dormancy_weeks",
                   "intra_annual_indiv_synchr_dormancy_weeks",
                   "corr_dormancy_precip",
                   "corr_dormancy_insol_JR",
                   "corr_dormancy_tmax_JR",
                   "percentage_site_years_with_leaf_turnover",
                   "total_nr_events_leaf_turnover",
                   "duration_turnover_weeks",
                   "median_onset_turnover_weeks",
                   "intra_species_synchr_turnover_weeks",
                   "intra_annual_indiv_synchr_turnover_weeks",
                   "corr_turnover_precip",
                   "corr_turnover_insol_JR",
                   "corr_turnover_tmax_JR"
                   )]

write.table(table2, "data/SI_table2.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")
#--------------------------------------------------------------------
table2_dormancy <- table2[,(names(table2) %in% c("dec_label",
                                            "deciduousness",
                                            "species_full",
                                            "nr_indiv",
                                            "site_years",
                                            "percentage_site_years_with_leaf_dormancy",
                                            "total_nr_events_leaf_dormancy",
                                            "duration_dormancy_weeks",
                                            "median_onset_dormancy_weeks",
                                            "intra_species_synchr_dormancy_weeks",
                                            "intra_annual_indiv_synchr_dormancy_weeks",
                                            "corr_dormancy_precip",
                                            "corr_dormancy_insol_JR",
                                            "corr_dormancy_tmax_JR"
                                            ))]

colnames(table2_dormancy) <- c("dec_label",
                      "Leaf phenology",
                      "Species",
                      "Nr ind in JR",
                      "Total site-years in JR",
                      "% site-years with events",
                      "Total nr of events",
                      "Length event (weeks)",
                      "Timing onset event (WOY)",
                      "DSI intra-species",
                      "DSI intra-annual",
                      "corr precip",
                      "corr sun",
                      "corr tmax")

write.table(table2_dormancy, "data/SI_table2_dormancy.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")
#--------------------------------------------------------------------
table2_turnover <- table2[,(names(table2) %in% c("dec_label",
                                            "deciduousness",
                                            "species_full",
                                            "nr_indiv",
                                            "site_years",
                                            "percentage_site_years_with_leaf_turnover",
                                            "total_nr_events_leaf_turnover",
                                            "duration_turnover_weeks",
                                            "median_onset_turnover_weeks",
                                            "intra_species_synchr_turnover_weeks",
                                            "intra_annual_indiv_synchr_turnover_weeks",
                                            "corr_turnover_precip",
                                            "corr_turnover_insol_JR",
                                            "corr_turnover_tmax_JR"))]

colnames(table2_turnover) <- c("dec_label",
                               "Leaf phenology",
                               "Species",
                               "Nr ind in JR",
                               "Total site-years in JR",
                               "% site-years with events",
                               "Total nr of events",
                               "Length event (weeks)",
                               "Timing onset event (WOY)",
                               "DSI intra-species",
                               "DSI intra-annual",
                               "corr precip",
                               "corr sun",
                               "corr temp")

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
                              "Strombosia pustulata",
                              "Synsepalum subcordatum",
                              "Dacryodes osika",
                              "Quassia silvestris",
                              "Petersianthus macrocarpus",
                              "Irvingia grandifolia",
                              "Erythrophleum suaveolens",
                              "Antrocaryon nannanii",
                              "Pericopsis elata")



table_manuscript <- overview %>%
  filter(species_full %in% species_table_manuscript)

table1_ms <- table_manuscript[,(names(table_manuscript) %in% c("dec_label",
                                            "full_name",
                                            "family_Meise",
                                            "deciduousness",
                                            "ecology",
                                            "dbh_max_cm_literature",
                                            "basal_area_percentage",
                                            "nr_indiv",
                                            "site_years"))]

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
table2_dormancy_manuscript <- table2_dormancy %>%
  filter(Species %in% species_table_manuscript)

table2_turnover_manuscript <- table2_turnover %>%
  filter(Species %in% species_table_manuscript)


write.table(table2_dormancy_manuscript, "data/table2_manuscript_dormancy.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")

write.table(table2_turnover_manuscript, "data/table2_manuscript_turnover.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")

