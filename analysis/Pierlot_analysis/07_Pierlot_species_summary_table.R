#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
graphics.off()
#----- load required packages --------------------------------------------------#
library(tidyverse)


#-----------------------------------------------------------------------
# read in summary data (file made in summary_statistics_elizabeth.r)
#-----------------------------------------------------------------------
overview <- read.csv("data/Pierlot_summ_species_pheno_characteristics.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)

# add full name and family
metadata <- read.csv("data-raw/phenology_archives_species_long_format_20200324.csv",
                     header = TRUE, sep = ",")
metadata$species_full <- paste(metadata$genus_Meise, metadata$species_Meise)
metadata$full_name <- paste(metadata$genus_Meise, metadata$species_Meise, metadata$author_Meise)
metadata <- metadata %>%
  dplyr::select(species_full,
                full_name,
                family_Meise) %>%
  dplyr::rename(family = family_Meise)

metadata <- metadata[!duplicated(metadata), ]
# test <- metadata %>%
#   filter(species_full %in% "Albizia adianthifolia")
metadata$species_full <- ifelse(metadata$full_name %in% "Albizia adianthifolia (De Wild.& T. Durand) Villiers","NA",metadata$species_full)
overview <- merge(overview, metadata, by = "species_full", all.x = TRUE)

# remove some columns you don't want in the overview table
overview <- overview %>%
  dplyr::select(-c(abundance_rel,
                   basal_area_species_all,
                   dbh_max_cm_literature,
                   height_m_literature,
                   Ref))

# round up some numbers
overview <- overview %>%
  mutate(percentage_site_years_with_phenophase = round(ratio_site_years_with_phenophase, digits=3) *100,
         basal_area_percentage = round(basal_area_weighted_site, digits=2)) %>%
  dplyr::select(-ratio_site_years_with_phenophase)


# print NAs as empty cells
overview$ecology <- ifelse(is.na(overview$ecology), "",overview$ecology)
overview$DBH_class_max_all <- ifelse(is.na(overview$DBH_class_max_all), "",overview$DBH_class_max_all)
# nr of events, if NA then zero
overview$total_nr_events <- ifelse(is.na(overview$total_nr_events), 0, overview$total_nr_events)
#----------------------------------
# combine means with SD
# and then remove the original columns to get an overview of what's left
overview$phenophase_length_weeks <- paste (round(overview$mean_phenophase_length_weeks, digits=1),
                                           "±",
                                           round(overview$sd_phenophase_length_weeks, digits=1))
overview$phenophase_length_weeks <- ifelse(overview$phenophase_length_weeks %in% "NA ± NA", "",overview$phenophase_length_weeks)
overview <- overview %>%
  dplyr::select(-c(mean_phenophase_length_weeks,
                   sd_phenophase_length_weeks))
#----------------------------------
overview$median_onset_weeks <- round(overview$median_intrasp_onset_weeks, digits=1)
overview$median_onset_weeks <- ifelse(is.na(overview$median_onset_weeks), "",overview$median_onset_weeks)
overview$sd_intrasp_onset_weeks <- round(overview$sd_intrasp_onset_weeks, digits=1)
overview$sd_intrasp_onset_weeks <- ifelse(is.na(overview$sd_intrasp_onset_weeks), "",overview$sd_intrasp_onset_weeks)
overview$sd_intrasp_onset_weeks <- ifelse(overview$sd_intrasp_onset_weeks %in% "0", "",overview$sd_intrasp_onset_weeks)

overview <- overview %>%
  dplyr::select(-c(median_intrasp_onset_weeks,
                   mean_intrasp_onset_weeks))
#----------------------------------
overview$intra_ind_onset_weeks <- paste(round(overview$mean_intra_ind_onset_weeks, digits=1),
                                        "±",
                                        round(overview$sd_intra_ind_onset_weeks, digits=1))
overview$intra_ind_onset_weeks <- ifelse(overview$intra_ind_onset_weeks %in% "NA ± NA", "",overview$intra_ind_onset_weeks)
overview$intra_ind_onset_weeks <- ifelse(overview$intra_ind_onset_weeks %in% "0 ± NA", "",overview$intra_ind_onset_weeks)
overview$intra_ind_onset_weeks <- ifelse(overview$intra_ind_onset_weeks %in% "0 ± 0", "",overview$intra_ind_onset_weeks)

overview <- overview %>%
  dplyr::select(-c(mean_intra_ind_onset_weeks,
                   sd_intra_ind_onset_weeks))

#--------------------------------------------------------------------
#--------------------------------------------------------------------
# now include the timeseries correlations with precipitation
#--------------------------------------------------------------------
#--------------------------------------------------------------------

tseries <- read.csv("data/Pierlot_timeseries_crosscorrelations.csv",
                    header = TRUE,
                    sep = ",",
                    stringsAsFactors = FALSE)

# round to 3 digits and add timing to corr

# dormancy x precip
tseries$corr_precip <- paste(round(tseries$corr_precip, digits=3),
                             " (t",tseries$corr_precip_timing,")", sep="")
tseries$corr_precip <- ifelse(tseries$corr_precip %in% "NA (tNA)", "",tseries$corr_precip)


# dormancy x sun
tseries$corr_sun <- paste(round(tseries$corr_sun, digits=3),
                          " (t",tseries$corr_sun_timing,")", sep="")
tseries$corr_sun <- ifelse(tseries$corr_sun %in% "NA (tNA)", "",tseries$corr_sun)


# dormancy x tmax
tseries$corr_temp <- paste(round(tseries$corr_temp, digits=3),
                           " (t",tseries$corr_temp_timing,")", sep="")
tseries$corr_temp <- ifelse(tseries$corr_temp %in% "NA (tNA)", "",tseries$corr_temp)


tseries = tseries[,!(names(tseries) %in% c("corr_precip_timing",
                                           "corr_sun_timing",
                                           "corr_temp_timing"
))]
#--------------------------------------------------------------------
overview <- merge(overview, tseries, by = c("species_full","phenophase"), all.x = TRUE)
rm(tseries)
#--------------------------------------------------------------------



#--------------------------------------------------------------------
#--------------------------------------------------------------------
# fourier analysis
#--------------------------------------------------------------------
#--------------------------------------------------------------------
four_dorm <- read.csv("data/Pierlot_fourier_leaf_dormancy.csv",
                      header = TRUE,
                      sep = ",",
                      stringsAsFactors = FALSE)
four_turn <- read.csv("data/Pierlot_fourier_leaf_turnover.csv",
                      header = TRUE,
                      sep = ",",
                      stringsAsFactors = FALSE)
dffour <- rbind(four_dorm, four_turn)
dffour <- dffour %>%
  select(-sig_95)

dffour$cycle_category <- ifelse(is.na(dffour$cycle_category), "",dffour$cycle_category)
dffour$cycle_dom_months <- round(dffour$cycle_dom_months, digits=1)
dffour$cycle_dom_months <- ifelse(is.na(dffour$cycle_dom_months), "",dffour$cycle_dom_months)
#--------------------------------------------------------------------
overview <- merge(overview, dffour, by = c("species_full","phenophase"), all.x = TRUE)
rm(four_dorm, four_turn, dffour)
#--------------------------------------------------------------------


#--------------------------------------------------------------------
# set phenological label
# to sort easily for making the eventual tables -. uniform label deciduous or evergreen
overview$dec_label <- ifelse(overview$deciduousness %in% c("deciduous","deciduous*","deciduous**"), "dec",
                             ifelse(overview$deciduousness %in% c("evergreen","evergreen*","evergreen**"), "ever", "NA"))


#--------------------------------------------------------------------
#--------------------------------------------------------------------
# write to file
#--------------------------------------------------------------------
#--------------------------------------------------------------------
table1 <- overview %>%
  filter(phenophase == "leaf_dormancy") %>%
  dplyr::select(dec_label,
                full_name,
                family,
                deciduousness,
                ecology,
                DBH_class_max_all,
                basal_area_percentage,
                nr_indiv,
                site_years,
                start_year,
                end_year)

write.table(table1, "data/Pierlot_SI_table1.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")

#--------------------------------------------------------------------
#--------------------------------------------------------------------
table2 <- overview %>%
  filter(phenophase == "leaf_dormancy") %>%
  dplyr::select(phenophase,
                dec_label,
                species_full,
                nr_indiv,
                site_years,
                site_years_with_phenophase, #percentage_site_years_with_phenophase,
                total_nr_events,
                cycle_category,
                cycle_dom_months,
                phenophase_length_weeks,
                median_onset_weeks,
                sd_intrasp_onset_weeks,
                intra_ind_onset_weeks,
                corr_precip,
                corr_sun,
                corr_temp,
                groups)

write.table(table2, "data/Pierlot_SI_table2_dormancy.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")
#--------------------------------------------------------------------
table3 <- overview %>%
  filter(phenophase == "leaf_turnover") %>%
  dplyr::select(phenophase,
                dec_label,
                species_full,
                nr_indiv,
                site_years,
                site_years_with_phenophase, #percentage_site_years_with_phenophase,
                total_nr_events,
                cycle_category,
                cycle_dom_months,
                phenophase_length_weeks,
                median_onset_weeks,
                sd_intrasp_onset_weeks,
                intra_ind_onset_weeks,
                corr_precip,
                corr_sun,
                corr_temp,
                groups)

write.table(table3, "data/Pierlot_SI_table3_turnover.csv",
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
                              "Allophylus africanus",
                              "Erythrophleum suaveolens",
                              "Pterocarpus soyauxii",
                              "Pericopsis elata")

table_a <- table2 %>%
  filter(species_full %in% species_table_manuscript)
table_b <- table3 %>%
  filter(species_full %in% species_table_manuscript)
table_ms <- rbind(table_a, table_b)

write.table(table_ms, "data/Pierlot_table1_manuscript.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")


