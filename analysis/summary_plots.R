

# read in summary data
overview <- read.csv("data/species_meta_data.csv",
                    header = TRUE,
                    sep = ",",
                    stringsAsFactors = FALSE)

# first remove if too few events to calculate index (minimum 5 events within a species)
overview$sd_intrasp_onset_leaf_dormancy_weeks <- ifelse(overview$site_years_with_leaf_dormancy < 5, NA, overview$sd_intrasp_onset_leaf_dormancy_weeks)
overview$sd_intrasp_onset_leaf_turnover_weeks <- ifelse(overview$site_years_with_leaf_turnover < 5, NA, overview$sd_intrasp_onset_leaf_turnover_weeks)
# first remove if too few events to calculate index (minimum 2 events per individual)
overview$mean_synchrony_individuals_onset_leaf_dormancy_weeks <- ifelse(overview$mean_nr_events_within_individuals_leaf_dormancy < 2, NA, overview$mean_synchrony_individuals_onset_leaf_dormancy_weeks)
overview$mean_synchrony_individuals_onset_leaf_turnover_weeks <- ifelse(overview$mean_synchrony_individuals_onset_leaf_turnover_weeks < 2, NA, overview$mean_synchrony_individuals_onset_leaf_turnover_weeks)


dfdec <- overview %>%
  filter(deciduousness == "deciduous" | deciduousness == "evergreen")
dfecol <- overview %>%
  filter(ecology == "shade" | ecology == "sun")

test <- dfecol %>%
  filter(ecology == "sun")


par(mfrow=c(4,4),
    mar=c(6,5,2,2))
# ratio eventyears:siteyears
boxplot(dfdec$ratio_site_years_with_leaf_dormancy ~ dfdec$deciduousness,ylab = "ratio years with canopy dormancy")
boxplot(dfdec$ratio_site_years_with_leaf_turnover ~ dfdec$deciduousness,ylab = "ratio years with canopy turnover")
boxplot(dfecol$ratio_site_years_with_leaf_dormancy ~ dfecol$ecology,ylab = "ratio years with canopy dormancy")
boxplot(dfecol$ratio_site_years_with_leaf_turnover ~ dfecol$ecology,ylab = "ratio years with canopy turnover")

# duration event
boxplot(dfdec$mean_duration_leaf_dormancy_weeks ~ dfdec$deciduousness,ylab = "mean duration canopy dormancy (weeks)")
boxplot(dfdec$mean_duration_leaf_turnover_weeks ~ dfdec$deciduousness,ylab = "mean duration canopy turnover (weeks)")
boxplot(dfecol$mean_duration_leaf_dormancy_weeks ~ dfecol$ecology,ylab = "mean duration canopy dormancy (weeks)")
boxplot(dfecol$mean_duration_leaf_turnover_weeks ~ dfecol$ecology,ylab = "mean duration canopy turnover (weeks)")

# intraspecific synchrony index
boxplot(dfdec$sd_intrasp_onset_leaf_dormancy_weeks ~ dfdec$deciduousness,ylab = "intraspecific synchrony index leaf dormancy (weeks)")
boxplot(dfdec$sd_intrasp_onset_leaf_turnover_weeks ~ dfdec$deciduousness,ylab = "intraspecific synchrony index leaf turnover (weeks)")
boxplot(dfecol$sd_intrasp_onset_leaf_dormancy_weeks ~ dfecol$ecology,ylab = "intraspecific synchrony index leaf dormancy (weeks)")
boxplot(dfecol$sd_intrasp_onset_leaf_turnover_weeks ~ dfecol$ecology,ylab = "intraspecific synchrony index leaf turnover (weeks)")

# intra-individual synchrony index
boxplot(dfdec$mean_synchrony_individuals_onset_leaf_dormancy_weeks ~ dfdec$deciduousness,ylab = "intra-individual synchrony index leaf dormancy (weeks)")
boxplot(dfdec$mean_synchrony_individuals_onset_leaf_turnover_weeks ~ dfdec$deciduousness,ylab = "intra-individual synchrony index leaf turnover (weeks)")
boxplot(dfecol$mean_synchrony_individuals_onset_leaf_dormancy_weeks ~ dfecol$ecology,ylab = "intra-individual synchrony index leaf dormancy (weeks)")
boxplot(dfecol$mean_synchrony_individuals_onset_leaf_turnover_weeks ~ dfecol$ecology,ylab = "intra-individual synchrony index leaf turnover (weeks)")

