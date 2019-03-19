# compile summary statistics
library(tidyverse)
library(ggplot2)
library(ggthemes)

#----------------------------------------------------------------------
#--------   Phenology data - species correction Meise   ---------------
#----------------------------------------------------------------------
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
df <- df[which(df$value != 0),]
df$join_id <- paste0("R",df$image,"-",df$image_row)

metadata <- read.csv("data/phenology_archives_species_long_format_20190314.csv",
                     header = TRUE, sep = ",")
metadata$join_id <- paste(metadata$image,metadata$row, sep = "-")

# test merge the two tables based upon the unique join ID
df <- merge(df, metadata, by = c("join_id","id"), all.x = TRUE)
df$species_full <- paste(df$genus_Meise, df$species_Meise)
#----------------------------------------------------------------------

# # read in census data and convert basal area
# census <- read.csv2("data/YGB_ForestPlotsNET.csv",
#                      header = TRUE,
#                      sep = ",",
#                      stringsAsFactors = FALSE)
#
# # calculate mean basal area across all species across all
# # mixed plots
# census <- census %>%
#   filter(grepl("MIX",Plot)) %>%
#   group_by(Plot, Species) %>%
#   summarize(basal_area = sum(pi*(C2DBH4/2000)^2, na.rm = TRUE)) %>%
#   group_by(Species) %>%
#   summarize(basal_area = mean(basal_area))
#
# # replace empty slots with Unknown for consistency
# df$genus <- ifelse(df$genus == "","Unknown",df$genus)
# df$species <- ifelse(df$species == "","Unknown",df$species)

# read in census data
census <- read.csv2("data/yangambi_mixed_forest_species_list.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)
colnames(census)[1] <- 'species_full'

# read in trait data from Steven Janssens
traits <- read.csv2("data/Dataset_traits_African_trees.csv",
                    header = TRUE,
                    sep = ",",
                    stringsAsFactors = FALSE)

# summary statistics: start year, end year, site years, site years with events, nr of individuals
# this is now per id
# counting unique year - id combinations for either
# the full or reduced (phenophase based) subsets
phenophase_selected <- c("leaf_turnover","leaf_dormancy")
# loop over all species (full name)
out <- lapply(unique(df$species_full), function(species_selected){

  # list of all observed data regardless of phenophase
  # get then number of unique years
  ss_full <- subset(df, species_full == species_selected,
                    select = c('year','id'))
  nr_years_full <- length(unique(paste(ss_full$year, ss_full$id)))
  nr_indiv <- length(unique(ss_full$id))

  start_year <- min(ss_full$year)
  end_year <- max(ss_full$year)

  # subset based on phenophase of interest + species
  # get the number of unique years (with observations)
  ss_phen <- subset(df,
                    species_full == species_selected &
                      phenophase == phenophase_selected,
                    select = c('year','id'))
  nr_years_phen <- length(unique(paste(ss_phen$year, ss_phen$id)))

  # return the ratio
  ratio <- nr_years_phen/nr_years_full

  return(data.frame(
    "species_full" = species_selected,
    # "phenophase"= phenophase_selected,
    "nr_indiv" = nr_indiv,
    "start_year" = start_year,
    "end_year"= end_year,
    "nr_years_full" = nr_years_full,
    "nr_years_phen" = nr_years_phen,
    "ratio" = ratio))
})

# bind everything row wise
out <- do.call("rbind", out)


# merge and drop merge column
out <- merge(out, census, by = "species_full", all.x = TRUE)
out <- merge(out, traits, by = "species_full", all.x = TRUE)
# test <- test %>% select(-Species)

# remove rows (species) not included in the Yangambi mixed forest census
final <- out[!(is.na(out$BAperc)),]

# write to file
write.table(final, "data/species_meta_data.csv",
          quote = FALSE,
          col.names = TRUE,
          row.names = FALSE,
          sep = ",")



# sum the basal area of the observed species
# if >30 good coverage
# message(sum(final$BA, na.rm = TRUE))

hist1 <- ggplot(final, aes(nr_years_full)) +
  geom_histogram(bins = 15) +
  labs(x = "Site Years",
       y = "Frequency") +
  theme_minimal() +
  theme(text = element_text(size=20))

hist2 <- ggplot(final, aes(nr_indiv)) +
  geom_histogram(bins = 15) +
  labs(x = "Individuals per species",
       y = "Frequency") +
  theme_minimal() +
  theme(text = element_text(size=20))

# ggsave("site_year_histogram.png", hist1)
# ggsave("individuals_histogram.png", hist2)


