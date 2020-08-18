#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
# graphics.off()
#----- load required packages --------------------------------------------------#
library(tidyverse)

#----------------------------------------------------------------------
#--------   Phenology data - species correction Meise   ---------------
#----------------------------------------------------------------------
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")

metadata <- read.csv("data/phenology_archives_species_long_format_20200324.csv",
                     header = TRUE, sep = ",")


# merge the two tables based upon the unique join ID
df$join_id <- paste0("R",df$image,"-",df$image_row)
metadata$join_id <- paste(metadata$image,metadata$row, sep = "-")
df <- df %>%
  select(!id) # remove id as updated in metadata (in metadata, empty ids are renamed to EK1, EK2, etc...))

data <- merge(df, metadata, by = c("join_id"), all.x = TRUE)

# get last updated species name
data$species_full <- paste(data$genus_Meise, data$species_Meise)

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
  select(image.y,
         image_row,
         image_col,
         species_full,
         id,
         phenophase,
         year,
         week,
         value)
data <- data %>%
  rename(image = image.y)
#----------------------------------------------------------------------
rm(df,metadata, empty_years)
#----------------------------------------------------------------------


# # ----------------------------------------------------------------------
# # ----------------------------------------------------------------------
# # Save as rds
# saveRDS(data, file = "data/jungle_rhythms_data_cleaned.rds")
# # ----------------------------------------------------------------------
# # ----------------------------------------------------------------------

#----------------------------------------------------------------------
#-------- census data  ------------------------------------------------
#-------- for manuscripts on leaf phenology and reproductive phenology
#-------- only work with species in 2012 Yangambi census
#----------------------------------------------------------------------
census <- read.csv("data/YGB_ForestPlotsNET_corrected_indet.csv",
                   header = TRUE,
                   sep = ",",
                   stringsAsFactors = FALSE)
census <- census %>%
  dplyr::rename("species_full" = Species)
census$species_full <- ifelse(census$species_full == "Unknown", NA, census$species_full)

# remove individuals without C1DBH4, these are new recruits from census2
# only use mixed plots
# remove those only at genus level
census_sp <- census %>%
  filter(!is.na(C1DBH4),
         grepl("MIX",Plot),
         !grepl("sp\\.",species_full))
# species lists census vs phenology datasets
census_species <- unique(census_sp$species_full)
phen_species <- unique(data$species_full)
species_list <- intersect(census_species, phen_species)
#----------------------------------------------------------------------
rm(census_species, phen_species, census_sp)
#----------------------------------------------------------------------
# cleaned dataset only for these species
data_sp <- data %>%
  filter(species_full %in% species_list)

# #----------------------------------------------------------------------
# #----------------------------------------------------------------------
# # Save as rds
# saveRDS(data_sp, file = "data/jungle_rhythms_data_manuscript_leaf_repro.rds")
# #----------------------------------------------------------------------
# #----------------------------------------------------------------------


