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


#----------------------------------------------------------------------
#----------------------------------------------------------------------
# # Save as rds
# saveRDS(data, file = "data/jungle_rhythms_data_cleaned.rds")
#----------------------------------------------------------------------
#----------------------------------------------------------------------
