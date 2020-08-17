#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
#----- load required libraries -------------------------------------------------#
library(tidyverse)
#-------------------------------------------------------------------------------#

#----------------------------------------------------------------------
#--------   Phenology Yangambi - species correction Meise   -----------
#----------------------------------------------------------------------
metadata <- read.csv("data/phenology_archives_species_long_format_20200324.csv",
                     header = TRUE, sep = ",")

metadata$species_full <- paste(metadata$genus_Meise, metadata$species_Meise)

# only select parameters you need, more clear structure to work with
metadata <- metadata %>%
  select(species_full,
         family_Meise,
         genus_Meise,
         species_Meise,
         subsp_Meise,
         author_Meise,
         status_Meise)

unique(metadata$status_Meise)

ygb_sp <- metadata[!duplicated(metadata), ]
ygb_sp <- ygb_sp %>%
  filter(status_Meise != c("Issues"))

#----------------------------------------------------------------------
#--------   Phenology Luki --------------------------------------------
#----------------------------------------------------------------------
luki <- read.csv("data/Phenology_Luki_species.csv",
                     header = TRUE, sep = ",")

luki <- luki %>%
  rename(species_full = Species_name)


luki <- merge(luki, ygb_sp, by = c("species_full"), all.x = TRUE)

# write.table(luki, "data/luki_corrections.csv",
#             quote = FALSE,
#             col.names = TRUE,
#             row.names = FALSE,
#             sep = ",")

#----------------------------------------------------------------------
#--------   Phenology data - Luki -------------------------------------
# actual species names used in the rds file
#----------------------------------------------------------------------
df <- readRDS("data/luki_weekly_annotations.rds")
metadata <- read.csv("data/luki_species_corrections.csv",
                     header = TRUE, sep = ",")
df <- df %>%
  rename(luki_original = species_full)
df <- merge(df, metadata, by = c("luki_original"), all.x = TRUE)

species_list <- as.data.frame(unique(df$luki_original))
colnames(species_list)[1] <- "luki_original"
species_list <- merge(species_list, metadata, by = c("luki_original"), all.x = TRUE)

# write.table(species_list, "data/luki_species_corrections_v2.csv",
#             quote = FALSE,
#             col.names = TRUE,
#             row.names = FALSE,
#             sep = ",")

