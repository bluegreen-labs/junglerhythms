# compile summary statistics
#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
# graphics.off()
#----- load required packages --------------------------------------------------#
library(tidyverse)
library(ggplot2)
library(ggthemes)
#-------------------------------------------------------------------------------#

# read in the weekly data
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
df <- df[which(df$value != 0),]

df$join_id <- paste0("R",df$image,"-",df$image_row)

metadata <- read.csv("data/phenology_archives_species_long_format_20190626.csv",
           header = TRUE, sep = ",")
metadata$join_id <- paste(metadata$image,metadata$row, sep = "-")

# test merge the two tables based upon the unique join ID
test <- merge(df, metadata, by = "join_id", all.x = TRUE)

# read in census data and convert basal area
census <- read.csv2("data/YGB_ForestPlotsNET_corrected_indet.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)

# calculate mean basal area across all species across all
# mixed plots
census <- census %>%
  filter(grepl("MIX",Plot)) %>%
  group_by(Plot, Species) %>%
  summarize(basal_area = sum(pi*(C2DBH4/2000)^2, na.rm = TRUE)) %>%
  group_by(Species) %>%
  summarize(basal_area = mean(basal_area))

# replace empty slots with Unknown for consistency
df$genus <- ifelse(df$genus == "","Unknown",df$genus)
df$species <- ifelse(df$species == "","Unknown",df$species)

# summarize site years and nr. individuals, start and end years
test <- df %>% group_by(genus, species) %>%
  summarize(start = min(year, na.rm = TRUE),
            end = max(year, na.rm = TRUE),
            site_years = length(unique(year)),
            individuals = length(unique(id)),
            Species = unique(paste(genus, species))) %>%
  filter(genus != "Unknown") %>%
  na.omit()

# merge and drop merge column
test <- merge(test, census, by = "Species", all.x = TRUE)
test <- test %>% select(-Species)

# write to file
write.table(test, "data/species_meta_data.csv",
          quote = FALSE,
          col.names = TRUE,
          row.names = FALSE,
          sep = ",")

# calculate the percentage of trees with a basal area measurement
pct_trees_BA <- length(which(!is.na(test$basal_area)))/nrow(test)
message(pct_trees_BA)

# sum the basal area of the observed species
# if >30 good coverage
message(sum(test$basal_area, na.rm = TRUE))

hist1 <- ggplot(test, aes(site_years)) +
  geom_histogram(bins = 10) +
  labs(x = "Site Years",
       y = "Frequency") +
  theme_minimal() +
  theme(text = element_text(size=20))

hist2 <- ggplot(test, aes(individuals)) +
  geom_histogram(bins = 10) +
  labs(x = "Individuals per species",
       y = "Frequency") +
  theme_minimal() +
  theme(text = element_text(size=20))

ggsave("site_year_histogram.png", hist1)
ggsave("individuals_histogram.png", hist2)


