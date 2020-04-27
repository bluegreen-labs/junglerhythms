#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
# graphics.off()
#----- load required packages --------------------------------------------------#
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(BBmisc)
library(dtw)
library(dtwclust)
# library(cluster)
#-------------------------------------------------------------------------------#

#----------------------------------------------------------------------
#--------   Phenology data - species correction Meise   ---------------
#----------------------------------------------------------------------
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
df$join_id <- paste0("R",df$image,"-",df$image_row)

metadata <- read.csv("data/phenology_archives_species_long_format_20200324.csv",
                     header = TRUE, sep = ",")
metadata$join_id <- paste(metadata$image,metadata$row, sep = "-")
# merge the two tables based upon the unique join ID
data <- merge(df, metadata, by = c("join_id"), all.x = TRUE)
data$species_full <- paste(data$genus_Meise, data$species_Meise)

# remove column id.x and rename id.y to id (--> in id.y, empty ids are renamed to EK1, EK2, etc...)
data$id <- as.character(data$id.y)
data = data[,!(names(data) %in% c("id.x","id.y"))]

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
  select(species_full,
         id,
         phenophase,
         year,
         week,
         value)
#----------------------------------------------------------------------
rm(df,metadata, empty_years)
#----------------------------------------------------------------------

# #----------------------------------------------------------------------
# #-------- get the species list ----------------------------------------
# #----------------------------------------------------------------------
overview <- read.csv("data/species_meta_data.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)

# merge number of events to data, so you can filter on minimum number of events
overview <- overview[,(names(overview) %in% c("species_full",
                                              "deciduousness",
                                              "basal_area_site"))]

# merge with data to filter on deciduousness
data <- merge(data, overview, by = "species_full", all.x = TRUE)

# get species_list of species
# - in basal area list
species_list_dormancy <- data %>%
  filter(basal_area_site > 0) %>%
  filter(phenophase == "leaf_dormancy") %>%
  group_by(species_full) %>%
  dplyr::summarise(count_events = sum(value), # total events
                   total_week = length(value)) # total weeks
species_list_dormancy <- species_list_dormancy %>%
  filter(count_events >=5)
species_list_dormancy <- species_list_dormancy$species_full
#----------------------------------------------------------------------

# species_specific week-based annual means
data_subset <- data %>%
  filter(species_full %in% species_list_dormancy) %>%
  filter(phenophase == "leaf_dormancy") %>%
  filter(grepl("deciduous",deciduousness)) %>%
  group_by(species_full, week) %>%
  dplyr::summarise(mean_value = mean(value, na.rm=TRUE))

# prep data-format for clustering
data_subset <- as.data.frame(spread(data_subset, key = "week", value = "mean_value")) # as data.frame because tibbles cant't have rownames
rownames(data_subset) <- data_subset$species_full
data_subset <- data_subset[,-1]
# normalise the data
dormancy.norm <- BBmisc::normalize(data_subset, method="standardize")


# cluster validation
# from https://stats.stackexchange.com/questions/351314/r-how-to-find-cluster-centroids-with-tsclust
cfg <- compare_clusterings_configs(
  "h",
  k = 2L:10L,
  controls = list(
    hierarchical = hierarchical_control(method = "all", symmetric = TRUE) # method = "average"
  ),
  distances = pdc_configs("d", dtw_basic = list()) #dtw_basic ; dtw2 ; gak ; lb_improved ; sdtw ; sbd
)

evaluator <- cvi_evaluators("CH")
comparison <- compare_clusterings(dormancy.norm, "h", cfg,
                                  score.clus = evaluator$score,
                                  pick.clus = evaluator$pick)

best <- repeat_clustering(dormancy.norm, comparison, comparison$pick$config_id)
plot(best)

sp_clust <- as.data.frame(best@cluster)
sp_clust$species_full <- rownames(sp_clust)
colnames(sp_clust)[1] <- "cluster_id"

# # hierarchical clustering using dynamic time warping (dtw) distance
# clust.hier <- dtwclust::tsclust(dormancy.norm, type = "h", k=6L, distance = "dtw")
# plot(clust.hier)

write.table(sp_clust, "data/cluster_test.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")


