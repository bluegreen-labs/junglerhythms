# join data for weighted mean
library(tidyverse)

data <- readRDS("./data/jungle_rhythms_weekly_annotations.rds")
data$species_full <- paste(data$genus, data$species)

sp_list <- read.table("data/yangambi_mixed_forest_species_list.csv",
                      header = TRUE,
                      sep = ",")

sp_list <- sp_list %>% rename("species_full" = Species)

data <- inner_join(data, sp_list, by = "species_full")

data_w <- data %>%
  filter(phenophase == "flowers") %>%
  na.omit()

data_w <- data_w %>%
  group_by(species_full, week) %>%
  summarise(mean_week = mean(value, na.rm = TRUE) * unique(BAperc))

final <- data_w %>%
  group_by(week) %>%
  summarise(ss = sum(mean_week, na.rm = TRUE))

plot(final$ss)
