library(tidyverse)
library(cowplot)

# process Luki data into compatible format
if(!exists("df")){
  df <- readxl::read_xlsx("data-raw/yangambi_validation_data.xlsx")

  week <- rep(1:48,22)
  year <- sort(rep(1937:1958,48))
  dates <- paste(year, week, sep = "-")
  colnames(df)[7:ncol(df)] <- dates
  colnames(df)[1:6] <- df[2,1:6]
  df <- df[-c(1:2),]
}

# rearrange
tmp <- df %>% gather(key="dates",
                      value = "value",
                      -c(1:6)) %>%
  mutate(year = as.numeric(unlist(lapply(strsplit(dates, "-"),"[[",1))),
         week = as.numeric(unlist(lapply(strsplit(dates, "-"),"[[",2))),
         value = as.numeric(value),
         image = gsub("R","",Image),
         genus = unlist(lapply(strsplit(Species," "),
                               function(x){if(length(x)==1){NA}else{x[1]}})),
         species = unlist(lapply(strsplit(Species," "),
                                 function(x){if(length(x)==1){NA}else{x[2]}}))
         ) %>%
  rename(id = "ID",
         phenophase = "Observation") %>%
  select(-c(
    "start year",
    "End year",
    "Species",
    "dates",
    "Image"
  )) %>%
  mutate(
    value = replace_na(value, 0),
    species_full = paste(genus, species, sep = " ")
    )

# inspect data
#tmp %>% filter(phenophase == "Fruit") %>% select(value) %>% table()

# recode labels
# What about the different numbers in the coding 2 and 4 in the
# leaf dormancy section (probably leaf_dormancy / leaf_turnover differences)
# check with original data

tmp <- tmp %>% mutate(phenophase = recode(phenophase,Flower = "flowers"),
              phenophase = recode(phenophase,Fruit = "fruit"),
              phenophase = recode(phenophase,Dissemination = "fruit_drop"),
              phenophase = recode(phenophase,dissemination = "fruit_drop"),
              phenophase = recode(phenophase,Defoliation = "leaf_dormancy"))

# data retention strategy, very strict, needs to align with
# other processing to make validation valid.
tmp <- tmp %>%
  group_by(species_full, id, year) %>%
  mutate(selection = if(all(value == 0)){
   NA
  }else{
    value
  }) %>%
  ungroup()  %>%
  na.omit()

# save the new tidy data
saveRDS(tmp, "data/yangambi_validation_annotations.rds")
