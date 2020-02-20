library(tidyverse)

# process Luki data into compatible format
df <- readxl::read_xlsx("data-raw/yangambi_validation_data.xlsx")

# create date
week <- rep(1:48,22)
year <- sort(rep(1937:1958,48))
dates <- paste(year, week, sep = "-")
colnames(df)[7:ncol(df)] <- dates
colnames(df)[1:6] <- df[2,1:6]
df <- df[-c(1:2),]

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
tmp <- tmp %>% mutate(phenophase = recode(phenophase,Flower = "flowering"),
              phenophase = recode(phenophase,Fruit = "fruiting"),
              phenophase = recode(phenophase,Dissemination = "fruit_drop"),
              phenophase = recode(phenophase,dissemination = "fruit_drop"),
              phenophase = recode(phenophase,Defoliation = "leaf_dormancy"))

data <- df %>%
  group_by(species_full, year) %>%
  mutate(selection = if(all(value == 0)){
   NA
  }else{
    value
  }) %>%
  na.omit() %>%
  ungroup()

# save the new tidy data
saveRDS(df, "data/luki_weekly_annotations.rds")
