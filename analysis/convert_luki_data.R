library(tidyverse)

# process Luki data into compatible format
df <- readxl::read_xlsx("data/Phenology Luki_1947-1958.xlsx")

# create date
week <- rep(1:36,12)
year <- sort(rep(1947:1958,36))
dates <- paste(year, week, sep = "-")
names(df)[10:ncol(df)] <- dates

# rearrange
df <- df %>% gather(key="dates",
                      value = "value",
                      -c(1:9)) %>%
  mutate(year = as.numeric(unlist(lapply(strsplit(dates, "-"),"[[",1))),
         week = as.numeric(unlist(lapply(strsplit(dates, "-"),"[[",2))),
         value = as.numeric(value),
         genus = unlist(lapply(strsplit(SpName," "),
                               function(x){if(length(x)==1){NA}else{x[1]}})),
         species = unlist(lapply(strsplit(SpName," "),
                                 function(x){if(length(x)==1){NA}else{x[2]}}))
         ) %>%
  rename(id = "Tree",
                phenophase = "X__3",
                notebook = "Notebook") %>%
  select(-c(
    "SpCode",
    "starting year",
    "X__1",
    "X__2",
    "dates",
    "Vernacular name",
    "SpName"
  )) %>%
  mutate(value = replace_na(value, 0),
         species_full = paste(genus, species, sep = " "))

# recode labels
df <- df %>% mutate(phenophase = recode(phenophase,FL = "flowering"),
              phenophase = recode(phenophase,FR = "fruiting"),
              phenophase = recode(phenophase,DISS = "fruit_drop"),
              phenophase = recode(phenophase,HIV = "leaf_dormancy"),
              phenophase = recode(phenophase,DEF = "leaf_turnover"))

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
