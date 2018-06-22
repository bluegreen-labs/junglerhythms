# plot things

library(tidyverse)
library(ggthemes)

#k <- readRDS("data/jungle_rhythms_weekly_annotations.rds")

my_species <- "Millettia laurentii"

d = k %>% filter(grepl(tolower(my_species),
                       tolower(paste(genus, species)))) %>%
  mutate(value=ifelse(value == 0, NA, 1)) %>%
  na.omit()

d$date <- as.Date(paste(d$year,round(d$week*7.6),sep="-"), "%Y-%j")

l <- length(unique(d$id))

species <- unique(d$species)

loc <- which(d$phenophase == "flowers")
d$value[loc] <- d$value[loc] * 2.5

loc <- which(d$phenophase == "fruit")
d$value[loc] <- d$value[loc] * 2

loc <- which(d$phenophase == "fruit_drop")
d$value[loc] <- d$value[loc] * 1.5

loc <- which(d$phenophase == "senescence")
d$value[loc] <- d$value[loc] * 1

p = ggplot(d, aes(x = date,
             y = value,
             colour = phenophase,
             group = phenophase)) +
  geom_point(stat = "identity") +
  theme_minimal() +
  scale_y_continuous(breaks = NULL) +
  labs(title = species, x = "Year") +
  facet_wrap( ~ id, nrow = l, strip.position = "left")
plot(p)
