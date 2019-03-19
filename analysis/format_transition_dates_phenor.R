# extract transition dates from JR
# weekly annotations to be ingested by the
# format_csv() function of phenor
library(zoo)
library(tidyverse)
library(circular)

data <- readRDS("./data/jungle_rhythms_weekly_annotations.rds")
data$species_full <- paste(data$genus, data$species)

# # first, find the date on which an event happens
# transition_dates <- data %>%
#   filter(phenophase == "leaf_dormancy") %>%
#   mutate(date = as.Date(sprintf("%s-%02d-1",year, week),
#                                                 "%Y-%W-%u")) %>%
#   group_by(id, year, genus, species) %>%
#   arrange(date) %>%
#   mutate(value = ifelse(value == 0, NA, 1)) %>%
#   mutate(value = na.approx(value, maxgap = 2, na.rm = FALSE)) %>%
#   mutate(value = ifelse(is.na(value), 0, 1)) %>%
#   mutate(diff_value = ifelse(c(value - lag(value)) == 1,
#                              1, 0)) %>%
#   summarize(week = week[which(diff_value == 1)[1]]) %>%
#   na.omit() %>%
#   ungroup()
#
# transition_dates <- data %>%
#   filter(phenophase == "leaf_dormancy") %>%
#   group_by(species_full, week) %>%


# then reformat those taking the year - date line into
# consideration
transition_dates <- transition_dates %>%
  mutate(degree = week * 360/48)

# approximate location of Yangambi
transition_dates$lat <- 0.804593
transition_dates$lon <- 24.452605

transition_dates$species_full <- paste(transition_dates$genus,
                                       transition_dates$species)

test <- data %>%
  filter(phenophase == "leaf_dormancy",
         value == 1) %>%
  #filter(species_full == "Samanea dinklagei") %>%
  mutate(degree = week * 360/48)


bla <- test %>%
  group_by(species_full) %>%
  filter(length(value) > 10) %>%
  summarise(
    mean_degree = mean.circular(
      circular(degree, units = "degrees")
    ),
    lower = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees")
      )$mu.ci[1],
    upper = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees")
    )$mu.ci[2]
    )

x <- circular(test$degree, units = "degrees")
print(mean(x))
ci <- mle.vonmises.bootstrap.ci(x)
print(ci)
rose.diag(ci$mu, bins=30, main=expression(mu))

bla <- bla %>%
  arrange(mean_degree) %>%
  mutate(y_value = 1:length(species_full))

# sort by dbh
# number by sorted dbh
# set y value by sorted dbh number
# set x end points by upper and lower limits
# plot mean as point
ggplot(data = bla) +
  geom_point(aes(x = mean_degree, y = y_value)) +
  geom_segment(aes(x = upper,
                   xend = lower,
                   y = y_value,
                   yend = y_value))

