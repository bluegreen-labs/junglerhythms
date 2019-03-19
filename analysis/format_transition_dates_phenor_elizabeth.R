# extract transition dates from JR
# weekly annotations to be ingested by the
# format_csv() function of phenor
library(zoo)
library(tidyverse)
library(circular)

#----------------------------------------------------------------------
#--------   Phenology data - species correction Meise   ---------------
#----------------------------------------------------------------------
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
df <- df[which(df$value != 0),]
df$join_id <- paste0("R",df$image,"-",df$image_row)

metadata <- read.csv("data/phenology_archives_species_long_format_20190319.csv",
                     header = TRUE, sep = ",")
metadata$join_id <- paste(metadata$image,metadata$row, sep = "-")

# test merge the two tables based upon the unique join ID
data <- merge(df, metadata, by = c("join_id","id"), all.x = TRUE)
data$species_full <- paste(data$genus_Meise, data$species_Meise)
#----------------------------------------------------------------------

# read in census data
census <- read.csv2("data/yangambi_mixed_forest_species_list.csv",
                    header = TRUE,
                    sep = ",",
                    stringsAsFactors = FALSE)
colnames(census)[1] <- 'species_full'


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------- This is onset of event    --------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
transition_dates <- data %>%
  filter(phenophase == "leaf_turnover" ) %>% #| phenophase == "leaf_turnover"
  group_by(species_full, id, phenophase) %>%
  do(event_length(.))

# transition_dates <- transition_dates %>%
#   mutate(date = as.Date(sprintf("%s-%02d-1",year_start, week_start),
#                         "%Y-%W-%u"))  # new column 'date' with year-month-day

# then reformat those taking the year - date line into
# consideration
transition_dates <- transition_dates %>%
  mutate(degree = week_start * 360/48)

bla <- transition_dates %>%
  filter(length(week_start) > 5) %>%    # only species with more than 10 events
  summarise(
    mean_degree = mean.circular(
      circular(degree, units = "degrees")
    ),
    lower = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees")
    )$mu.ci[1],
    upper = mle.vonmises.bootstrap.ci(
      circular(degree, units = "degrees")
    )$mu.ci[2],
    nr_events = length(week_start)
  )

# merge with census and remove rows (species) not included in the Yangambi mixed forest census
bla <- merge(bla, census, by = "species_full", all.x = TRUE)
bla <- bla[!(is.na(bla$BAperc)),]

# bla$lower <- ifelse(bla$nr_event == 1, NA, bla$lower)
# bla$upper <- ifelse(bla$nr_event == 1, NA, bla$upper)

# sort by BA
# number by sorted BA
# set y value by sorted BA number
bla <- bla %>%
  arrange(BA) %>%
  mutate(y_value = (1:length(species_full)))

# rescaling between 0 and 360 degrees
bla$mean_rescaled <- ifelse(bla$mean_degree < 0, bla$mean_degree +360, bla$mean_degree)
bla$upper_rescaled <- ifelse (bla$upper < 0, bla$upper + 360, bla$upper)
bla$upper_rescaled <- circular(bla$upper_rescaled, units = "degrees")
bla$lower_rescaled <- ifelse (bla$lower < 0, bla$lower + 360, bla$lower)
bla$lower_rescaled <- circular(bla$lower_rescaled, units = "degrees")


bla_subset <- bla %>%
  filter(mean_rescaled < lower_rescaled | mean_rescaled > upper_rescaled)
bla_subset$seg1_x <- 0
bla_subset$seg2_x <- 360

bla_ok <- bla %>%
  filter(mean_rescaled > lower_rescaled & mean_rescaled < upper_rescaled)

# bla_single <- bla %>%
#   filter(nr_events == 1)




# set x end points by upper and lower limits
# plot mean as point

ggplot(data = bla_ok) +
  geom_point(aes(x = mean_rescaled, y = y_value)) +
  geom_segment(aes(x = lower_rescaled,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value)) +
  geom_point(data=bla_subset, aes(x = mean_rescaled, y = y_value)) +
  geom_segment(data=bla_subset, aes(x = upper_rescaled,
                                    xend = seg1_x,
                                    y = y_value,
                                    yend = y_value)) +
  geom_segment(data=bla_subset, aes(x = lower_rescaled,
                                    xend = seg2_x,
                                    y = y_value,
                                    yend = y_value)) +
  # geom_point(data=bla_single, aes(x = mean_rescaled, y = y_value)) +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0,359,30),
                     labels = month.abb) +
  scale_y_continuous(limits = c(-20,length(bla_ok$mean_rescaled) )) + #staring at '-10' so that points are squeezed in the center
  annotate("rect", xmin = 330, xmax = 360, ymin = 0, ymax = length(bla_ok$mean_rescaled) , alpha = .2) + #Dec
  annotate("rect", xmin = 0, xmax = 60, ymin = 0, ymax = length(bla_ok$mean_rescaled) , alpha = .2) + # jan - feb
  annotate("rect", xmin = 150, xmax = 210, ymin = 0, ymax = length(bla_ok$mean_rescaled) , alpha = .2) + # jun-jul
  coord_polar() +
  labs(x="",
       y="",
       title = "Leaf turnover") +
  theme(panel.grid.major.x = element_line(colour = "grey75",
                                          size = 0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11),
        strip.text = element_text(face = "italic", size = 13),
        strip.background = element_rect(fill="white"),
        # legend.position = "right",
        plot.margin=unit(c(1,0,0,0),"cm")
  )

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------- This is mean event time, not onset of event    -----------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
test <- data %>%
  filter(phenophase == "leaf_turnover",
         value == 1) %>%
  # filter(species_full == "Prioria balsamifera") %>%
  mutate(degree = week * 360/48)


bla <- test %>%
  group_by(species_full) %>%
  filter(length(value) > 10) %>%     # we only look at species with more than 10 events
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

# merge with census and remove rows (species) not included in the Yangambi mixed forest census
bla <- merge(bla, census, by = "species_full", all.x = TRUE)
bla <- bla[!(is.na(bla$BAperc)),]

# sort by BA
# number by sorted BA
# set y value by sorted BA number
bla <- bla %>%
  arrange(BA) %>%
  mutate(y_value = 1:length(species_full))

# rescaling between 0 and 360 degrees
bla$mean_rescaled <- ifelse(bla$mean_degree < 0, bla$mean_degree +360, bla$mean_degree)
bla$upper_rescaled <- ifelse (bla$upper < 0, bla$upper + 360, bla$upper)
bla$upper_rescaled <- circular(bla$upper_rescaled, units = "degrees")
bla$lower_rescaled <- ifelse (bla$lower < 0, bla$lower + 360, bla$lower)
bla$lower_rescaled <- circular(bla$lower_rescaled, units = "degrees")

bla_subset <- bla %>%
  filter(mean_rescaled < lower_rescaled | mean_rescaled > upper_rescaled)
bla_subset$seg1_x <- 0
bla_subset$seg2_x <- 360

bla_ok <- bla %>%
  filter(mean_rescaled > lower_rescaled & mean_rescaled < upper_rescaled)


# set x end points by upper and lower limits
# plot mean as point

ggplot(data = bla_ok) +
  geom_point(aes(x = mean_rescaled, y = y_value)) +
  geom_segment(aes(x = lower_rescaled,
                   xend = upper_rescaled,
                   y = y_value,
                   yend = y_value)) +
  geom_point(data=bla_subset, aes(x = mean_rescaled, y = y_value)) +
  geom_segment(data=bla_subset, aes(x = upper_rescaled,
                   xend = seg1_x,
                   y = y_value,
                   yend = y_value)) +
  geom_segment(data=bla_subset, aes(x = lower_rescaled,
                                    xend = seg2_x,
                                    y = y_value,
                                    yend = y_value)) +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0,359,30),
                     labels = month.abb) +
  scale_y_continuous(limits = c(-20,length(bla$mean_rescaled))) + #staring at '-10' so that points are squeezed in the center

  annotate("rect", xmin = 330, xmax = 360, ymin = 0, ymax = length(bla$mean_rescaled), alpha = .2) + #Dec
  annotate("rect", xmin = 0, xmax = 60, ymin = 0, ymax = length(bla$mean_rescaled), alpha = .2) + # jan - feb
  annotate("rect", xmin = 150, xmax = 210, ymin = 0, ymax = length(bla$mean_rescaled), alpha = .2) + # jun-jul
  coord_polar() +
  labs(x="",
       y="",
       title = "Leaf turnover") +
  theme(panel.grid.major.x = element_line(colour = "grey75",
                                          size = 0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11),
        strip.text = element_text(face = "italic", size = 13),
        strip.background = element_rect(fill="white"),
        # legend.position = "right",
        plot.margin=unit(c(1,0,0,0),"cm")
  )
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
