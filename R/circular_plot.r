#o <- readRDS("data/jungle_rhythms_weekly_annotations.rds")

circle_plot <- function(data,
                        species_name = "Millettia laurentii",
                        threshold = 0.5){

  data$species_full <- paste(data$genus, data$species)

k <- data %>% filter(grepl(tolower(species_name),
               tolower(species_full))) %>%
    filter(phenophase != "fruit_drop")

k <- k %>% group_by(week, phenophase) %>%
  summarise(mean_value = mean(value, na.rm=TRUE))

k <- k %>% group_by(phenophase) %>%
  mutate(mean_value = scales::rescale(mean_value))

k <- k %>%
  mutate(mean_value = ifelse(mean_value < threshold,NA,1)) %>%
  na.omit()

k$month <- months(as.Date(paste("1952",
                                round(k$week*7.625),sep="-"), "%Y-%j"))

# locate positions on a y-axis
k$pos <- NA
k$pos[k$phenophase == "senescence"] <- 1
k$pos[k$phenophase == "flowers"] <- 2
k$pos[k$phenophase == "fruit"] <- 3

print(unique(k$month))

b2 <- ggplot(data = k,
             aes(
               x = week,
               xend = week + 1,
               y = pos,
               yend = pos,
               group = phenophase,
               colour = phenophase
             )) +
  geom_segment(size = 10) +
  scale_colour_economist() +
  scale_x_continuous(limits = c(1,49),
                     breaks = seq(1,48,4),
                     labels = month.name) +
  scale_y_continuous(limits = c(0,3.2)) +
  coord_polar() +
  labs(x="",
       y="",
       title = species_name) +
  theme(panel.grid.major.x = element_line(colour = "grey89", size=0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(family = "Times",
                                   colour = "black",
                                   size = 12),
        legend.position="right",
        plot.margin=unit(c(1,1,1,1),"cm")
  )

plot(b2)

}

circle_plot(data = weekly_data, threshold = 0.5)
