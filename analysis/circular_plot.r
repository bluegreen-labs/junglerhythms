#' Create circular plot of phenology data
#'
#' @param data junglerhythms data file
#' @param species_name end coordinate fo a line
#' @param threshold threshold to use when summarizing weekly data
#' @export
#' @return ggplot object

circle_plot <- function(data,
                        species_name = "Millettia laurentii",
                        threshold = 0.5){

  rescale <- function(x){
    if(all(x == 0)){
      return(x)
    } else {
      scales::rescale(x)
    }
  }

  # juggling with data
  data$species_full <- paste(data$genus, data$species)

  # split out a particular species
  data_subset <- data %>% filter(grepl(tolower(species_name),
                             tolower(species_full)))

  # group by week and take the mean value
  data_subset <- data_subset %>% group_by(week, phenophase) %>%
    summarise(mean_value = mean(value, na.rm=TRUE))

  # rescale and threshold (majority vote)
  data_subset <- data_subset %>% group_by(phenophase) %>%
    mutate(mean_value = rescale(mean_value)) %>%
    mutate(week = ifelse(mean_value <= threshold, NA, week))

  # add a month column
  data_subset$month <- months(as.Date(paste("1952",
                      round(data_subset$week*7.625),sep="-"), "%Y-%j"))

  # locate positions on a y-axis
  data_subset$pos <- NA
  data_subset$pos[data_subset$phenophase == "senescence"] <- 4
  data_subset$pos[data_subset$phenophase == "flowers"] <- 3
  data_subset$pos[data_subset$phenophase == "fruit"] <- 2
  data_subset$pos[data_subset$phenophase == "fruit_drop"] <- 1

  # create plot and return ggplot plot object
  p <- ggplot(data = data_subset,
               aes(
                 x = week,
                 xend = week + 1,
                 y = pos,
                 yend = pos,
                 group = phenophase,
                 colour = phenophase
               )) +
    geom_segment(size = 10) +
    scale_colour_ptol(labels = c("flowers",
                                 "fruit",
                                 "fruit drop",
                                 "senescence"),
                      name = "") +
    scale_x_continuous(limits = c(1,49),
                       breaks = seq(1,48,4),
                       labels = month.abb) +
    scale_y_continuous(limits = c(0,4.4)) +
    coord_polar() +
    labs(x="",
         y="",
         title = "") +
    theme(panel.grid.major.x = element_line(colour = "grey89", size=0.3),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.background = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 20),
          legend.position="bottom",
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          plot.margin=unit(c(1,1,1,1),"cm")
    )
  return(p)
}
