#' Create linear plot of phenology data
#'
#' @param data junglerhythms data file
#' @param species_name end coordinate fo a line
#' @export
#' @return ggplot object

library(tidyverse)
library(viridis)


circular_linear_plot <- function(
  data,
  species_name = "Millettia laurentii",
  drought_mm = 140,
  leg_pos = c(1,0.1),
  leg_gradient = c(0,0.2,1),
  title_name = "(a) evergreen"
){

  #------------------------------------------------------------------------
  # data for circular plots
  #------------------------------------------------------------------------

  # split out a particular species
  data_subset_circ <- data %>%
    filter(grepl(tolower(species_name),
                 tolower(species_full)))

  # group by week and take the mean value
  data_subset_circ <- data_subset_circ %>%
    group_by(species_full, week, phenophase) %>%
    dplyr::summarise(mean_value = mean(value, na.rm=TRUE))#,

  # add a month column
  data_subset_circ$month <- months(as.Date(paste("1952",
                                                 round((data_subset_circ$week*7.625)-7),sep="-"), "%Y-%j"))

  # locate positions on a y-axis
  data_subset_circ$pos <- NA
  data_subset_circ$pos[data_subset_circ$phenophase == "flowers"] <- 2
  data_subset_circ$pos[data_subset_circ$phenophase == "fruit"] <- 1.5
  data_subset_circ$pos[data_subset_circ$phenophase == "fruit_drop"] <- 1
  data_subset_circ <- na.omit(data_subset_circ)



  #------------------------------------------------------------------------
  # linear flowering
  #------------------------------------------------------------------------
  # subset on species
  data_subset_lin_flower <- data %>%
    filter(grepl(tolower(species_name), tolower(species_full)),
           phenophase == "flowers")

  species_nr <- length(unique(data_subset_lin_flower$species_full))
  if (species_nr == 1){
    species_nr = 4
  }

  # convert date
  data_subset_lin_flower$date <- as.Date(
    paste(data_subset_lin_flower$year,
          round((data_subset_lin_flower$week*7.6)-7),sep="-"), "%Y-%j")

  # average by week
  data_subset_lin_flower <- data_subset_lin_flower %>%
    group_by(species_full, date) %>%
    dplyr::summarise(mean_value = mean(value),
                     scaled_value = ifelse(any(value > 0), 1, 0))

  #------------------------------------------------------------------------
  # linear fruiting
  #------------------------------------------------------------------------
  # subset on species
  data_subset_lin_fruit <- data %>%
    filter(grepl(tolower(species_name), tolower(species_full)),
           phenophase == "fruit")
  species_nr2 <- length(unique(data_subset_lin_fruit$species_full))
  if (species_nr == 1){
    species_nr = 4
  }
  # convert date
  data_subset_lin_fruit$date <- as.Date(
    paste(data_subset_lin_fruit$year,
          round((data_subset_lin_fruit$week*7.6)-7),sep="-"), "%Y-%j")
  # average by week
  data_subset_lin_fruit <- data_subset_lin_fruit %>%
    group_by(species_full, date) %>%
    dplyr::summarise(mean_value = mean(value),
                     scaled_value = ifelse(any(value > 0), 1, 0))

  #------------------------------------------------------------------------
  # linear fruit drop
  #------------------------------------------------------------------------
  # subset on species
  data_subset_lin_dispersal <- data %>%
    filter(grepl(tolower(species_name), tolower(species_full)),
           phenophase == "fruit_drop")
  species_nr2 <- length(unique(data_subset_lin_dispersal$species_full))
  if (species_nr == 1){
    species_nr = 4
  }
  # convert date
  data_subset_lin_dispersal$date <- as.Date(
    paste(data_subset_lin_dispersal$year,
          round((data_subset_lin_dispersal$week*7.6)-7),sep="-"), "%Y-%j")
  # average by week
  data_subset_lin_dispersal <- data_subset_lin_dispersal %>%
    group_by(species_full, date) %>%
    dplyr::summarise(mean_value = mean(value),
                     scaled_value = ifelse(any(value > 0), 1, 0))

  #------------------------------------------------------------------------
  # circular plot
  #------------------------------------------------------------------------
  p_circ <- ggplot(data = data_subset_circ,
                   aes(
                     x = week,
                     xend = week + 1,
                     y = pos,
                     yend = pos,
                     colour = mean_value
                   )) +
    scale_colour_distiller(palette = "YlOrBr",
                           direction = 1,
                           name = "Freq.",
                           values = leg_gradient) + #c(0,0.2,1)
    annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 2.4, alpha = .2) + #jan - feb
    annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 2.4, alpha = .2) + # jun-jul
    annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 2.4, alpha = .2) + # dec
  geom_segment(size = 4) +
    scale_x_continuous(limits = c(1,49),
                       breaks = seq(1,48,4),
                       labels = month.abb) +
    scale_y_continuous(limits = c(0,2.4)) +
    coord_polar() +
    labs(x="",
         y="",
         title = "") +
    theme(panel.grid.major.x = element_line(colour = "grey89",
                                            size = 0.3),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.background = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 11),
          strip.text = element_blank(),
          strip.background = element_rect(fill="white"),
          legend.position = leg_pos,
          legend.box.margin=margin(c(50,50,50,150)),
          plot.margin = unit(c(0,3,0.7,0.5),"cm")
    ) +
    facet_wrap(~ species_full,ncol=1)

  #------------------------------------------------------------------------
  # linear plot
  #------------------------------------------------------------------------
  p_lin <- ggplot(data_subset_lin_flower, aes(x = date,
                                          y = mean_value)) +
    geom_line(stat = "identity",aes(color="1Flowering")) +
    geom_line(data=data_subset_lin_fruit,aes(color="2Fruiting")) +
    geom_line(data=data_subset_lin_dispersal,aes(color="3Dispersal")) +
    scale_colour_manual(values = c("blue","red","orange")) +
    theme_minimal() +
    labs(title = title_name,
         y = "Freq. phenological event",
         x = "Year",
         color = "Event") +
    scale_x_date(date_breaks = "1 years",
                 date_labels = "%Y",
                 limits = as.Date(c('1937-01-01','1956-12-31'))) +
    scale_y_continuous(limits = c(0,1),
                       breaks = c(0,0.5,1)) +
    theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = 'white', colour = 'white'),
          plot.title = element_text(hjust = -0.25),
          strip.text = element_text(hjust = 0, size = 13),
          axis.line.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.ticks.x = element_blank(),
          legend.position="none",
          plot.margin = unit(c(0.1,0,1,1),"cm")
    ) +
    facet_wrap( ~ species_full, ncol = 1)


  return(grid.arrange(p_lin, p_circ, widths = c(2,1.5)))
}


# p_repro <- circular_linear_plot(data,
#                                 species_name = "Entandrophragma angolense",
#                                 title_name = "")
