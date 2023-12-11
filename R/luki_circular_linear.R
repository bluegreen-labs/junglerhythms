#' Create linear plot of phenology data
#'
#' @param data junglerhythms data file
#' @param species_name end coordinate fo a line
#' @export
#' @return ggplot object

library(tidyverse)
library(viridis)
library(gridExtra)
library(ggplot2)

luki_circular_linear_plot <- function(
  data,
  species_name = "Millettia laurentii",
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
    dplyr::summarise(percent_value = mean(value, na.rm=TRUE) *100)#,

  # add a month column
  data_subset_circ$month <- months(as.Date(paste("1952",
                                                 round((data_subset_circ$week*10.167)-9),sep="-"), "%Y-%j"))
  # locate positions on a y-axis
  data_subset_circ$pos <- NA
  data_subset_circ$pos[data_subset_circ$phenophase == "leaf_dormancy"] <- 2
  data_subset_circ$pos[data_subset_circ$phenophase == "leaf_turnover"] <- 1.5
  data_subset_circ <- na.omit(data_subset_circ)



  #------------------------------------------------------------------------
  # linear leaf dormancy
  #------------------------------------------------------------------------
  # subset on species
  data_subset_lin_LD <- data %>%
    filter(grepl(tolower(species_name), tolower(species_full)),
           phenophase == "leaf_dormancy")

  ids <- data_subset_lin_LD %>% # to add in facet_wrap titles (n = xx)
    group_by(species_full) %>%
    dplyr::summarise(ids = length(unique(id)))

  # average by date
  data_subset_lin_LD$date <- as.Date(paste(data_subset_lin_LD$year,
                                           round((data_subset_lin_LD$week*10.167)-9),sep="-"), "%Y-%j")
  data_subset_lin_LD <- data_subset_lin_LD %>%
    group_by(species_full, date) %>%
    dplyr::summarise(percent_value = mean(value, na.rm=TRUE) *100,
                     scaled_value = ifelse(any(value > 0), 1, 0))

  data_subset_lin_LD <- merge(data_subset_lin_LD, ids, by = c("species_full"), all.x = TRUE)
  data_subset_lin_LD$sp_title <- paste(data_subset_lin_LD$species_full, " (n = ", data_subset_lin_LD$ids, ")", sep = "")


  #------------------------------------------------------------------------
  # linear leaf turnover
  #------------------------------------------------------------------------
  # subset on species
  data_subset_lin_LT <- data %>%
    filter(grepl(tolower(species_name), tolower(species_full)),
           phenophase == "leaf_turnover")

  # average by date
  data_subset_lin_LT$date <- as.Date(paste(data_subset_lin_LT$year,
                                           round((data_subset_lin_LT$week*10.167)-9),sep="-"), "%Y-%j")
  data_subset_lin_LT <- data_subset_lin_LT %>%
    group_by(species_full, date) %>%
    dplyr::summarise(percent_value = mean(value, na.rm=TRUE) *100,
                     scaled_value = ifelse(any(value > 0), 1, 0))

  data_subset_lin_LT <- merge(data_subset_lin_LT, ids, by = c("species_full"), all.x = TRUE)
  data_subset_lin_LT$sp_title <- paste(data_subset_lin_LT$species_full, " (n = ", data_subset_lin_LT$ids, ")", sep = "")


  #------------------------------------------------------------------------
  # circular plot
  #------------------------------------------------------------------------
  p_circ <- ggplot(data = data_subset_circ,
                   aes(
                     x = week,
                     xend = week + 1,
                     y = pos,
                     yend = pos,
                     colour = percent_value
                   )) +
    scale_colour_distiller(palette = "YlOrBr",
                           direction = 1,
                           name = "mean annual \n% individuals\nwith events",
                           values = leg_gradient) +
    annotate("rect", xmin = 16, xmax = 28, ymin = 0, ymax = 2.4, alpha = .2) + #jun - sept

    # annotate("text", x = 1, y = 0.8, label = "LD", col = "grey50") +
    # annotate("text", x = 25, y = 0.8, label = "SD", col = "grey50") +
    # annotate("text", x = 37, y = 0.8, label = "LW", col = "grey50") +
    # annotate("text", x = 13, y = 0.8, label = "SW", col = "grey50") +

  geom_segment(size = 4) +
    scale_x_continuous(limits = c(1,37),
                       breaks = seq(1,36,3),
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
  p_lin <- ggplot(data_subset_lin_LD, aes(x = date,
                                          y = percent_value)) +
    geom_line(aes(color="Canopy dormancy")) +
    geom_line(data=data_subset_lin_LT,aes(color="Canopy turnover")) +
    scale_colour_manual(values = c("black","grey60")) +
    theme_minimal() +
    labs(title = title_name,
         y = "% of individuals with events",
         x = "Year",
         color = "Event") +
    scale_x_date(date_breaks = "1 years",
                 date_labels = "%Y",
                 limits = as.Date(c('1946-01-01','1958-12-31')),
                 expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,100),
                       breaks = c(0,50,100)) +
    theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = 'white', colour = 'white'),
          plot.title = element_text(hjust = -0.25),
          strip.text = element_text(hjust = 0, size = 13, face = "italic"),
          axis.line.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0.1,0,1,1),"cm")
    ) +
    facet_wrap( ~ sp_title, ncol = 1)

  return(grid.arrange(p_lin, p_circ, widths = c(2,1.5)))
}



# circular_plot <- function(
#   data,
#   species_name = "Millettia laurentii",
#   leg_pos = c(1,0.1),
#   leg_gradient = c(0,0.2,1),
#   title_name = "(a) evergreen"
# ){
#   #------------------------------------------------------------------------
#   # data for circular plots
#   #------------------------------------------------------------------------
#   # split out a particular species
#   data_subset_circ <- data %>%
#     filter(grepl(tolower(species_name),
#                  tolower(species_full)))
#
#   # group by week and take the mean value
#   data_subset_circ <- data_subset_circ %>%
#     group_by(species_full, week, phenophase, dummy) %>%
#     dplyr::summarise(percent_value = mean(value, na.rm=TRUE) *100)
#
#   # add a month column
#   data_subset_circ$month <- months(as.Date(paste("1952",
#                                                  round((data_subset_circ$week*7.625)-7),sep="-"), "%Y-%j"))
#   # locate positions on a y-axis
#   data_subset_circ$pos <- NA
#   data_subset_circ$pos[data_subset_circ$phenophase == "leaf_dormancy"] <- 2
#   data_subset_circ$pos[data_subset_circ$phenophase == "leaf_turnover"] <- 1.5
#   data_subset_circ <- na.omit(data_subset_circ)
#
#   #------------------------------------------------------------------------
#   # circular plot
#   #------------------------------------------------------------------------
#   p_circ <- ggplot(data = data_subset_circ,
#                    aes(
#                      x = week,
#                      xend = week + 1,
#                      y = pos,
#                      yend = pos,
#                      colour = percent_value
#                    )) +
#     scale_colour_distiller(palette = "YlOrBr",
#                            direction = 1,
#                            name = "mean annual \n% individuals\nwith events",
#                            values = leg_gradient) +
#     annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 2.4, alpha = .2) + #jan - feb
#     annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 2.4, alpha = .2) + # jun-jul
#     annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 2.4, alpha = .2) + # dec
#     annotate("text", x = 1, y = 0.8, label = "LD", col = "grey50", size = 3.5) +
#     annotate("text", x = 25, y = 0.8, label = "SD", col = "grey50", size = 3.5) +
#     annotate("text", x = 37, y = 0.8, label = "LW", col = "grey50", size = 3.5) +
#     annotate("text", x = 13, y = 0.8, label = "SW", col = "grey50", size = 3.5) +
#
#     geom_segment(size = 4) +
#     scale_x_continuous(limits = c(1,49),
#                        breaks = seq(1,48,4),
#                        labels = month.abb) +
#     scale_y_continuous(limits = c(0,2.4)) +
#     coord_polar() +
#     labs(x="",
#          y="",
#          title = title_name) +
#     theme(panel.grid.major.x = element_line(colour = "grey89",
#                                             size = 0.3),
#           panel.grid.minor.x = element_blank(),
#           panel.grid.major.y = element_blank(),
#           panel.background = element_blank(),
#
#           plot.background = element_rect(fill = 'white', colour = 'white'),
#           plot.title = element_text(hjust = 0),
#
#           axis.ticks.y = element_blank(),
#           axis.text.y = element_blank(),
#           axis.text.x = element_text(size = 10),
#           strip.text = element_text(hjust = 0, size = 13, face = "italic"),
#           strip.background = element_rect(fill="white"),
#           legend.position = leg_pos,
#           legend.box.margin=margin(c(50,50,50,150)),
#           plot.margin = unit(c(0.5,3,0,0.5),"cm")
#     ) +
#     facet_wrap(~ dummy, ncol=2, dir="v")
#
#
#   plot(p_circ)
# }

