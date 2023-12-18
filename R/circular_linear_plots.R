#' Create species-specific linear and/or circular plot of phenology data
#' Linear plots show full time-series
#' Circular plots show annual means
#'
#' @param data junglerhythms data file
#' @param species_name list of species
#' @param leg_position legend position (might need to change depending on how many species called)
#' @param title_name title of the plot
#' @export
#' @return ggplot object

library(ggnewscale)

circular_linear_plot <- function(
  data = data,
  species_name = "Millettia laurentii",
  leg_pos = c(1,0.2),
  title_name = "(a) evergreen"
){

  #------------------------------------------------------------------------
  # data - circular plots
  #------------------------------------------------------------------------

  # split out a particular species
  data_subset_circ <- data %>%
    filter(grepl(tolower(species_name),
                 tolower(species_full)))

  # group by week and take the mean value
  data_subset_circ <- data_subset_circ %>%
    group_by(species_full, week, phenophase) %>%
    dplyr::summarise(percent_value = mean(value, na.rm=TRUE) *100)

  # add a month column
  data_subset_circ$month <- months(as.Date(paste("1952",
                                            round((data_subset_circ$week*7.625)-7),sep="-"), "%Y-%j"))

  # locate positions on the y-axis for the circular plots
  data_subset_circ$pos <- NA
  data_subset_circ$pos[data_subset_circ$phenophase == "leaf_dormancy"] <- 2
  data_subset_circ$pos[data_subset_circ$phenophase == "leaf_turnover"] <- 1.5
  data_subset_circ <- na.omit(data_subset_circ)

  # separate dormancy and turnover, to get different scale colors in the fig
  data_subset_circ_turn <- data_subset_circ %>%
    filter(phenophase %in% 'leaf_turnover')
  data_subset_circ_dorm <- data_subset_circ %>%
    filter(phenophase %in% 'leaf_dormancy')

  #------------------------------------------------------------------------
  # data leaf dormancy - linear plots
  #------------------------------------------------------------------------
  # subset on species
  data_subset_lin_LD <- data %>%
    filter(grepl(tolower(species_name), tolower(species_full)),
           phenophase == "leaf_dormancy")

  ids <- data_subset_lin_LD %>% # to add in facet_wrap titles (n = xx)
    group_by(species_full) %>%
    dplyr::summarise(ids = length(unique(id)))

  # average by date
  data_subset_lin_LD <- data_subset_lin_LD %>%
    group_by(species_full, date) %>%
    dplyr::summarise(percent_value = mean(value, na.rm=TRUE) *100,
              scaled_value = ifelse(any(value > 0), 1, 0))

  data_subset_lin_LD <- merge(data_subset_lin_LD, ids, by = c("species_full"), all.x = TRUE)
  data_subset_lin_LD$sp_title <- paste(data_subset_lin_LD$species_full, " (n = ", data_subset_lin_LD$ids, ")", sep = "")

  #------------------------------------------------------------------------
  # data leaf turnover - linear plots
  #------------------------------------------------------------------------
  # subset on species
  data_subset_lin_LT <- data %>%
    filter(grepl(tolower(species_name), tolower(species_full)),
           phenophase == "leaf_turnover")

  # average by date
  data_subset_lin_LT <- data_subset_lin_LT %>%
    group_by(species_full, date) %>%
    dplyr::summarise(percent_value = mean(value, na.rm=TRUE) *100,
              scaled_value = ifelse(any(value > 0), 1, 0))

  data_subset_lin_LT <- merge(data_subset_lin_LT, ids, by = c("species_full"), all.x = TRUE)
  data_subset_lin_LT$sp_title <- paste(data_subset_lin_LT$species_full, " (n = ", data_subset_lin_LT$ids, ")", sep = "")


  #------------------------------------------------------------------------
  # circular plot
  #------------------------------------------------------------------------
  p_circ <- ggplot() +
    annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 2.4, alpha = .2) + #jan - feb
    annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 2.4, alpha = .2) + # jun-jul
    annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 2.4, alpha = .2) + # dec
    annotate("text", x = 1, y = 0.8, label = "LD", col = "grey50", size = 3.5) +
    annotate("text", x = 25, y = 0.8, label = "SD", col = "grey50", size = 3.5) +
    annotate("text", x = 37, y = 0.8, label = "LW", col = "grey50", size = 3.5) +
    annotate("text", x = 13, y = 0.8, label = "SW", col = "grey50", size = 3.5) +
    geom_segment(data = data_subset_circ_dorm,
                 aes(
                   x = week,
                   xend = week + 1,
                   y = pos,
                   yend = pos,
                   colour = percent_value
                 ),
                 size = 4) +
    scale_colour_gradientn(colours = c("#ffffcc", "#DFC27D" ,"#BF812D", "#8C510A", "#543005",'#662506'),
                           name = "Senescence:\nmean annual \n% individuals\nwith events") +
    new_scale_color() +
    geom_segment(data = data_subset_circ_turn,
                 aes(
                   x = week,
                   xend = week + 1,
                   y = pos,
                   yend = pos,
                   colour = percent_value
                 ),
                 size = 4) +
    scale_colour_gradientn(colours = c("#ffffcc", "#66c2a4", "#41ae76", "#238b45", "#006d2c", "#00441b"),
                           name = "Turnover:\nmean annual \n% individuals\nwith events") +
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
                 limits = as.Date(c('1937-01-01','1956-12-31')),
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



circular_plot <- function(
  data = data,
  species_name = "Millettia laurentii",
  leg_pos = c(1,0.2),
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
    group_by(species_full, week, phenophase, dummy) %>%
    dplyr::summarise(percent_value = mean(value, na.rm=TRUE) *100)

  # add a month column
  data_subset_circ$month <- months(as.Date(paste("1952",
                                                 round((data_subset_circ$week*7.625)-7),sep="-"), "%Y-%j"))
  # locate positions on the y-axis
  data_subset_circ$pos <- NA
  data_subset_circ$pos[data_subset_circ$phenophase == "leaf_dormancy"] <- 2
  data_subset_circ$pos[data_subset_circ$phenophase == "leaf_turnover"] <- 1.5
  data_subset_circ <- na.omit(data_subset_circ)

  # separate dormancy and turnover, to get different scale colors in the fig
  data_subset_circ_turn <- data_subset_circ %>%
    filter(phenophase %in% 'leaf_turnover')
  data_subset_circ_dorm <- data_subset_circ %>%
    filter(phenophase %in% 'leaf_dormancy')

  #------------------------------------------------------------------------
  # circular plot
  #------------------------------------------------------------------------
  p_circ <- ggplot() +
    annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 2.4, alpha = .2) + #jan - feb
    annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 2.4, alpha = .2) + # jun-jul
    annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 2.4, alpha = .2) + # dec
    annotate("text", x = 1, y = 0.8, label = "LD", col = "grey50", size = 3.5) +
    annotate("text", x = 25, y = 0.8, label = "SD", col = "grey50", size = 3.5) +
    annotate("text", x = 37, y = 0.8, label = "LW", col = "grey50", size = 3.5) +
    annotate("text", x = 13, y = 0.8, label = "SW", col = "grey50", size = 3.5) +
    geom_segment(data = data_subset_circ_dorm,
                 aes(
                   x = week,
                   xend = week + 1,
                   y = pos,
                   yend = pos,
                   colour = percent_value
                 ),
                 size = 4) +
    scale_colour_gradientn(colours = c("#ffffcc", "#DFC27D" ,"#BF812D", "#8C510A", "#543005",'#662506'),
                           name = "Senescence:\nmean annual \n% individuals\nwith events") +
    new_scale_color() +
    geom_segment(data = data_subset_circ_turn,
                 aes(
                   x = week,
                   xend = week + 1,
                   y = pos,
                   yend = pos,
                   colour = percent_value
                 ),
                 size = 4) +
    scale_colour_gradientn(colours = c("#ffffcc", "#66c2a4", "#41ae76", "#238b45", "#006d2c", "#00441b"),
                           name = "Turnover:\nmean annual \n% individuals\nwith events") +
    scale_x_continuous(limits = c(1,49),
                       breaks = seq(1,48,4),
                       labels = month.abb) +
    scale_y_continuous(limits = c(0,2.4)) +
    coord_polar() +
    labs(x="",
         y="",
         title = title_name) +
    theme(panel.grid.major.x = element_line(colour = "grey89",
                                            size = 0.3),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.background = element_blank(),

          plot.background = element_rect(fill = 'white', colour = 'white'),
          plot.title = element_text(hjust = 0),

          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 10),
          strip.text = element_text(hjust = 0, size = 13, face = "italic"),
          strip.background = element_rect(fill="white"),
          legend.position = leg_pos,
          legend.box.margin=margin(c(50,50,50,150)),
          plot.margin = unit(c(0.5,3,0,0.5),"cm")
    ) +
    facet_wrap(~ dummy, ncol=2, dir="v")


  plot(p_circ)
}
