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
  # viridis_rescaling = 0.15,
  leg_pos = c(1,0.1),
  leg_gradient = c(0,0.2,1),
  title_name = "(a) evergreen"
){

  #------------------------------------------------------------------------
  # data for circular plots
  #------------------------------------------------------------------------
  # add site years, + individuals
  # count sites / years
  # data <- data %>%
  #   group_by(species_full, id) %>%
  #   mutate(S = length(unique(year))) %>%
  #   ungroup()
  #
  # data <- data %>%
  #   group_by(species_full) %>%
  #   mutate(N = length(unique(id)),
  #          S = sum(S))
  #
  # data$species_full <- paste0(data$species_full,"\n",
  #                             "(N = ",data$N,", S = ",
  #                             data$S,")")

  # split out a particular species
  data_subset_circ <- data %>%
    filter(grepl(tolower(species_name),
                 tolower(species_full)))

  # group by week and take the mean value
  data_subset_circ <- data_subset_circ %>%
    group_by(species_full, week, phenophase) %>%
    dplyr::summarise(mean_value = mean(value, na.rm=TRUE))#,
              # N = unique(N),
              # S = unique(S))

  # # rescale
  # rescale <- function(x){
  #   if(all(x == 0)){
  #     return(x)
  #   } else {
  #     scales::rescale(x)
  #   }
  # }
  # data_subset_circ <- data_subset_circ %>%
  #   group_by(species_full, phenophase) %>%
  #   mutate(mean_value = rescale(mean_value))

  # add a month column
  data_subset_circ$month <- months(as.Date(paste("1952",
                                            round(data_subset_circ$week*7.625),sep="-"), "%Y-%j"))

  # locate positions on a y-axis
  data_subset_circ$pos <- NA
  data_subset_circ$pos[data_subset_circ$phenophase == "leaf_dormancy"] <- 2
  data_subset_circ$pos[data_subset_circ$phenophase == "leaf_turnover"] <- 1.5
  data_subset_circ <- na.omit(data_subset_circ)


  #------------------------------------------------------------------------
  # climate data
  #------------------------------------------------------------------------
  # # read climate data
  # climate <- readRDS("data/yangambi_monthly_climate.rds")
  # climate <- climate %>% filter(year < 1960,
  #                               measurement == "precip") %>%
  #   mutate(col = ifelse(value < drought_mm,"dry","wet"))

  #------------------------------------------------------------------------
  # linear leaf dormancy
  #------------------------------------------------------------------------
  # subset on species
  data_subset_lin_LD <- data %>%
    filter(grepl(tolower(species_name), tolower(species_full)),
           phenophase == "leaf_dormancy")

  species_nr <- length(unique(data_subset_lin_LD$species_full))
  if (species_nr == 1){
    species_nr = 4
  }

  # convert date
  data_subset_lin_LD$date <- as.Date(
    paste(data_subset_lin_LD$year,
          round(data_subset_lin_LD$week*7.6),sep="-"), "%Y-%j")

  # average by week
  data_subset_lin_LD <- data_subset_lin_LD %>%
    group_by(species_full, date) %>%
    dplyr::summarise(mean_value = mean(value),
              scaled_value = ifelse(any(value > 0), 1, 0))

  #------------------------------------------------------------------------
  # linear leaf turnover
  #------------------------------------------------------------------------
  # subset on species
  data_subset_lin_LT <- data %>%
    filter(grepl(tolower(species_name), tolower(species_full)),
           phenophase == "leaf_turnover")
  species_nr2 <- length(unique(data_subset_lin_LT$species_full))
  if (species_nr == 1){
    species_nr = 4
  }
  # convert date
  data_subset_lin_LT$date <- as.Date(
    paste(data_subset_lin_LT$year,
          round(data_subset_lin_LT$week*7.6),sep="-"), "%Y-%j")
  # average by week
  data_subset_lin_LT <- data_subset_lin_LT %>%
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
    # scale_color_viridis_c(option="plasma",
    #                       name = "Freq.",
    #                       direction = 1,
    #                       begin = 0,
    #                       end = 1) +
    # scale_color_viridis_c(option = "plasma", #,"viridis"
    #                       name = "Freq.",
    #                       direction = 1,
    #                       rescaler = function(x, to = c(0, 1), from = NULL) {
    #   ifelse(x < viridis_rescaling,
    #          scales::rescale(x,
    #                          to = to,
    #                          from = c(min(x, na.rm = TRUE), viridis_rescaling)),
    #          1)}) +
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
          strip.text = element_blank(), #element_text(face = "italic", size = 13),
          strip.background = element_rect(fill="white"),
          legend.position = leg_pos,
          legend.box.margin=margin(c(50,50,50,150)),
          plot.margin = unit(c(0,3,0.7,0.5),"cm") #unit(c(0,0.5,0.7,0.5)
    ) +
    facet_wrap(~ species_full,ncol=1)

  #------------------------------------------------------------------------
  # linear plot
  #------------------------------------------------------------------------
  p_lin <- ggplot(data_subset_lin_LD, aes(x = date,
                               y = mean_value)) +
    geom_line(stat = "identity",aes(color="Canopy dormancy")) +
    geom_line(data=data_subset_lin_LT,aes(color="Canopy turnover")) +
    scale_colour_manual(values = c("black","grey60")) +
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
          plot.margin = unit(c(0.1,0,1,1),"cm")#,
          # panel.margin.y = unit(2, "lines")
    ) +
    facet_wrap( ~ species_full, ncol = 1)

  # p_bottom <- ggplot() +
  #   geom_col(data = climate, aes(x = date,
  #                                y = value,
  #                                colour = col,
  #                                fill = col)) +
  #   scale_colour_manual(values = c("lightcoral","lightblue"),
  #                       aesthetics = c("fill","colour")) +
  #   labs(y = "precip. (mm/month)",
  #        x = "Year") +
  #   scale_x_date(date_breaks = "1 year",
  #                date_labels = "%Y",
  #                limits = as.Date(c('1935-01-01','1960-01-01'))) +
  #   theme_minimal() +
  #   theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
  #         panel.grid.minor.x =  element_blank(),
  #         panel.grid.minor.y = element_blank(),
  #         panel.background = element_blank(),
  #         plot.background = element_rect(fill = 'white', colour = 'white'),
  #         strip.text = element_text(hjust = 0),
  #         axis.line.x = element_blank(),
  #         axis.text.x = element_text(angle = 90, hjust = 1),
  #         legend.position="none"
  #   )

  return(grid.arrange(p_lin, p_circ, widths = c(2,1.5)))
}


# p_deciduous <- circular_linear_plot(data,
#                                     species_name = "Entandrophragma angolense",
#                                     viridis_rescaling = 0.15,
#                                     title_name = "(b) Deciduous")
# pdf("~/Desktop/deciduous.pdf",8.4,12)
# plot(p_deciduous)
# dev.off()

# p_evergreen <- circular_linear_plot(data,
#                                     species_name = query_evergreen,
#                                     viridis_rescaling = 0.05,
#                                     title_name = "(a) Evergreen")
# pdf("~/Desktop/evergreen_viridis.pdf",8.4,14.4)
# plot(p_evergreen)
# dev.off()
