#' Plot weekly annotations
#'
#' Visualizes the weekly annotations
#'
#' @param df annotations for a particular yearly section
#' @export
#'

plot_weekly_annotations <- function(data){

  # reorder stuff
  loc <- which(data$phenophase == "flowers")
  data$value[loc] <- data$value[loc] * 2.5

  loc <- which(data$phenophase == "fruit")
  data$value[loc] <- data$value[loc] * 2

  loc <- which(data$phenophase == "fruit_drop")
  data$value[loc] <- data$value[loc] * 1.5

  loc <- which(data$phenophase == "senescence")
  data$value[loc] <- data$value[loc] * 1

  # plot things
  p <- ggplot(data, aes(x = week,
                        y = value,
                        colour = phenophase,
                        group = as.factor(phenophase))) +
    geom_point(stat = "identity") +
    theme_minimal() +
    labs(x = "week", y = "") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.background = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.x = element_line(colour = "grey89", size = 0.3),
          strip.text.y = element_text(angle = 180),
          legend.position="right",
          plot.title = element_text(size=22),
          text = element_text(size=20)) +
    scale_x_continuous(breaks = seq(0,48,4),
                       limits = c(0,49),
                       labels = ) +
    geom_vline(xintercept = seq(0,47,4) + 0.5, colour = "grey") +
    geom_vline(xintercept = c(0.5,24.5,48.5), colour = "black") +
    scale_y_continuous(breaks = NULL,
                       limits = c(1,2.5)) +
    scale_colour_ptol(labels = c("flowers",
                                 "fruit",
                                 "fruit drop",
                                 "leaf dormancy",
                                 "leaf turnover"),
                      name = "")
  return(p)
}

