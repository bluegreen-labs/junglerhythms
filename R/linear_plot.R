#' Create linear plot of phenology data
#'
#' @param data junglerhythms data file
#' @param species_name end coordinate fo a line
#' @export
#' @return ggplot object

linear_plot <- function(
  data,
  species_name = "Millettia laurentii"
  ){

  data_subset = data %>% filter(grepl(tolower(species_name),
                                      tolower(paste(genus, species)))) %>%
    mutate(value = ifelse(value == 0, 0, 1))

  data_subset$date <- as.Date(paste(data_subset$year,
                                    round(data_subset$week*7.6),sep="-"), "%Y-%j")

  # number of individuals
  l <- length(unique(data_subset$id))

  # reorder stuff
  loc <- which(data_subset$phenophase == "flowers")
  data_subset$value[loc] <- data_subset$value[loc] * 2.5

  loc <- which(data_subset$phenophase == "fruit")
  data_subset$value[loc] <- data_subset$value[loc] * 2

  loc <- which(data_subset$phenophase == "fruit_drop")
  data_subset$value[loc] <- data_subset$value[loc] * 1.5

  loc <- which(data_subset$phenophase == "senescence")
  data_subset$value[loc] <- data_subset$value[loc] * 1

  p = ggplot(data_subset, aes(x = date,
                              y = value,
                              colour = phenophase,
                              group = phenophase)) +
    geom_point(stat = "identity") +
    theme_minimal() +
    labs(y = "Individual",
         x = "Year",
         title = species_name) +
    theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.background = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.x = element_line(colour = "grey89", size = 0.3),
          strip.text.y = element_text(angle = 180),
          legend.position="none",
          plot.title = element_text(size=22),
          text = element_text(size=20)) +
    scale_y_continuous(breaks = NULL,
                       limits = c(1,2.5)) +
    scale_colour_ptol(labels = c("flowers",
                                 "fruit",
                                 "fruit drop",
                                 "senescence"),
                      name = "") +
    facet_wrap( ~ id, nrow = l, strip.position = "left")
  return(p)
}
