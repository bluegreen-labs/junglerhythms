#' Create circular plot of phenology data
#'
#' @param data junglerhythms data file
#' @param species_name end coordinate fo a line
#' @param threshold threshold to use when summarizing weekly data
#' @export
#' @return ggplot object

circle_plot <- function(data,
                        species_name = "Millettia laurentii"){

  rescale <- function(x){
    if(all(x == 0)){
      return(x)
    } else {
      scales::rescale(x)
    }
  }

  # add site years, + individuals

  # juggling with data
  data$species_full <- paste(data$genus, data$species)

  # count sites / years
  data <- data %>%
    group_by(species_full, id) %>%
    mutate(S = length(unique(year))) %>%
    ungroup()

  data <- data %>%
    group_by(species_full) %>%
    mutate(N = length(unique(id)),
           S = sum(S))

  data$species_full <- paste0(data$species_full,"\n",
                             "(N = ",data$N,", S = ",
                             data$S,")")

  # split out a particular species
  data_subset <- data %>%
    filter(grepl(tolower(species_name),
                             tolower(species_full)))

  # group by week and take the mean value
  data_subset <- data_subset %>%
    group_by(species_full, week, phenophase) %>%
    summarise(mean_value = mean(value, na.rm=TRUE),
              N = unique(N),
              S = unique(S))

  # rescale
  # data_subset <- data_subset %>%
  #   group_by(species_full, phenophase) %>%
  #   mutate(mean_value = rescale(mean_value))

  # add a month column
  data_subset$month <- months(as.Date(paste("1952",
                      round(data_subset$week*7.625),sep="-"), "%Y-%j"))

  # locate positions on a y-axis
  data_subset$pos <- NA
  data_subset$pos[data_subset$phenophase == "leaf_dormancy"] <- 2
  data_subset$pos[data_subset$phenophase == "leaf_turnover"] <- 1.5
  data_subset <- na.omit(data_subset)

  # create plot and return ggplot plot object
  p <- ggplot(data = data_subset,
              aes(
                x = week,
                xend = week + 1,
                y = pos,
                yend = pos,
                colour = mean_value
              )) +
    scale_colour_distiller(palette = "Spectral",
                           name = "fractional intensity") +
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
          strip.text = element_text(face = "italic", size = 13),
          strip.background = element_rect(fill="white"),
          legend.position = "right",
          plot.margin=unit(c(1,1,1,1),"cm")
    ) +
    facet_wrap(~ species_full)
  return(p)
}

library(tidyverse)

data <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
data$species_full <- paste(data$genus, data$species)
sp <- na.omit(read.csv2("data/Specieslist.csv",
                         header = TRUE, sep = ",",
                         stringsAsFactors = FALSE))
#query <- paste(sp$Species[1:42],collapse = "|")

p <- circle_plot(data, species_name = query)

pdf("~/Desktop/test.pdf",25,25)
plot(p)
dev.off()
