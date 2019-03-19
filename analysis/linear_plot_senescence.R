#' Create linear plot of phenology data
#'
#' @param data junglerhythms data file
#' @param species_name end coordinate fo a line
#' @export
#' @return ggplot object

library(ggplot2)
library(ggthemes)
library(gridExtra)
library(tidyverse)

linear_plot_senescence <- function(
  data,
  species_name = "Millettia laurentii",
  drought_mm = 140
  ){

  # read climate data
  climate <- readRDS("data/yangambi_monthly_climate.rds")
  climate <- climate %>% filter(year < 1960,
                                measurement == "precip") %>%
    mutate(col = ifelse(value < drought_mm,"dry","wet"))

  # subset on species
  data_subset <- data %>%
    filter(grepl(tolower(species_name), tolower(species_full)),
    phenophase == "leaf_dormancy")

  species_nr <- length(unique(data_subset$species_full))
  if (species_nr == 1){
    species_nr = 4
  }

  # convert date
  data_subset$date <- as.Date(
    paste(data_subset$year,
    round(data_subset$week*7.6),sep="-"), "%Y-%j")

  # average by week
  data_subset <- data_subset %>%
    group_by(species_full, date) %>%
    summarize(mean_value = mean(value),
              scaled_value = ifelse(any(value > 0), 1, 0))

  p <- ggplot(data_subset, aes(x = date,
                              y = mean_value)) +
    geom_line(stat = "identity") +
    theme_minimal() +
    labs(y = "% senescence (across individuals)",
         x = "Year") +
    scale_x_date(date_breaks = "1 years",
                 date_labels = "%Y",
                 limits = as.Date(c('1935-01-01','1960-01-01'))) +
    scale_y_continuous(breaks = c(0,0.5,1)) +
    theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = 'white', colour = 'white'),
          strip.text = element_text(hjust = 0),
          axis.line.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
          ) +
  facet_wrap( ~ species_full, ncol = 1)

  p_bottom <- ggplot() +
    geom_col(data = climate, aes(x = date,
                                 y = value,
                                 colour = col,
                                 fill = col)) +
    scale_colour_manual(values = c("lightcoral","lightblue"),
                        aesthetics = c("fill","colour")) +
    labs(y = "precip. (mm/month)",
         x = "Year") +
    scale_x_date(date_breaks = "1 year",
                 date_labels = "%Y",
                 limits = as.Date(c('1935-01-01','1960-01-01'))) +
    theme_minimal() +
    theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
          panel.grid.minor.x =  element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = 'white', colour = 'white'),
          strip.text = element_text(hjust = 0),
          axis.line.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position="none"
    )

  return(grid.arrange(p,p_bottom, heights = c((species_nr-2)/species_nr,
                                              2/species_nr)))
}

# read data
data <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
data$species_full <- paste(data$genus, data$species)
#
# # # test
# data <- data %>%
#   group_by(species_full, year) %>%
#   mutate(selection = if(all(value == 0)){
#    NA
#   }else{
#     value
#   }) %>%
#   na.omit() %>%
#   ungroup()

# query
sp <- na.omit(read.csv2("data/Specieslist.csv",
                         header = TRUE, sep = ",",
                         stringsAsFactors = FALSE))#[1:10,]
query <- paste(sp$Species,collapse = "|")


# create plot
p <- linear_plot_senescence(data = data, species_name = "Croton macrostachyus")

#pdf("~/Desktop/test.pdf",25,25)
plot(p)
#dev.off()


