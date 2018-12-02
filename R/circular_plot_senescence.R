#' Create circular plot of phenology data
#'
#' @param data junglerhythms data file
#' @param species_name end coordinate fo a line
#' @param threshold threshold to use when summarizing weekly data
#' @export
#' @return ggplot object

circle_plot_senescence <- function(data,
                        species_name = "Millettia laurentii"){

  rescale <- function(x){
    if(all(x == 0 | all(is.na(x)))){
      return(x)
    } else {
      scales::rescale(x)
    }
  }

  # juggling with data
  data$species_full <- paste(data$genus, data$species)

  # split out a particular species
  data_subset <- data %>% filter(
    grepl(tolower(species_name),
          tolower(species_full)),
    phenophase == "senescence")

  # group by week and take the mean value
  data_subset <- data_subset %>%
    group_by(species_full, week) %>%
    summarise(mean_value = mean(value, na.rm = TRUE))

  #data_subset <- data_subset %>%
  #  mutate(mean_value = rescale(mean_value)) %>%
  #  ungroup()

  # locate positions on a y-axis
  data_subset$pos <- 3
  data_subset$variable <- "senescence"

  # add a month column
  data_subset$month <- months(
    as.Date(paste("1952",
                  round(data_subset$week*7.625),sep="-"), "%Y-%j"))

  # data season
  data_season <- data_subset
  data_season <- data_season %>%
    mutate(sel = ifelse(week > 44 | week <= 8, 1, NA)) %>%
    mutate(week = ifelse((week > 24 & week <= 28) | sel == 1, week, NA))

  data_season$variable <- "dry season"
  data_season$pos <- 2
  data_season$mean_value <- 1
  data_season$col <- "green"

  # create plot and return ggplot plot object
  p <- ggplot() +
    geom_segment(data = data_season,
                 aes(
                   x = week,
                   xend = week + 1,
                   y = pos,
                   yend = pos
                 ), size = 10, colour = "grey") +
    geom_segment(data = data_subset,
                 aes(
                   x = week,
                   xend = week + 1,
                   y = pos,
                   yend = pos,
                   colour = mean_value
                 ),
                 size = 10) +
    scale_colour_gradient("% Senescence", low = "white", high = "red") +
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
          #axis.text.x = element_text(size = 20),
          legend.position="right"
          #legend.text = element_text(size = 15),
          #legend.title = element_text(size = 15)
    ) +
    facet_wrap(~species_full)
  return(p)
}

# # read data
# data <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
# data$species_full <- paste(data$genus, data$species)
#
# # test
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
sp <- na.omit(read.table("~/Downloads/Specieslist.csv",header = TRUE, sep = ",",
                 stringsAsFactors = FALSE))
query <- paste(sp$Species,collapse = "|")

# create plot
p <- circle_plot_senescence(data = data,
                            species_name = query)

pdf("~/Desktop/test.pdf",25,25)
  plot(p)
dev.off()
