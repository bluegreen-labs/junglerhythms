#' Create circular plot of phenology data
#'
#' @param data junglerhythms data file
#' @param species_name end coordinate fo a line
#' @param threshold threshold to use when summarizing weekly data
#' @export
#' @return ggplot object

circle_plot_senescence <- function(data,
                        species_name = "Millettia laurentii",
                        weeks = 48){

  rescale <- function(x){
    if(all(x == 0 | all(is.na(x)))){
      return(x)
    } else {
      scales::rescale(x)
    }
  }

  if (weeks == 48){
    multiplier <- 7.625
  } else {
    multiplier <- 10.13889
  }

  # juggling with data
  data$species_full <- paste(data$genus, data$species)

  # split out a particular species
  data_subset <- data %>% filter(
    grepl(tolower(species_name),
          tolower(species_full)),
    phenophase == "leaf_dormancy" | phenophase == "leaf_turnover")

  # group by week and take the mean value
  # data_subset <- data_subset %>%
  #   group_by(species_full, week, phenophase) %>%
  #   summarise(mean_value = mean(value, na.rm = TRUE)) %>%
  #   ungroup()

  data_subset <- data_subset %>%
    group_by(species_full, week, phenophase) %>%
    summarise(mean_value = mean(value, na.rm = TRUE)) %>%
    mutate(mean_value = rescale(mean_value),
         pos = ifelse(phenophase == "leaf_turnover", 3, 2),
         month = months(
           as.Date(paste("1952",
                         round(week * multiplier),sep="-"),
                   "%Y-%j"))) %>%
    ungroup()

  leaf_dormancy <- data_subset %>%
    filter(phenophase == "leaf_dormancy") %>%
    mutate(mean_value = rescale(mean_value),
           pos = 3,
           month = months(
             as.Date(paste("1952",
                           round(week * multiplier),sep="-"),
                     "%Y-%j"))) %>%
    ungroup()

  leaf_turnover <- data_subset %>%
    filter(phenophase == "leaf_turnover") %>%
    mutate(mean_value = rescale(mean_value),
           pos = 2,
           month = months(
             as.Date(paste("1952",
                           round(week * multiplier),sep="-"),
                     "%Y-%j"))) %>%
    ungroup()

  # data season
  if (weeks == 48){
    x_vals <- c(33:36,1:6,NA,21:24,NA)
  }else{
    x_vals <- c(33:36,1:6,NA,21:24,NA)
  }

  # create plot and return ggplot plot object
  p <- ggplot() +
    geom_line(aes(x = x_vals,
                   #xend = week + 1,
                   y = rep(2.1,length(x_vals))
                   #yend = pos
                 ), colour = "grey",
                  na.rm = FALSE) +
    geom_segment(data = leaf_turnover,
                 aes(
                   x = week,
                   xend = week,
                   y = 0, # + 1
                   yend = mean_value
                   #fill = mean_value
                 )) +
    #scale_fill_gradient("% leaf turnover", low = "white", high = "green") +
    geom_line(data = leaf_dormancy,
                 aes(
                   x = week,
                   #xmax = week + 1,
                   y = mean_value # + 1
                   #ymax = pos + 1,
                   #colour = mean_value
                 ), linetype = "dashed") +
    #scale_colour_gradient("% leaf dormancy", low = "white", high = "red") +
    scale_x_continuous(limits = c(1,49),
                       breaks = seq(1,48,4),
                       labels = month.abb) +
    scale_y_continuous(limits = c(0,2.4)) +
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
          legend.position="right",
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)
    )

  #+
  #  facet_wrap(~ species_full)
  return(p)
}

data <- readRDS("data/jungle_rhythms_weekly_annotations.rds")

# create plot
p <- circle_plot_senescence(data = data,
                            species_name = "Carapa procera",
                            weeks = 36)

#pdf("~/Desktop/test.pdf",25,25)
plot(p)
#dev.off()
