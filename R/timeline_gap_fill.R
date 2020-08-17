#' First function:
#' Growing timelines of individuals during observation years with no single phenophase event.
#' We keep the startyear as the first year during which any phenophase event is recorded
#' and the endyear as the last year during which any phenophase event is recorded.
#' Years within this timeline without phenophase events observed have been removed as missing.
#' However, it is not unusual for these years without phenophase events to exist.
#' Therefore we can still assume these years as observed but without events.
#' To still be cautious about filling in these years, we only fill years part of 2 years of consecutive years without phenophase events.
#' If this is longer, we keep NA
#'
#' @param data junglerhythms data file
#' @param species_name
#' @param pheno only one phenophase
#' @export
#' @return datasets with filled gapyears


missing_year_gaps <- function(
  data = data,
  species_name = "Afzelia bipindensis",
  pheno = "leaf_turnover",
  gapfill_missingyears = 2){

  timelines_output <- data.frame()

  for (j in 1:length(species_name)){

    data_subset <- data %>%
      filter(species_full %in% species_name[j]) %>%
      filter(phenophase %in% pheno)

# data_subset <- data2 %>%
#   filter(species_full %in% "Afzelia bipindensis") %>%
#   filter(phenophase == "leaf_dormancy") %>%
#   filter(id %in% "2702")
    # convert date
    data_subset$date <- as.Date(paste(data_subset$year,
                                      round((data_subset$week*7.6)-7),sep="-"), "%Y-%j")
    # sort dataframe according to date
    data_subset <- data_subset %>%
      dplyr::arrange(date)

    #---
    # grow dataset to full range:
    # if consecutive years of missing data is 2 years -> fill with zero -> this needs to be done at ID level first
    #---
    individuals <- unique(data_subset$id)
    data_grow <- data.frame()
    for (i in 1:length(individuals)){
      data_ind_grow <- data_subset %>%
        filter(id %in% individuals[i])
      # get list of initial years in reduced dataset
      initial_years <- unique(data_ind_grow$year)
      # grow reduced dataset to full range based on min and max years
      years <- sort(rep(min(data_ind_grow$year):max(data_ind_grow$year), 48))
      days <- rep(round((1:48 * 7.6) - 7),length(unique(years)))
      dates_full <- data.frame(date = as.Date(paste(years, days, sep = "-"), "%Y-%j"))
      data_ind_grow <- merge(data_ind_grow, dates_full, by = "date", all.y = TRUE)
      # fill 'grown datasets' with species name, phenophase, id and year, in empty years
      data_ind_grow$species_full <- unique(na.omit(data_ind_grow$species_full))
      data_ind_grow$phenophase <- unique(na.omit(data_ind_grow$phenophase))
      data_ind_grow$id <- unique(na.omit(data_ind_grow$id))
      data_ind_grow$year <- format.Date(data_ind_grow$date, "%Y")
      data_ind_grow$week <-  ifelse(is.na(data_ind_grow$week), c(1:48), data_ind_grow$week)

      # get list of years in grown dataset
      full_range_years <- as.numeric(unique(data_ind_grow$year))

      # get list of years that were missing in the initial dataset and group those that are consecutive
      missing_years <- setdiff(full_range_years, initial_years)
      # if missing years found, run statement to fill in with zero if consecutive period is limited to 2 years
      if(length(missing_years) > 0){
        missing_years_consec <- cumsum(c(1, abs(missing_years[-length(missing_years)] - missing_years[-1]) > 1))
        missing_years <- as.data.frame(cbind(missing_years, missing_years_consec))
        colnames(missing_years)[1] <- 'year'
        missing_years_length <- missing_years %>%
          group_by(missing_years_consec) %>%
          dplyr::summarise(lgh_consec_years = length(missing_years_consec))
        missing_years <- merge(missing_years, missing_years_length, by = "missing_years_consec", all.x = TRUE)
        # merge with full range dataset
        data_ind_grow <- merge(data_ind_grow, missing_years, by = "year", all.x = TRUE)

        # if only 2 consecutive years of missing data -> fill with zeros. If longer, keep NA
        data_ind_grow$value <- ifelse(is.na(data_ind_grow$value) & data_ind_grow$lgh_consec_years <= gapfill_missingyears , 0, data_ind_grow$value) #< 3
      } else {
        data_ind_grow$missing_years_consec <- "NA"
        data_ind_grow$lgh_consec_years <- "NA"
      }
      data_grow <- rbind(data_grow, data_ind_grow)
    }

    # #-----------------------------
    # # overview plots (from Bush)
    # #-----------------------------
    # # timeseries plot of raw data for each individual
    # grown_plots <- ggplot(data_grow,
    #                       aes(x = date,
    #                           y = value))+
    #   geom_line()+
    #   scale_y_continuous(labels=NULL)+
    #   ylab("Observations") +
    #   xlab("") +
    #   # ggtitle(paste(species_name[j], "-", pheno)) + #species_name[j]
    #   theme_light(base_size = 14) +
    #   theme(legend.position = "none") +
    #   facet_wrap(~id,
    #              ncol = 1,
    #              strip.position = "left")
    #
    # plot(grown_plots)

    timelines_output <- rbind(timelines_output, data_grow)


  }


    return(timelines_output)

}


# firsttest <- two_year_gaps(data,
#                             species_name = c("Albizia adianthifolia","Allanblackia floribunda","Afrostyrax lepidophyllus"), #"Erythrophleum suaveolens",Irvingia grandifolia, Pericopsis elata
#                             pheno = "leaf_turnover") #flowers, leaf_turnover, leaf_dormancy, fruit, fruit_drop
#
# data$date <- as.Date(paste(data$year,
#                                   round((data$week*7.6)-7),sep="-"), "%Y-%j")
# alb1 <- data %>%
#   filter(species_full %in% "Albizia adianthifolia") %>%
#   filter(phenophase %in% "leaf_turnover") %>%
#   group_by(date) %>%
#   dplyr::summarise(mean_value = mean(value))
# ggplot(alb1,
#        aes(x = date,
#            y = mean_value))+
#   geom_line()+
#   scale_y_continuous(labels=NULL)+
#   ylab("Observations") +
#   xlab("") +
#   theme_light(base_size = 14) +
#   theme(legend.position = "none")
#
# alb2 <- firsttest %>%
#   filter(species_full %in% "Albizia adianthifolia") %>%
#   filter(phenophase %in% "leaf_turnover") %>%
#   group_by(date) %>%
#   dplyr::summarise(mean_value = mean(value))
# ggplot(alb2,
#        aes(x = date,
#            y = mean_value))+
#   geom_line()+
#   scale_y_continuous(labels=NULL)+
#   ylab("Observations") +
#   xlab("") +
#   theme_light(base_size = 14) +
#   theme(legend.position = "none")
