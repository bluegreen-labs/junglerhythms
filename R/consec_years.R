#' for fourier analysis, no missing years can be incorporated.
#' If missing years within a timeline, then only the longest consecutive timeline is kept.
#'
#' @param data junglerhythms data file
#' @param species_name list of species
#' @param pheno only one phenophase
#' @export
#' @return dataframe timelines at species-level, no gaps


consecutive_timeline_sp <- function(
  data = data,
  species_name = "Afzelia bipindensis",
  pheno = "leaf_turnover"){

  timelines_output <- data.frame()

  for (j in 1:length(species_name)){

    data_subset <- data %>%
      filter(species_full %in% species_name[j]) %>%
      filter(phenophase %in% pheno)

    #-------------------
    # SPECIES LEVEL
    # average by date
    #-------------------
    data_sp <- data_subset %>%
      group_by(date) %>%
      dplyr::summarise(mean_value = mean(value,na.rm=T))

    data_sp$mean_value[is.nan(data_sp$mean_value)] <- NA
    data_sp$scaled_value <- ifelse(data_sp$mean_value > 0, 1,
                                   ifelse(data_sp$mean_value == 0, 0, NA))

    # sort dataframe according to date and add a month and year column
    data_sp <- data_sp %>%
      dplyr::arrange(date)
    data_sp$month <- format.Date(data_sp$date, "%m")
    data_sp$year <- format.Date(data_sp$date, "%Y")

    #-------------------
    # for fourier at SPECIES LEVEL
    # if timeline not continuous
    # only keep longest timeline
    #-------------------
    ### if longer periods of missing years are found, the longest strech of consecutive timelines for an individual is kept
    timeline <- data_sp[!(is.na(data_sp$mean_value)),]
    years_timeline <- as.numeric(unique(timeline$year))
    consec_timelines <- cumsum(c(1, abs(years_timeline[-length(years_timeline)] - years_timeline[-1]) > 1))
    years_timeline <- as.data.frame(cbind(years_timeline,consec_timelines))
    colnames(years_timeline)[1] <- 'year'
    years_timeline_length <- years_timeline %>%
      group_by(consec_timelines) %>%
      dplyr::summarise(lgh_consec_timelines = length(consec_timelines))
    years_timeline <- merge(years_timeline, years_timeline_length, by = "consec_timelines", all.x = TRUE)
    timeline <- merge(timeline, years_timeline, by = "year", all.x = TRUE)
    # only keep the years which together form the longest consecutive timeline
    timeline$mean_value <- ifelse(timeline$lgh_consec_timelines == max(timeline$lgh_consec_timelines) , timeline$mean_value, NA)
    # if two or more timelines non-consecutive but equally long
    if(length(unique(timeline$lgh_consec_timelines)) == 1 & length(unique(timeline$consec_timelines)) > 1){
      test <- timeline %>%
        group_by(consec_timelines) %>%
        dplyr::summarise(observations = sum(mean_value))
      pos_max <- which.max(test$observations)
      timeline$mean_value <- ifelse(timeline$consec_timelines ==  pos_max, timeline$mean_value, NA)
    }
    # remove all years with NA value
    timeline <- timeline[!(is.na(timeline$mean_value)),]

    # add species name and phenophase
    timeline$species_full <- species_name[j]
    timeline$phenophase <- pheno

    timelines_output <- rbind(timelines_output, timeline)
  }

  return(timelines_output)

}

