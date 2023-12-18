#' plot-level annual phenological signal
#' species-specific annual signals are weighted by their basal area at the plot level
#'
#' @param data junglerhythms data file
#' @param census_plot yangambi census data at plot level
#' @param total_basal_area_plot total basal area of a plot
#' @param species_name list of species
#' @param pheno only one phenophase
#' @param minimum_siteyears species not included if fewer observation-years overall (across all individuals)
#' @export
#' @return dataframe


standlevel_phen_plotlevel <- function(
    data = data,
    census_plot = census_plot,
    total_basal_area_plot = total_basal_area_plot,
    species_list_dorm = dorm_sp1,
    species_list_turn = turn_sp1,
    minimum_siteyears = minimum_siteyears
){


  #-----------------------------------------------------------------------
  #-- each phenophase is done seperately
  #-- because different species can be sourced per phenophase
  #-- based on group classifications
  #-----------------------------------------------------------------------

  #-----------------------------------------------------------------------
  #------------ Leaf turnover  -------------------------------------------
  #-----------------------------------------------------------------------
  if(!is_empty(species_list_turn)){
    timelines_turn <- data %>%
      filter(phenophase == "leaf_turnover",
             species_full %in% species_list_turn,
             !is.na(value))
    # for each species, get summaries of each week (makes full year per species)
    data_LT <- timelines_turn %>%
      group_by(species_full, week) %>%
      dplyr::summarise(mean_week = mean(value),
                       total_week = length(value))

    # filter out species with too few total observation years across individuals
    # to have a meaningfull average year (set by minimum_siteyears)
    data_LT <- data_LT %>%
      filter(total_week >= minimum_siteyears)

    # add species-level basal area (across the 5-ha plots)
    data_LT_plot <- inner_join(data_LT, census_plot, by = c("species_full"))


    final_LT_plot <- data_LT_plot %>%
      group_by(Plot, week) %>%
      dplyr::summarise(ss = sum(mean_week*basal_area_plot, na.rm = TRUE))
    final_LT_plot <- inner_join(final_LT_plot, total_basal_area_plot, by = c("Plot"))
    final_LT_plot <- final_LT_plot %>%
      mutate(ss_turn = ss/total_basal_area_plot*100)
    final_LT_plot <- final_LT_plot %>%
      dplyr::select(Plot, week, ss_turn)
  } else {
    final_LT_plot <- data.frame(Plot = rep(c("Inventory 20","Inventory 21","Inventory 23"),48),
                                week = rep(c(1:48),2),
                                ss_turn = NA)
  }
  #-----------------------------------------------------------------------

  #-----------------------------------------------------------------------
  #------------ Leaf dormancy  -------------------------------------------
  #-----------------------------------------------------------------------
  if(!is_empty(species_list_dorm)){
    timelines_dorm <- data %>%
      filter(phenophase == "leaf_dormancy",
             species_full %in% species_list_dorm,
             !is.na(value))

    data_LD <- timelines_dorm %>%
      group_by(species_full, week) %>%
      dplyr::summarise(mean_week = mean(value),
                       total_week = length(value))
    data_LD <- data_LD %>%
      filter(total_week >= minimum_siteyears)

    data_LD_plot <- inner_join(data_LD, census_plot, by = c("species_full"))

    final_LD_plot <- data_LD_plot %>%
      group_by(Plot, week) %>%
      dplyr::summarise(ss = sum(mean_week*basal_area_plot, na.rm = TRUE))
    final_LD_plot <- inner_join(final_LD_plot, total_basal_area_plot, by = c("Plot"))
    final_LD_plot <- final_LD_plot %>%
      mutate(ss_dorm = ss/total_basal_area_plot*100)
    final_LD_plot <- final_LD_plot %>%
      dplyr::select(Plot, week, ss_dorm)
  } else {
    final_LD_plot <- data.frame(Plot = rep(c("Inventory 20","Inventory 21","Inventory 23"),48),
                                week = rep(c(1:48),2),
                                ss_dorm = NA)
  }
  #-----------------------------------------------------------------------


  # #-----------------------------------------------------------------------
  # #------------ Leaf flushing  -------------------------------------------
  # #-----------------------------------------------------------------------
  # if(!is_empty(species_list_dorm)) {
  #   #------
  #   # onset of flushing is first calculated as first week end of dormancy events
  #   # full timeline for flushing is first created for each individual
  #   #------
  #   # just to avoid confusion, work with new name
  #   timelines_flush <- timelines_dorm
  #   # new column 'flushing_date' to merge with later on
  #   timelines_flush$flushing_date <- paste(timelines_flush$year, timelines_flush$week, sep = "-")
  #
  #   # function event_length gets
  #   # for each individual of the requested species
  #   # start/end year&week of each event
  #   flushing_timing <- event_length(data = timelines_flush,
  #                                   species_name = species_list_dorm,
  #                                   pheno = "leaf_dormancy")
  #   # individuals without events give NA in event_length, so remove
  #   flushing_timing <- flushing_timing %>%
  #     filter(!is.na(year_start))
  #
  #   # date of onset flushing = week after end dormancy event
  #   flushing_timing$flushing_date <- paste(flushing_timing$year_end, flushing_timing$week_end +1, sep = "-")
  #
  #   flushing_timing <- flushing_timing %>%
  #     select(species_full,
  #            id,
  #            flushing_date)
  #   flushing_timing$flushing_value <- 1
  #
  #   # merge using flushing_date for each id/species
  #   # flushing_value is
  #   data_flush <- merge(timelines_flush, flushing_timing, by = c("species_full","id","flushing_date"), all.x = TRUE)
  #   data_flush$flushing_value <- ifelse(is.na(data_flush$flushing_value),"0", data_flush$flushing_value)
  #   # if value for dormancy was NA, year was not observed (full timelines were needed for using 'event_length')
  #   # so remove year for flushing
  #   data_flush$flushing_value <- ifelse(is.na(data_flush$value),NA, data_flush$flushing_value)
  #   data_flush$flushing_value <- as.numeric(data_flush$flushing_value)
  #
  #   #------
  #   ## now continue the same as done for turnover and dormancy
  #   #------
  #   data_flush <- data_flush %>%
  #     group_by(species_full, week) %>%
  #     dplyr::summarise(mean_week = mean(flushing_value, na.rm = TRUE),
  #                      total_week = length(flushing_value))
  #
  #   data_flush <- data_flush %>%
  #     filter(total_week >= minimum_siteyears)
  #
  #   data_flush_plot <- inner_join(data_flush, census_plot, by = c("species_full"))
  #
  #   final_flush_plot <- data_flush_plot %>%
  #     group_by(Plot, week) %>%
  #     dplyr::summarise(ss = sum(mean_week*basal_area_plot, na.rm = TRUE))
  #   final_flush_plot <- inner_join(final_flush_plot, total_basal_area_plot, by = c("Plot"))
  #   final_flush_plot <- final_flush_plot %>%
  #     mutate(ss_flush = ss/total_basal_area_plot*100)
  #   final_flush_plot <- final_flush_plot %>%
  #     dplyr::select(Plot, week, ss_flush)
  # } else {
  #   final_flush_plot <- data.frame(Plot = rep(c("Inventory 20","Inventory 23"),48),
  #                                  week = rep(c(1:48),2),
  #                                  ss_flush = NA)
  # }
  # #-----------------------------------------------------------------------

  #-----------------------------------------------------------------------
  #------  output dataframe  ---------------------------------------------
  #-----------------------------------------------------------------------
  # # add full_range to make sure all plots are included, even if a species is not present in the plot
  # combined <-  data.frame(Plot = rep(c("Inventory 20","Inventory 23"),48),
  #                         week = rep(c(1:48),2),
  #                         full_range = NA)
  # combined <- merge(combined, final_LT_plot, by = c("Plot","week"), all.x = TRUE)
  # combined <- merge(combined, final_LD_plot, by = c("Plot","week"), all.x = TRUE)
  # combined <- merge(combined, final_flush_plot, by = c("Plot","week"), all.x = TRUE)
  # # both dormancy and turnover represent main period of senescence
  # combined$ss_senescence <- combined$ss_turn + combined$ss_dorm
  # combined <- combined %>%
  #   dplyr::select(!full_range)

  combined <- left_join(final_LD_plot,final_LT_plot)

  return(combined)
}
