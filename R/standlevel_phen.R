#' stand-level annual phenological signal
#' species-specific annual signals are weighted by their basal area
#'
#' @param data junglerhythms data file
#' @param census_site yangambi census data at site level
#' @param total_basal_area_site
#' @param species_name list of species
#' @param pheno only one phenophase
#' @param minimum_siteyears species not included if fewer observation-years overall (across all individuals)
#' @export
#' @return dataframe


standlevel_phen <- function(
  data = data,
  census_site = census_site,
  total_basal_area_site = total_basal_area_site,
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
                       # count_week = sum(value),
                       total_week = length(value))

    # filter out species with too few total observation years across individuals
    # to have a meaningfull average year (set by minimum_siteyears)
    data_LT <- data_LT %>%
      filter(total_week >= minimum_siteyears)

    # add species-level basal area (across the 5-ha plots)
    data_LT <- inner_join(data_LT, census_site, by = c("species_full"))

    # weighted mean for standlevel phenological pattern
    final_LT <- data_LT %>%
      group_by(week) %>%
      dplyr::summarise(ss_turn = sum(mean_week*basal_area_site, na.rm = TRUE)/total_basal_area_site *100)
    # ss2 = wtd.mean(mean_week, basal_area_site, na.rm = TRUE)*100,
    # weighted_sd = sqrt(wtd.var(mean_week, basal_area_site, na.rm = TRUE))*100,
    # weighted_CI = ss2 + 1.96*weighted_sd)
  } else {
    final_LT <- data.frame(week = c(1:48),
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

    data_LD <- inner_join(data_LD, census_site, by = c("species_full"))

    final_LD <- data_LD %>%
      group_by(week) %>%
      dplyr::summarise(ss_dorm = sum(mean_week*basal_area_site, na.rm = TRUE)/total_basal_area_site *100)
  } else {
    final_LD <- data.frame(week = c(1:48),
                           ss_dorm = NA)
  }

  #-----------------------------------------------------------------------


  #-----------------------------------------------------------------------
  #------------ Leaf flushing  -------------------------------------------
  #-----------------------------------------------------------------------
  if(!is_empty(species_list_dorm)) {
    #------
    # onset of flushing is first calculated as first week end of dormancy events
    # full timeline for flushing is first created for each individual
    #------
    # just to avoid confusion, work with new name
    timelines_flush <- timelines_dorm
    # new column 'flushing_date' to merge with later on
    timelines_flush$flushing_date <- paste(timelines_flush$year, timelines_flush$week, sep = "-")

    # function event_length gets
    # for each individual of the requested species
    # start/end year&week of each event
    flushing_timing <- event_length(data = timelines_flush,
                                    species_name = species_list_dorm,
                                    pheno = "leaf_dormancy")
    # individuals without events give NA in event_length, so remove
    flushing_timing <- flushing_timing %>%
      filter(!is.na(year_start))

    # date of onset flushing = week after end dormancy event
    flushing_timing$flushing_date <- paste(flushing_timing$year_end, flushing_timing$week_end +1, sep = "-")

    flushing_timing <- flushing_timing %>%
      select(species_full,
             id,
             flushing_date)
    flushing_timing$flushing_value <- 1

    # merge using flushing_date for each id/species
    # flushing_value is
    data_flush <- merge(timelines_flush, flushing_timing, by = c("species_full","id","flushing_date"), all.x = TRUE)
    data_flush$flushing_value <- ifelse(is.na(data_flush$flushing_value),"0", data_flush$flushing_value)
    # if value for dormancy was NA, year was not observed (full timelines were needed for using 'event_length')
    # so remove year for flushing
    data_flush$flushing_value <- ifelse(is.na(data_flush$value),NA, data_flush$flushing_value)
    data_flush$flushing_value <- as.numeric(data_flush$flushing_value)

    #------
    ## now continue the same as done for turnover and dormancy
    #------
    data_flush <- data_flush %>%
      group_by(species_full, week) %>%
      dplyr::summarise(mean_week = mean(flushing_value, na.rm = TRUE),
                       total_week = length(flushing_value))
    data_flush <- data_flush %>%
      filter(total_week >= minimum_siteyears)

    data_flush <- inner_join(data_flush, census_site, by = c("species_full"))

    final_flush <- data_flush %>%
      group_by(week) %>%
      dplyr::summarise(ss_flush = sum(mean_week*basal_area_site, na.rm = TRUE)/total_basal_area_site *100)
  } else {
    final_flush <- data.frame(week = c(1:48),
                           ss_flush = NA)
  }

  #-----------------------------------------------------------------------


  #-----------------------------------------------------------------------
  #------  output dataframe  ---------------------------------------------
  #-----------------------------------------------------------------------
  combined <- merge(final_LT, final_LD, by = c("week"), all.x = TRUE)
  combined <- merge(combined, final_flush, by = c("week"), all.x = TRUE)
  # both dormancy and turnover represent main period of senescence
  combined$ss_senescence <- combined$ss_turn + combined$ss_dorm

  return(combined)
}
