standlevel_phen_onset <- function(
  data = data,
  species_list_dorm = dorm_sp1,
  species_list_turn = turn_sp1
){

  combined <- data.frame()

  timelines_dorm <- missing_year_gaps(data = data,
                                      species_name = species_list_dorm,
                                      pheno = "leaf_dormancy",
                                      gapfill_missingyears = 0)
  timelines_turn <- missing_year_gaps(data = data,
                                      species_name = species_list_turn,
                                      pheno = "leaf_turnover",
                                      gapfill_missingyears = 0)

  #-----------------------------------------------------------------------
  #------------ Leaf turnover  -------------------------------------------
  #-----------------------------------------------------------------------
  timelines_turn <- missing_year_gaps(data = data,
                                     species_name = species_list_turn,
                                     pheno = "leaf_turnover",
                                     gapfill_missingyears = 0)


  timelines_turn$turnover_date <- paste(timelines_turn$year, timelines_turn$week, sep = "-")

  turnover_timing <- event_length(data = timelines_turn,
                                    species_name = species_list_turn,
                                    pheno = "leaf_turnover")
  turnover_timing <- turnover_timing %>%
    filter(!is.na(year_start))

  turnover_timing$turnover_date <- paste(turnover_timing$year_start, turnover_timing$week_start, sep = "-")

  turnover_timing <- turnover_timing %>%
    select(species_full,
           id,
           turnover_date)
  turnover_timing$turnover_value <- 1

  data_turn <- merge(timelines_turn, turnover_timing, by = c("species_full","id","turnover_date"), all.x = TRUE)
  data_turn$turnover_value <- ifelse(is.na(data_turn$turnover_value),"0", data_turn$turnover_value)
  data_turn$turnover_value <- ifelse(is.na(data_turn$value),NA, data_turn$turnover_value)
  data_turn$turnover_value <- as.numeric(data_turn$turnover_value)

  data_turn <- data_turn %>%
    group_by(species_full, date) %>%
    dplyr::summarise(mean_week = mean(turnover_value, na.rm = TRUE),
                     count_week = sum(turnover_value),
                     total_week = length(turnover_value))

  # # filter out species with too few total siteyears to have a meaningfull average year
  # data_flush <- data_flush %>%
  #   filter(total_week >= minimum_siteyears)

  # data_flush$mean_week <- ifelse(data_flush$mean_week < minimum_event_frequency, 0, data_flush$mean_week)

  data_turn_site <- inner_join(data_turn, census_site, by = c("species_full"))


  final_turn_site <- data_turn_site %>%
    group_by(date) %>%
    dplyr::summarise(ss_turn = sum(mean_week*basal_area_site, na.rm = TRUE)/total_basal_area_site *100)
  #-----------------------------------------------------------------------


  #-----------------------------------------------------------------------
  #------------ Leaf senescence  -----------------------------------------
  #-----------------------------------------------------------------------
  timelines_sen <- missing_year_gaps(data = data,
                                     species_name = species_list_dorm,
                                     pheno = "leaf_dormancy",
                                     gapfill_missingyears = 0)


  timelines_sen$senescence_date <- paste(timelines_sen$year, timelines_sen$week, sep = "-")

  senescence_timing <- event_length(data = timelines_sen,
                                  species_name = species_list_dorm,
                                  pheno = "leaf_dormancy")
  senescence_timing <- senescence_timing %>%
    filter(!is.na(year_start))

  senescence_timing$senescence_date <- paste(senescence_timing$year_start, senescence_timing$week_start, sep = "-")

  senescence_timing <- senescence_timing %>%
    select(species_full,
           id,
           senescence_date)
  senescence_timing$senescence_value <- 1

  data_sen <- merge(timelines_sen, senescence_timing, by = c("species_full","id","senescence_date"), all.x = TRUE)
  data_sen$senescence_value <- ifelse(is.na(data_sen$senescence_value),"0", data_sen$senescence_value)
  data_sen$senescence_value <- ifelse(is.na(data_sen$value),NA, data_sen$senescence_value)
  data_sen$senescence_value <- as.numeric(data_sen$senescence_value)

  data_sen <- data_sen %>%
    group_by(species_full, date) %>%
    dplyr::summarise(mean_week = mean(senescence_value, na.rm = TRUE),
                     count_week = sum(senescence_value),
                     total_week = length(senescence_value))

  # # filter out species with too few total siteyears to have a meaningfull average year
  # data_flush <- data_flush %>%
  #   filter(total_week >= minimum_siteyears)

  # data_flush$mean_week <- ifelse(data_flush$mean_week < minimum_event_frequency, 0, data_flush$mean_week)

  data_sen_site <- inner_join(data_sen, census_site, by = c("species_full"))


  final_sen_site <- data_sen_site %>%
    group_by(date) %>%
    dplyr::summarise(ss_sen = sum(mean_week*basal_area_site, na.rm = TRUE)/total_basal_area_site *100)

  #-----------------------------------------------------------------------

  #-----------------------------------------------------------------------
  #------------ Leaf flushing  -------------------------------------------
  #-----------------------------------------------------------------------
  timelines_flush <- missing_year_gaps(data = data,
                                       species_name = species_list_dorm,
                                       pheno = "leaf_dormancy",
                                       gapfill_missingyears = 0)


  timelines_flush$flushing_date <- paste(timelines_flush$year, timelines_flush$week, sep = "-")

  flushing_timing <- event_length(data = timelines_flush,
                                  species_name = species_list_dorm,
                                  pheno = "leaf_dormancy")
  flushing_timing <- flushing_timing %>%
    filter(!is.na(year_start))

  flushing_timing$flushing_date <- paste(flushing_timing$year_end, flushing_timing$week_end +1, sep = "-")

  flushing_timing <- flushing_timing %>%
    select(species_full,
           id,
           flushing_date)
  flushing_timing$flushing_value <- 1

  data_flush <- merge(timelines_flush, flushing_timing, by = c("species_full","id","flushing_date"), all.x = TRUE)
  data_flush$flushing_value <- ifelse(is.na(data_flush$flushing_value),"0", data_flush$flushing_value)
  data_flush$flushing_value <- ifelse(is.na(data_flush$value),NA, data_flush$flushing_value)
  data_flush$flushing_value <- as.numeric(data_flush$flushing_value)

  data_flush <- data_flush %>%
    group_by(species_full, date) %>%
    dplyr::summarise(mean_week = mean(flushing_value, na.rm = TRUE),
                     count_week = sum(flushing_value),
                     total_week = length(flushing_value))

  # # filter out species with too few total siteyears to have a meaningfull average year
  # data_flush <- data_flush %>%
  #   filter(total_week >= minimum_siteyears)

  # data_flush$mean_week <- ifelse(data_flush$mean_week < minimum_event_frequency, 0, data_flush$mean_week)

  data_flush_site <- inner_join(data_flush, census_site, by = c("species_full"))


  final_flush_site <- data_flush_site %>%
    group_by(date) %>%
    dplyr::summarise(ss_flush = sum(mean_week*basal_area_site, na.rm = TRUE)/total_basal_area_site *100) #,
  #-----------------------------------------------------------------------


  #-----------------------------------------------------------------------
  #------------ plots -------  -------------------------------------------
  #-----------------------------------------------------------------------
  combined <- merge(final_turn_site, final_sen_site, by = c("date"), all.x = TRUE)
  combined <- merge(combined, final_flush_site, by = c("date"), all.x = TRUE)


  return(combined)
}
