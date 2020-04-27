standlevel_phen <- function(
  data = data,
  species_list_dorm = dorm_sp1,
  species_list_turn = turn_sp1
  ){

  combined <- data.frame()

  timelines_dorm <- two_year_gaps(data = data,
                                  species_name = species_list_dorm,
                                  pheno = "leaf_dormancy")
  timelines_turn <- two_year_gaps(data = data,
                                  species_name = species_list_turn,
                                  pheno = "leaf_turnover")

  #-----------------------------------------------------------------------
  #------------ Leaf turnover  -------------------------------------------
  #-----------------------------------------------------------------------
  timelines_turn <- timelines_turn %>%
    filter(!is.na(value))

  data_LT <- timelines_turn %>%
    group_by(species_full, week) %>%
    dplyr::summarise(mean_week = mean(value),
              count_week = sum(value),
              total_week = length(value))

  # filter out species with too few total siteyears across individuals to have a meaningfull average year
  data_LT <- data_LT %>%
    filter(total_week >= minimum_siteyears)

  data_LT_site <- inner_join(data_LT, census_site, by = c("species_full"))
  data_LT_plot <- inner_join(data_LT, census_plot, by = c("species_full"))

  # # check to see if total basal area per week is constant
  # test <- data_LT_plot %>%
  #   group_by(Plot, week) %>%
  #   summarise(total_basal_area_week = sum(basal_area_plot)) #,

  # ### reference website for confidence intervals
  # ### https://stats.stackexchange.com/questions/224029/calculating-weighted-cis-and-interpretation
  final_LT_site <- data_LT_site %>%
    group_by(week) %>%
    dplyr::summarise(ss = sum(mean_week*basal_area_site, na.rm = TRUE)/total_basal_area_site *100) #,
  # ss2 = wtd.mean(mean_week, basal_area_site, na.rm = TRUE)*100,
  # weighted_sd = sqrt(wtd.var(mean_week, basal_area_site, na.rm = TRUE))*100,
  # weighted_CI = ss2 + 1.96*weighted_sd)

  final_LT_plot <- data_LT_plot %>%
    group_by(Plot, week) %>%
    dplyr::summarise(ss = sum(mean_week*basal_area_plot, na.rm = TRUE))
  final_LT_plot <- inner_join(final_LT_plot, total_basal_area_plot, by = c("Plot"))
  final_LT_plot <- final_LT_plot %>%
    mutate(ss = ss/total_basal_area_plot*100)
  #-----------------------------------------------------------------------

  #-----------------------------------------------------------------------
  #------------ Leaf dormancy  -------------------------------------------
  #-----------------------------------------------------------------------
  timelines_dorm <- timelines_dorm %>%
    filter(!is.na(value))

  data_LD <- timelines_dorm %>%
    group_by(species_full, week) %>%
    dplyr::summarise(mean_week = mean(value),
              count_week = sum(value),
              total_week = length(value))

  data_LD <- data_LD %>%
    filter(total_week >= minimum_siteyears)


  data_LD_site <- inner_join(data_LD, census_site, by = c("species_full"))
  data_LD_plot <- inner_join(data_LD, census_plot, by = c("species_full"))

  final_LD_site <- data_LD_site %>%
    group_by(week) %>%
    dplyr::summarise(ss = sum(mean_week*basal_area_site, na.rm = TRUE)/total_basal_area_site *100) #,

  final_LD_plot <- data_LD_plot %>%
    group_by(Plot, week) %>%
    dplyr::summarise(ss = sum(mean_week*basal_area_plot, na.rm = TRUE))
  final_LD_plot <- inner_join(final_LD_plot, total_basal_area_plot, by = c("Plot"))
  final_LD_plot <- final_LD_plot %>%
    mutate(ss = ss/total_basal_area_plot*100)
  #-----------------------------------------------------------------------


  #-----------------------------------------------------------------------
  #------------ Leaf flushing  -------------------------------------------
  #-----------------------------------------------------------------------
  timelines_flush <- two_year_gaps(data = data,
                                   species_name = species_list_dorm,
                                   pheno = "leaf_dormancy")


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
    group_by(species_full, week) %>%
    dplyr::summarise(mean_week = mean(flushing_value, na.rm = TRUE),
              count_week = sum(flushing_value),
              total_week = length(flushing_value))

  # filter out species with too few total siteyears to have a meaningfull average year
  data_flush <- data_flush %>%
    filter(total_week >= minimum_siteyears)

  # data_flush$mean_week <- ifelse(data_flush$mean_week < minimum_event_frequency, 0, data_flush$mean_week)

  data_flush_site <- inner_join(data_flush, census_site, by = c("species_full"))
  data_flush_plot <- inner_join(data_flush, census_plot, by = c("species_full"))


  final_flush_site <- data_flush_site %>%
    group_by(week) %>%
    dplyr::summarise(ss = sum(mean_week*basal_area_site, na.rm = TRUE)/total_basal_area_site *100) #,

  final_flush_plot <- data_flush_plot %>%
    group_by(Plot, week) %>%
    dplyr::summarise(ss = sum(mean_week*basal_area_plot, na.rm = TRUE))
  final_flush_plot <- inner_join(final_flush_plot, total_basal_area_plot, by = c("Plot"))
  final_flush_plot <- final_flush_plot %>%
    mutate(ss = ss/total_basal_area_plot*100)
  #-----------------------------------------------------------------------


  #-----------------------------------------------------------------------
  #------------ plots -------  -------------------------------------------
  #-----------------------------------------------------------------------
  combined <- final_LT_site
  colnames(combined)[2] <- "ss_turn"
  combined <- merge(combined, final_LD_site, by = c("week"), all.x = TRUE)
  combined$ss_senescence <- combined$ss_turn + combined$ss

  colnames(final_flush_site)[2] <- "ss_flush"
  combined <- merge(combined, final_flush_site, by = c("week"), all.x = TRUE)


  return(combined)
}
