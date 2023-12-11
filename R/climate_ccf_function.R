#' cross-correlations species pheno data with climate data
#'
#' @param data junglerhythms data file
#' @param climate yangambi climate data
#' @param climate.variable precip, sun or temp (one variable at a time)
#' @param species_name list of species
#' @param pheno only one phenophase
#' @export
#' @return dataframe significant correlation values and timing



climate_ccf <- function(
  data = data,
  climate = climate,
  climate.variable = "precip",
  species_name = "Scorodophoeus zenkeri",
  pheno = "leaf_turnover"){

  ccf_output_clim <- data.frame(matrix(0,nrow=length(species_list),ncol=17))
  colnames(ccf_output_clim)[1] <- "species_full"
  colnames(ccf_output_clim)[2:14] <- c("t-6","t-5","t-4","t-3","t-2","t-1","t","t+1","t+2","t+3","t+4","t+5","t+6")
  colnames(ccf_output_clim)[15] <- "ci"
  colnames(ccf_output_clim)[16] <- paste("corr", climate.variable, sep = "_")
  colnames(ccf_output_clim)[17] <- paste("corr", climate.variable, "timing", sep = "_")

  for (j in 1:length(species_name)){
    data_subset <- data %>%
      filter(species_full %in% species_name[j]) %>%
      filter(phenophase %in% pheno)


    # data_subset <- timelines_sp_consec_turn %>%
    #   filter(species_full %in% "Afrostyrax lepidophyllus")

    # convert date
    data_subset$date_monthly <- format(data_subset$date, "%Y-%m")
    data_subset$year <- format(data_subset$date, "%Y")

    # average by year_month
    data_sp_monthly <- data_subset %>%
      group_by(total_nr_events,
               year, date_monthly) %>%
      dplyr::summarise(monthly_value = mean(scaled_value)) %>%
      dplyr::arrange(date_monthly)

    # depending on species timeline, change range climate timeline
    first_year <- as.numeric(data_sp_monthly$year[1])
    last_year <- as.numeric(last(data_sp_monthly$year))
    climate_range <- climate %>%
      filter(year >= first_year,
             year <= last_year) %>%
      dplyr::arrange(date)


    # merge with climate data
    data_sp_monthly <- merge(data_sp_monthly, climate_range, by = c("date_monthly","year"), all.x = TRUE)

    # select climate variable needed
    if(climate.variable == "precip"){
      climate_ts <- data_sp_monthly$precip
    } else if(climate.variable == "sun"){
      climate_ts <- data_sp_monthly$insol_JR
    } else if(climate.variable == "temp"){
      climate_ts <- data_sp_monthly$tmax_JR
    }

    # nr of events
    nr_events <- max(data_sp_monthly$total_nr_events)

    ci <- 0.95
    # cross-correlation + confidence interval
    if(max(data_sp_monthly$monthly_value, na.rm = TRUE) > 0 & nr_events >= 5){
      corr.clim <- ccf(climate_ts, data_sp_monthly$monthly_value, lag = 6, pl = FALSE, na.action = na.pass)
      # plot(corr.clim, main = paste(species_name[j], pheno, climate.variable, sep = "-"))
      ci_value <- qnorm((1 + ci)/2)/sqrt(corr.clim$n.used)
      ccf_output_clim[j,1] <- species_name[j]
      ccf_output_clim[j,2:14] <- as.numeric(corr.clim$acf)
      ccf_output_clim[j,15] <- as.numeric(ci_value)
    } else {
      ccf_output_clim[j,1] <- species_name[j]
      ccf_output_clim[j,2:14] <- NA
      ccf_output_clim[j,15] <- NA
    }


  }

  ccf_output_clim$corr <- ifelse(abs(ccf_output_clim[,8]) > ccf_output_clim$ci, ccf_output_clim[,8], # 0
                                 ifelse(abs(ccf_output_clim[,7]) > ccf_output_clim$ci, ccf_output_clim[,7], # t-1
                                        ifelse(abs(ccf_output_clim[,6]) > ccf_output_clim$ci, ccf_output_clim[,6], # t-2
                                               ifelse(abs(ccf_output_clim[,5]) > ccf_output_clim$ci, ccf_output_clim[,5], # t-3
                                                      NA))))
  ccf_output_clim$corr.timing <- ifelse(abs(ccf_output_clim[,8]) > ccf_output_clim$ci, "0", # 0
                                        ifelse(abs(ccf_output_clim[,7]) > ccf_output_clim$ci, "-1", # t-1
                                               ifelse(abs(ccf_output_clim[,6]) > ccf_output_clim$ci, "-2", # t-2
                                                      ifelse(abs(ccf_output_clim[,5]) > ccf_output_clim$ci, "-3", # t-3
                                                             NA))))

  ccf_output_clim <- ccf_output_clim %>%
    dplyr::select(species_full,
                  corr,
                  corr.timing)

  colnames(ccf_output_clim)[2] <- paste("corr", climate.variable, sep = "_")
  colnames(ccf_output_clim)[3] <- paste("corr", climate.variable, "timing", sep = "_")
  ccf_output_clim$phenophase <- pheno

  return(ccf_output_clim)

}
