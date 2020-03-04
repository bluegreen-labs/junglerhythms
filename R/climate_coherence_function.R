#' coherence analysis (co-fourier analysis) with climate data
#'
#' @param data timelines of consecutive years, already at species level
#' @param climate
#' @param species_name
#' @param pheno
#' @export
#' @return phase difference between climate and phenology time series, coherence values and significance


spans_lookup_monthly <- list(observations = c(24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 180, 192, 204, 216, 228, 240),
                             # spans_smooth = list(NULL, NULL, NULL, NULL, # smoothed spectrum of data -> smoothed periodogram bandwith to approx. 0.05
                             #                     NULL, NULL, NULL,
                             #                     3,3,3,3,
                             #                     3,3,3,
                             #                     c(3,3), c(3,3), c(3,3), c(3,3), c(3,3)),
                             spans_smooth=list(NULL, NULL, NULL, NULL, # smoothed spectrum of data -> smoothed periodogram bandwith to approx. 0.1
                                               3,3,3,
                                               c(3,3), c(3,3), c(3,3), c(3,3),
                                               c(3,5), c(3,5), c(3,5),
                                               c(5,5), c(5,5), c(5,5),
                                               c(5,7), c(5,7)),
                             spans_super_smooth = list(c(5,7), c(9,9), c(11,11), c(13,13),
                                                       c(15,17), c(19,21), c(21,21), c(23,23),
                                                       c(25,27), c(29,29), c(29,31), c(33,35),
                                                       c(37,39), c(37,39), c(39,41), c(45,45),
                                                       c(45,45), c(49,51), c(49,51)))
# spectrum function
spec_cross_fun <- function(x)spec.pgram(x, spans = spans_lookup_monthly$spans_smooth[[which.min(abs(spans_lookup_monthly$observations-length(x)))]],
                                        demean = T, detrend = T, plot = F)

#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------

climate_coherence <- function(
  data = data,
  climate = climate,
  species_name = "Scorodophoeus zenkeri",
  pheno = "leaf_turnover"){

  coherence_output <- data.frame()

  for (j in 1:length(species_name)){

    data_subset <- data %>%
      filter(species_full %in% species_name[j]) %>%
      filter(phenophase %in% pheno)

# data_subset <- timelines_sp_consec %>%
#   filter(species_full %in% "Afzelia bipindensis")

    # convert date
    data_subset$date_monthly <- format(data_subset$date, "%Y-%m")

    # average by year_month
    data_sp_monthly <- data_subset %>%
      group_by(year, date_monthly) %>%
      dplyr::summarise(monthly_value = mean(scaled_value)) %>%
      dplyr::arrange(date_monthly)

    # data as timeseries, species_level
    first_year <- as.numeric(data_sp_monthly$year[1])
    last_year <- as.numeric(last(data_sp_monthly$year))
    data_ts <- ts(data = data_sp_monthly$monthly_value, start = first_year, frequency = 12)

    # depending on species timeline, change range climate timeseries
    climate_range <- climate %>%
      filter(year >= first_year,
             year <= last_year) %>%
      dplyr::arrange(date)
    # climate time series
    climate_insol_ts <- ts(data = climate_range$insol_JR, start = first_year, frequency = 12)
    climate_tmax_ts <- ts(data = climate_range$tmax_JR, start = first_year, frequency = 12)
    climate_precip_ts <- ts(data = climate_range$precip, start = first_year, frequency = 12)
    climate_precip_neg_ts <- ts(data = climate_range$precip*(-1), start = first_year, frequency = 12)

    ##########################
    # coherence analysis
    ##########################

    #---- precip -----------------------------------------------------------------------------------
    ts_cross_precip <- ts.intersect(data_ts, climate_precip_ts)
    # the dominant frequency of phen: spec[,1] is from the phen data
    freq_phen <- spec_cross_fun(ts_cross_precip)$freq[which.max(spec_cross_fun(ts_cross_precip)$spec[,1])]/12
    # dominant cycle in months
    cycl_phen <- 1/freq_phen
    # the phase difference (in radians) at the dominant frequency
    phase_radians_precip <- spec_cross_fun(ts_cross_precip)$phase[which.max(spec_cross_fun(ts_cross_precip)$spec[,1])]
    phase_degree_precip <- (phase_radians_precip * 180) / (pi) # transform to degrees
    phase_diff_precip <- phase_degree_precip/360*cycl_phen # transform to months
    # the coherence at the dominant frequency
    coh_precip <- spec_cross_fun(ts_cross_precip)$coh[which.max(spec_cross_fun(ts_cross_precip)$spec[,1])]
    # significance
    sr <- spec_cross_fun(ts_cross_precip)
    f = qf(.95, 2, sr$df-2)
    L = (sr$df * sr$n.used)/(sr$orig.n * 2)
    C = f/(L-1+f) # coherence large then this is significant (p < 0.05)
    signif <- ifelse(coh_precip > C, TRUE, FALSE)
    # plot(sr, plot.type = "coh", ci.lty = 2)
    # abline(h = C)

    # #---- precip inverse test -----------------------------------------------------------------------------------
    # ts_cross_precip_inv <- ts.intersect(climate_precip_ts, data_ts)
    # # the dominant frequency of the climate variable: spec[,1] is from the climate timeseries
    # freq_clim <- spec_cross_fun(ts_cross_precip_inv)$freq[which.max(spec_cross_fun(ts_cross_precip_inv)$spec[,1])]/12
    # # dominant cycle in months
    # cycl_clim <- 1/freq_clim
    # # the phase difference (in radians) at the dominant frequency
    # phase_radians_precip_inv <- spec_cross_fun(ts_cross_precip_inv)$phase[which.max(spec_cross_fun(ts_cross_precip_inv)$spec[,1])]
    # phase_degree_precip_inv <- (phase_radians_precip_inv * 180) / (pi) # transform to degrees
    # phase_diff_precip_inv <- phase_degree_precip_inv/360*cycl_clim # transform to months
    # # the coherence at the dominant frequency
    # coh_precip_inv <- spec_cross_fun(ts_cross_precip_inv)$coh[which.max(spec_cross_fun(ts_cross_precip_inv)$spec[,1])]
    # # significance
    # sr_inv <- spec_cross_fun(ts_cross_precip_inv)
    # f_inv = qf(.95, 2, sr_inv$df-2)
    # L_inv = (sr_inv$df * sr_inv$n.used)/(sr_inv$orig.n * 2)
    # C_inv = f_inv/(L_inv-1+f_inv) # coherence large then this is significant (p < 0.05)
    # signif_inv <- ifelse(coh_precip_inv > C_inv, TRUE, FALSE)

    #---- precip negative -----------------------------------------------------------------------------------
    ts_cross_precip_neg <- ts.intersect(data_ts, climate_precip_neg_ts)
    # the dominant frequency of the climate variable: spec[,1] is from the climate timeseries
    freq_phen_neg <- spec_cross_fun(ts_cross_precip_neg)$freq[which.max(spec_cross_fun(ts_cross_precip_neg)$spec[,1])]/12
    # dominant cycle in months
    cycl_phen_neg <- 1/freq_phen_neg
    # the phase difference (in radians) at the dominant frequency
    phase_radians_precip_neg <- spec_cross_fun(ts_cross_precip_neg)$phase[which.max(spec_cross_fun(ts_cross_precip_neg)$spec[,1])]
    phase_degree_precip_neg <- (phase_radians_precip_neg * 180) / (pi) # transform to degrees
    phase_diff_precip_neg <- phase_degree_precip_neg/360*cycl_phen_neg # transform to months
    # the coherence at the dominant frequency
    coh_precip_neg <- spec_cross_fun(ts_cross_precip_neg)$coh[which.max(spec_cross_fun(ts_cross_precip_neg)$spec[,1])]
    # significance
    sr_neg <- spec_cross_fun(ts_cross_precip_neg)
    f_neg = qf(.95, 2, sr_neg$df-2)
    L_neg = (sr_neg$df * sr_neg$n.used)/(sr_neg$orig.n * 2)
    C_neg = f_neg/(L_neg-1+f_neg) # coherence large then this is significant (p < 0.05)
    signif_neg <- ifelse(coh_precip_neg > C_neg, TRUE, FALSE)


    coherence_output_species <- data.frame(species = species_name[j],
      phenophase = pheno,
      cycl_phen = cycl_phen,
      phase_diff_precip = phase_diff_precip,
      coh_precip,
      signif = signif,
      # cycl_clim = cycl_clim,
      # phase_diff_precip_inv = phase_diff_precip_inv,
      # coh_precip_inv,
      # signif_inv = signif_inv,
      cycl_phen_neg = cycl_phen_neg,
      phase_diff_precip_neg = phase_diff_precip_neg,
      coh_precip_neg,
      signif_neg = signif_neg
    )
    coherence_output <- rbind(coherence_output, coherence_output_species)

  }

  return(coherence_output)
}






# ##########################
# # find cff max lag = cross correlation
# ##########################
# # merge with climate data
# data_sp_monthly <- merge(data_sp_monthly, climate_range, by = "date_monthly", all.x = TRUE)
# # cross-correlation + confidence interval
# corr.precip <- ccf(data_sp_monthly$precip, data_sp_monthly$monthly_value, lag = 6, pl = FALSE)
# plot(corr.precip)
# corr.insol_JR <- ccf(data_sp_monthly$insol_JR, data_sp_monthly$monthly_value, lag = 6, pl = FALSE)
# plot(corr.insol_JR)
# ##########################


#     # #-------------------------------------------------
#     # # coherence <- ccoh(data_ts, climate_insol_ts, channel = c(1,1), wl = 24, ovlp = 1)
#     # library(seewave)
#     # coherence <- coh(data_ts, climate_insol_ts, f = 12) #, xlim = c(0,12/1000))
#     # lags      <- coherence[,1]
#     # vals      <- coherence[,2]
#     # val.max   <- vals[which.max(vals)]
#     # lag.max   <- lags[which.max(vals)]
#     #
#     # freq_dom = (spec_fun(data_ts)$freq[which.max(spec_fun(data_ts)$spec)])/12 #frequency of the dominant peak
#     # cycle_dom = 1/freq_dom
#
#
#     ##########################
#     # find cff max lag = cross correlation
#     ##########################
#
#     # merge with climate data
#     data_grow_monthly <- merge(data_grow_monthly, climate, by = "yr_month", all.x = TRUE)
#     plot(data_grow_monthly$monthly_value)
#     plot(data_grow_monthly$precip)
#     # cross-correlation + confidence interval
#     corr.precip <- ccf(data_grow_monthly$precip, data_grow_monthly$monthly_value, lag = 6, pl = FALSE)
#     plot(corr.precip)
#     corr.insol_JR <- ccf(data_grow_monthly$insol_JR, data_grow_monthly$monthly_value, lag = 6, pl = FALSE)
#     plot(corr.insol_JR)
#     ##########################
#
#
#     ##########################
#     # fourier transform
#     ##########################
#     # Function for co-Fourier analysis of a simulated cosine curve against the empirical data
#     spec_cross_fun <- function(x)spec.pgram(x, spans = c(3,3), #spans=spans_lookup$spans_smooth[[which.min(abs(spans_lookup$months-length(x)))]],
#                                             demean=T,detrend=T,plot=F)
#
# # #Function to extract data from co-Fourier anaysis regarding phase difference in radians
# # # ts_cross_insol <- ts.intersect(data_ts, climate_insol_ts)
# # # ts_cross_tmax <- ts.intersect(data_ts, climate_tmax_ts)
# # ts_cross_precip <- ts.intersect(data_ts, climate_precip_ts)
# #
# # #the phase difference (in radians) at the dominant frequency
# # # phase_radians_insol <- spec_fun(ts_cross_insol)$phase[which.max(spec_fun(ts_cross_insol)$spec[,1])]
# # # coh_insol <- spec_fun(ts_cross_insol)$coh[which.max(spec_fun(ts_cross_insol)$spec[,1])]
# # # phase_degree_insol <- (phase_radians_insol * 180) / (pi)
# #
# #
# # freq_phen <- (spec_fun(data_ts)$freq[which.max(spec_fun(data_ts)$spec)])/12
# # cycl_phen <- 1/freq_phen
# # # dominant cycle of first timeseries
# # freq_phen1 <- spec_fun(ts_cross_precip)$freq[which.max(spec_fun(ts_cross_precip)$spec[,1])]/12
# # cycl_phen1 <- 1/freq_phen1
# # # dominant cycle of second timeseries
# # freq_precip <- spec_fun(ts_cross_precip)$freq[which.max(spec_fun(ts_cross_precip)$spec[,2])]/12
# # cycl_precip <- 1/freq_precip
# #
# # #the phase difference (in radians) at the dominant frequency of phen
# # phase_radians_precip <- spec_fun(ts_cross_precip)$phase[which.max(spec_fun(ts_cross_precip)$spec[,1])]
# # coh_precip <- spec_fun(ts_cross_precip)$coh[which.max(spec_fun(ts_cross_precip)$spec[,1])]
# # phase_degree_precip <- (phase_radians_precip * 180) / (pi)
# # #the phase difference (in radians) at the dominant frequency of precip
# # phase_radians_precip2 <- spec_fun(ts_cross_precip)$phase[which.max(spec_fun(ts_cross_precip)$spec[,2])]
# # coh_precip2 <- spec_fun(ts_cross_precip)$coh[which.max(spec_fun(ts_cross_precip)$spec[,2])]
# # phase_degree_precip2 <- (phase_radians_precip2 * 180) / (pi)
# #
# #
# # # reverse order intersect
# # ts_cross_precip_reverse <- ts.intersect(climate_precip_ts, data_ts)
# # # dominant cycle of first timeseries
# # freq_phen1 <- spec_fun(ts_cross_precip_reverse)$freq[which.max(spec_fun(ts_cross_precip_reverse)$spec[,1])]/12
# # cycl_phen1 <- 1/freq_phen1
# # # dominant cycle of second timeseries
# # freq_precip <- spec_fun(ts_cross_precip_reverse)$freq[which.max(spec_fun(ts_cross_precip_reverse)$spec[,2])]/12
# # cycl_precip <- 1/freq_precip
# # #the phase difference (in radians) at the dominant frequency of precip
# # phase_radians_precip <- spec_fun(ts_cross_precip_reverse)$phase[which.max(spec_fun(ts_cross_precip_reverse)$spec[,1])]
# # coh_precip <- spec_fun(ts_cross_precip_reverse)$coh[which.max(spec_fun(ts_cross_precip_reverse)$spec[,1])]
# # phase_degree_precip <- (phase_radians_precip * 180) / (pi)
# # #the phase difference (in radians) at the dominant frequency of phen
# # phase_radians_precip2 <- spec_fun(ts_cross_precip_reverse)$phase[which.max(spec_fun(ts_cross_precip_reverse)$spec[,2])]
# # coh_precip2 <- spec_fun(ts_cross_precip_reverse)$coh[which.max(spec_fun(ts_cross_precip_reverse)$spec[,2])]
# # phase_degree_precip2 <- (phase_radians_precip2 * 180) / (pi)
# # ############ phase difference the same with reversed signs
# # #----------> coherence is 1, because periodograms are not smoothed! You can only do cohorence on smoothed
# #
# #  #---------
# #  sr=spec.pgram(ts_cross_precip, spans=rep(3,5))#kernel("daniell",9),taper=0,plot=FALSE) #cbind(soi,rec)
# #  sr$df
# #  f = qf(.999, 2, sr$df-2)
# #  C = f/(18+f)
# #  plot(sr, plot.type = "coh", ci.lty = 2)
# #  plot(sr, plot.type = "mar", ci.lty = 2)
# #  abline(h = C)
# #
# #  max(sr$coh)
# #  which.max(sr$coh)
# #  sr$phase[36]
# #  which.max(sr$spec)
# #  sr$phase[177]
# #
# #  plot(sr$freq, sr$coh)
# #
# #  #-----------
# #   # https://www.rdocumentation.org/packages/IRISSeismic/versions/1.5.2/topics/crossSpectrum
# #  library(IRISSeismic)
# #  # Calculate the cross spectrum
# #  DF <- crossSpectrum(ts.union(data_ts, climate_insol_ts),spans=c(3,5,7,9))
# #
# #  max(DF$coh)
# #  which.max(DF$coh)
# #  DF$phase[43]
# #  which.max(DF$spec)
# #
# #  # Calculate the transfer function
# #  transferFunction <- DF$Pxy / DF$Pxx
# #  transferAmp <- Mod(transferFunction)
# #  transferPhase <- pracma::mod(Arg(transferFunction) * 180/pi,360)
# #
# #  # 2 rows
# #  layout(matrix(seq(2)))
# #
# #  # Plot
# #  plot(1/DF$freq,transferAmp,type='l',log='x',
# #       xlab="Period (sec)",
# #       main="Transfer Function Amplitude")
# #
# #  plot(1/DF$freq,transferPhase,type='l',log='x',
# #       xlab="Period (sec)", ylab="degrees",
# #       main="Transfer Function Phase")
# #
# #  # #-----------
#  # test <- spectrum(ts.union(data_ts, climate_insol_ts),spans=c(3,5,7,9))
#  # max(test$coh)
#  # which.max(test$coh)
#  # test$phase[43]
#  # which.max(test$spec)
#  #
#  #
#  # #----------------
#  # #https://stats.stackexchange.com/questions/21020/how-can-i-estimate-the-phase-difference-between-two-periodic-time-series
#  #  datos <- ts.union(climate_precip_ts,data_ts)
#  # datos <- window(datos,
#  #                 start=c(1938,1),
#  #                 end=c(1959,12))
#  # sp <- spectrum(datos,
#  #                main="Petróleo e IPC",
#  #                spans=rep(3,5))
#  # par(mfrow=c(2,1))
#  # plot(sp,plot.type="coh")
#  # plot(sp,plot.type="phase")
#
#
#
# # fourier transform
# fourier_precip <- data.frame(freq = spec_fun(climate_precip_ts)$freq/12,
#                           spec = spec_fun(climate_precip_ts)$spec,
#                           spec_norm = spec_fun(climate_precip_ts)$spec*(1/mean(spec_fun(climate_precip_ts)$spec)),
#                           freq_null = spec_null_fun(climate_precip_ts)$freq/12,
#                           spec_null = spec_null_fun(climate_precip_ts)$spec,
#                           spec_norm_null = spec_null_fun(climate_precip_ts)$spec*(1/mean(spec_null_fun(climate_precip_ts)$spec)))
# a <- ggplot(fourier_precip) +
#    geom_line(aes(x = freq,
#                  y = spec_norm)) +
#    geom_line(aes(x = freq_null,
#                  y = spec_norm_null),
#              color = "blue") +
#    scale_y_continuous(labels=NULL) +
#    ylab("Power (standardised spectrum)") +
#    xlab("Cycle frequency") +
#    ggtitle("") +
#    theme_light(base_size=14)+
#    theme(legend.position = "none")
#
# # fourier transform
# fourier_phen <- data.frame(freq = spec_fun(data_ts)$freq/12,
#                              spec = spec_fun(data_ts)$spec,
#                              spec_norm = spec_fun(data_ts)$spec*(1/mean(spec_fun(data_ts)$spec)),
#                              freq_null = spec_null_fun(data_ts)$freq/12,
#                              spec_null = spec_null_fun(data_ts)$spec,
#                              spec_norm_null = spec_null_fun(data_ts)$spec*(1/mean(spec_null_fun(data_ts)$spec)))
# b <- ggplot(fourier_phen) +
#   geom_line(aes(x = freq,
#                 y = spec_norm)) +
#   geom_line(aes(x = freq_null,
#                 y = spec_norm_null),
#             color = "blue") +
#   scale_y_continuous(labels=NULL) +
#   ylab("Power (standardised spectrum)") +
#   xlab("Cycle frequency") +
#   ggtitle("") +
#   theme_light(base_size=14)+
#   theme(legend.position = "none")
# grid.arrange(a, b, heights = c(1,1))
#
#
#
