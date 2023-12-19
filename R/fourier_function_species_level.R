#' Peak detection and cyclicity using fourier transform
#'
#' @param data junglerhythms data file
#' @param species_name list of species
#' @param pheno only one phenophase
#' @param perio_plot TRUE/FALSE plot periodogram
#' @export
#' @return timing and stats for peaks

peak_detection <- function(
  data = data,
  species_name = "Scorodophoeus zenkeri",
  pheno = "leaf_turnover",
  perio_plot = TRUE){

  # following Bush et al. 2017
  # Bush, E. R., Abernethy, K. A., Jeffery, K., Tutin, C., White, L., Dimoto, E., ...,
  # Bunnefeld, N. (2017). Fourier analysis to detect phenological cycles using tropical
  # field data and simulations. Methods in Ecology and Evolution, 8, 530-540.
  # -> user-specified widths of the Daniell kernel smoother
  # -> spectral estimate is smoothed to the specified band-width
  # -> this depends on the lenght of the time-series
  # -> for various ts lenghts (with observation frequenties of 48 per year) the span is provided
  spans_lookup <- list(observations = c(96,144, 192, 240, 288, 336, 384, 432, 480, 528,576,624,672,720,768,816,864,912,960, 1008),
                       # smoothed spectrum of data -> smoothed periodogram bandwith to approx. 0.1
                       spans_smooth = list(NULL, NULL, NULL, NULL,
                                           3,3,3,3,
                                           c(3,3), c(3,3), c(3,3),
                                           c(3,5), c(3,5), c(3,5),
                                           c(5,5), c(5,5), c(5,5),
                                           c(5,7), c(5,7), c(5,7)),
                       # super-smoothed spectrum of data -> smoothed periodogram bandwith to approx. 1
                       spans_super_smooth = list(c(5,7), c(7,9), c(11,11), c(13,13),
                                                 c(15,17), c(19,19), c(19,21), c(23,23),
                                                 c(25,27), c(27,29), c(29,31),
                                                 c(33,33), c(35,37), c(37,39),
                                                 c(39,41), c(45,45), c(45,47),
                                                 c(49,51), c(49,51), c(49,51)))

  # spectrum function for general smoother periodogram, to identify dominant peaks
  spec_fun <- function(x) spectrum(x, spans = spans_lookup$spans_smooth[[which.min(abs(spans_lookup$observations - length(x)))]],
                                   plot = F, demean = T, detrend = T)

  # spectrum function for null hypothesis spectrum (super-smoothed spectrum of data -> smoothed periodogram bandwith to approx. 1)
  spec_null_fun <- function(x) spectrum(x, spans = spans_lookup$spans_super_smooth[[which.min(abs(spans_lookup$observations - length(x)))]],
                                        plot = F, demean = T, detrend = T)

  spectrum_output <- data.frame()

  for (j in 1:length(species_name)){

    data_subset <- data %>%
      filter(species_full %in% species_name[j]) %>%
      filter(phenophase %in% pheno)

    #----------------------------------------------------------------------
    #----- fourier transform + periodogram
    #----------------------------------------------------------------------
    # time series at the species level
    first_year <- as.numeric(format.Date(data_subset$date[1], "%Y"))
    data_ts <- ts(data = data_subset$mean_value, start = c(first_year,1), frequency = 48)

    # fourier transform -> spectrum
    fourier_df <- data.frame(freq = spec_fun(data_ts)$freq/48,
                             spec = spec_fun(data_ts)$spec,
                             spec_norm = spec_fun(data_ts)$spec*(1/mean(spec_fun(data_ts)$spec)),
                             freq_null = spec_null_fun(data_ts)$freq/48,
                             spec_null = spec_null_fun(data_ts)$spec,
                             spec_norm_null = spec_null_fun(data_ts)$spec*(1/mean(spec_null_fun(data_ts)$spec)))

    # # extract key variables from spectrum outputs
    #----- Periodogram analysis with confidence intervals to find frequency and phase
    freq_dom = (spec_fun(data_ts)$freq[which.max(spec_fun(data_ts)$spec)])/48 # frequency of the dominant peak
    cycle_dom = 1/freq_dom
    spec_dom = max(spec_fun(data_ts)$spec) # spectrum of dominant peak
    dfree = spec_fun(data_ts)$df # degrees of freedom of the smoothed periodogram

    lower_ci_95 = (dfree*spec_dom)/(qchisq(c(0.975),dfree)) # lower CI of dominant peak of smoothed periodogram
    spec_null_dom <- spec_null_fun(data_ts)$spec[which(abs((spec_null_fun(data_ts)$freq)/48-freq_dom)==min(abs((spec_null_fun(data_ts)$freq)/48-freq_dom)))] # the corresponding dominant peak from the null continuum
    sig_95 = lower_ci_95 > spec_null_dom # If the null spectrum is lower than the lower confidence interval, dominant peak can be considered "signficantly different" to null hypothesis
    cycle_category <- as.factor(ifelse(1/freq_dom > length(data_ts)/2,"No cyclicity","cyclicity"))
    # Is the dominant cycle both significant and represent a regular cycle?
    sig_cycle_95 <- as.factor(ifelse(cycle_category!="No cyclicity" & sig_95==TRUE,TRUE,FALSE))

    lower_ci_99 = (dfree*spec_dom)/(qchisq(c(0.995),dfree)) # lower CI of dominant peak of smoothed periodogram
    sig_99 = lower_ci_99 > spec_null_dom # If the null spectrum is lower than the lower confidence interval, dominant peak can be considered "signficantly different" to null hypothesis
    # Is the dominant cycle both significant and represent a regular cycle?
    sig_cycle_99 <- as.factor(ifelse(cycle_category!="No cyclicity" & sig_99==TRUE,TRUE,FALSE))

    spec_zero = spec_fun(data_ts)$spec[1]
    spec_ratio = spec_zero/spec_dom


    spectrum_output_species <- data.frame(species_full = species_name[j],
      phenophase = pheno,
      freq_dom = freq_dom,
      cycle_dom = cycle_dom,
      spec_dom = spec_dom,
      spec_zero = spec_zero,
      spec_ratio = spec_ratio,
      dfree = dfree,
      lower_ci_95 = lower_ci_95,
      spec_null_dom = spec_null_dom,
      sig_95 = sig_95,
      cycle_category = cycle_category,
      sig_cycle_95 = sig_cycle_95,

      lower_ci_99 = lower_ci_99,
      sig_99 = sig_99,
      sig_cycle_99 = sig_cycle_99)

    spectrum_output <- rbind(spectrum_output, spectrum_output_species)


    #-----------------------------
    # Plot overview plots and periodograms for visual inspection
    #-----------------------------
    # timeseries plot
    timeline_plot <- ggplot(data_subset,
                             aes(x = date,
                                 y = mean_value))+
      geom_line()+
      scale_y_continuous(labels=NULL)+
      ylab("Observations") +
      xlab("") +
      ggtitle(paste(species_name[j], "-", pheno)) +
      theme_light(base_size = 14) +
      theme(legend.position = "none")


    ##Plot periodograms for visual inspection
    spectrum_output_species$xloc <- 0.4
    spectrum_output_species$yloc <- max(fourier_df$spec_norm, na.rm=T)*0.75
    spectrum_output_species$yloc2 <- max(fourier_df$spec_norm, na.rm=T)*0.5
    spectrum_output_species$cycle_months <- ifelse(spectrum_output_species$sig_95 %in% "TRUE",
                                                   format(round(spectrum_output_species$cycle_dom/4,1)), NA)

    periodogram <- ggplot(fourier_df) +
      geom_line(aes(x = freq,
                    y = spec_norm)) +
      geom_line(aes(x = freq_null,
                    y = spec_norm_null),
                color = "blue") +
      scale_y_continuous(labels=NULL) +
      ylab("Power (standardised spectrum)") +
      xlab("Cycle frequency") +
      ggtitle("") +
      theme_light(base_size=14)+
      theme(legend.position = "none") +

      geom_text(data = spectrum_output_species,
                aes(x = xloc,
                    y = yloc,
                    label = cycle_category)) +
      geom_text(data = spectrum_output_species,
                aes(x = xloc,
                    y = yloc2,
                    label = paste("cycle = ",cycle_months, "months")))


    # Periodograms
    if(perio_plot == TRUE && sig_95 == TRUE){
      plot_name <- paste("~/Desktop/FFT/",species_name[j], "_", pheno,".png",sep = "")
      png(plot_name, width = 925, height = 700)
      grid.arrange(timeline_plot, periodogram, widths = c(1,1))
      dev.off()
    }

  }

  return(spectrum_output)
}
