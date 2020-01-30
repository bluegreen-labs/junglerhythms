#' Peak detection and cyclicity using fourier transform
#'
#' @param data junglerhythms data file
#' @param species_name end coordinate fo a line
#' @export
#' @return timing and stats for peaks

library(TSA)
library(tidyverse)
library(plyr)
library(ggplot2)
library(gtools)
library(grid)
library(gridExtra)
library(circular)
library(fields)
library(DescTools)
library(truncnorm)

# Main spectrum function from Bush
spec_fun <- function(x) spectrum(x, #spans = c(3,3), # adding spans changes significance of dominant frequencies found!!!
                                 plot=F, demean=T, detrend=T)
# spectrum function for null hypothesis spectrum
spec_null_fun <- function(x) spectrum(x, spans=c(55,57), # a large span from Bush c(55,57)
                                      plot=F, demean=T, detrend=T)
# Function for co-Fourier analysis of a simulated cosine curve against the empirical data
spec_cross_fun <- function(x) spec.pgram(x,#spans=
                                         plot=F, demean=T, detrend=T)


spans_lookup <- list(obs=c(96,144, 192, 240, 288, 336, 384, 432, 480, 528,576,624,672,720,768,816,864,912,960, 1008),
                     # spans_smooth=list(NULL, NULL, NULL, NULL,
                     #                   NULL, NULL, NULL, NULL,
                     #                   3,3,3,3,
                     #                   3,3,3,
                     #                   c(3,3), c(3,3), c(3,3), c(3,3)),
                   spans_smooth=list(NULL, NULL, NULL, NULL,
                                     3,3,3,3,
                                     c(3,3), c(3,3), c(3,3),
                                     c(3,5), c(3,5), c(3,5),
                                     c(5,5), c(5,5), c(5,5),
                                     c(5,7), c(5,7), c(5,7)),
                   spans_super_smooth=list(c(5,7), c(7,9), c(11,11), c(13,13),
                                           c(15,17), c(19,19), c(19,21), c(23,23),
                                           c(25,27), c(27,29), c(29,31),
                                           c(33,33), c(35,37), c(37,39),
                                           c(39,41), c(45,45), c(45,47),
                                           c(49,51), c(49,51), c(49,51)))

#Function to calculate smoothed spectrum - > same function as before
spec_fun <- function(x) spectrum(x,spans=spans_lookup$spans_smooth[[which.min(abs(spans_lookup$obs-length(x)))]],
                               plot=F,demean=T,detrend=T) #spectrum function for normal smoother periodogram

#Function to calculate the null continuum (super-smoothed spectrum of data) for null hypothesis test
spec_null_fun<-function(x) spectrum(x,spans=spans_lookup$spans_super_smooth[[which.min(abs(spans_lookup$obs-length(x)))]],
                                    plot=F,demean=T,detrend=T) #spectrum function for null hypothesis spectrum


# #----
# # span detection
# #-----
# species_name = "Scorodophloeus zenkeri"
# pheno = "flowers"
# data_subset <- data %>%
#   filter(species_full %in% species_name) %>%
#   filter(phenophase %in% pheno)
#
# data_subset$date <- as.Date(paste(data_subset$year,
#                                   round((data_subset$week*7.6)-7),sep="-"), "%Y-%j")
# data_subset$month <- format.Date(data_subset$date, "%m")
# # average by date
# data_subset <- data_subset %>%
#   group_by(species_full, date) %>%
#   dplyr::summarise(mean_value = mean(value),
#                    scaled_value = ifelse(any(value > 0), 1, 0))
# # sort dataframe according to date
# data_subset <- data_subset %>%
#   dplyr::arrange(date)
#
# data_subset <- data_subset %>%
#   filter(date > '1950-12-31')
#
# first_year <- as.numeric(format.Date(data_subset$date[1], "%Y"))
# last_year <- as.numeric(format.Date(last(data_subset$date), "%Y"))
# data_ts <- ts(data = data_subset$mean_value, start = c(first_year,1), frequency = 48) #, end = c(last_year,48)
#
#
# spec_fun <- function(x) spectrum(x, spans = NULL,#c(3,3),
#                                  plot=F,demean=T,detrend=T)
# spec_fun(data_ts)$bandwidth

#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
peak_detection <- function(
  data,
  species_name,
  pheno,
  # species_name = "Anonidium mannii",
  # pheno = "flowers",
  perio_plot = TRUE
){

  #----------------------------------------------------------------------
  #----- subset data as timeseries
  #----------------------------------------------------------------------
  # # get out specific individual
  if(missing(species_name)){
    data_subset <- data
  } else {
    data_subset <- data %>%
      filter(species_full %in% species_name) %>%
      filter(phenophase %in% pheno)
  }
  # species_name = "Irvingia grandifolia" #Erythrophleum suaveolens
  # pheno = "flowers"
  # data_subset <- data %>%
  #   filter(species_full %in% species_name) %>%
  #   filter(phenophase %in% pheno)

  # convert date
  data_subset$date <- as.Date(paste(data_subset$year,
                                    round((data_subset$week*7.6)-7),sep="-"), "%Y-%j")
  # add a month column
  data_subset$month <- format.Date(data_subset$date, "%m")

  # # average by date
  # data_subset <- data_subset %>%
  #   group_by(species_full, date) %>%
  #   dplyr::summarise(mean_value = mean(value),
  #                    scaled_value = ifelse(any(value > 0), 1, 0))

  # sort dataframe according to date
  data_subset <- data_subset %>%
    dplyr::arrange(id,date)

  #-----------------------------
  # overview plots (from Bush)
  #-----------------------------

  # Summarise phenology scores across sample for each month and year
  # (to give proportions for boxplots)
  data_df_sample_summary <- ddply(data_subset,.(year,month), summarize,
                                propPhenophase = length(value[value>0])/length(id))

  data_df_sample_summary$month <- mapvalues(data_df_sample_summary$month,
                                            from=c("01","02","03","04","05","06","07","08","09","10","11","12"),
                                            to=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
  data_df_sample_summary$month <- factor(data_df_sample_summary$month,
                                         ordered=T,
                                         levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

  # Write circular boxplot graphic
  Phenology_circular_boxplot <- ggplot(data_df_sample_summary,
                                       aes(x=factor(month),
                                           y=propPhenophase)) +
    geom_boxplot(fill="orange",outlier.colour="lightgrey") +
    coord_polar(start=0) +
    theme_minimal(base_size=16) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.ticks = element_blank()) +
    ylab("")

  # Figure 1b - timeseries plot of raw data for each individual
  raw_plots <- ggplot(data_subset,
                    aes(x = date,
                        y = value))+
    geom_line()+
    scale_y_continuous(labels=NULL)+
    ylab("Observations") +
    xlab("") +
    ggtitle(paste(species_name, "-", pheno)) +
    theme_light(base_size = 14) +
    theme(legend.position = "none") +
    #       # strip.background = element_blank(),
    #       strip.text = element_blank()) +
    facet_wrap(~id,
               ncol = 1,
               strip.position = "left")

  # grid.arrange(raw_plots, Phenology_circular_boxplot, ncol=2, widths = c(2,1))

  #----------------------------------------------------------------------
  #----- fourier transform + periodogram
  #----------------------------------------------------------------------
  fourier_df <- data.frame()
  spectrum_output <- data.frame()

  individuals <- unique(data_subset$id)

  # time series at the individual level

  for (i in 1:length(individuals)){
  data_ind <- data_subset %>%
    filter(id %in% individuals[i]) #
  first_year <- as.numeric(format.Date(data_ind$date[1], "%Y"))
  last_year <- as.numeric(format.Date(last(data_ind$date), "%Y"))

  data_ts <- ts(data = data_ind$value, start = c(first_year,1), frequency = 48) #, end = c(last_year,48)

  # fourier transform
  fourier_df_ind <- data.frame(id = individuals[i],
                               freq = spec_fun(data_ts)$freq/48,
                               spec = spec_fun(data_ts)$spec,
                               spec_norm = spec_fun(data_ts)$spec*(1/mean(spec_fun(data_ts)$spec)),
                               freq_null = spec_null_fun(data_ts)$freq/48,
                               spec_null = spec_null_fun(data_ts)$spec,
                               spec_norm_null = spec_null_fun(data_ts)$spec*(1/mean(spec_null_fun(data_ts)$spec)))

  fourier_df <- rbind(fourier_df,fourier_df_ind)


  # # extract key variables from spectrum outputs
  #----- Periodogram analysis with confidence intervals to find frequency and phase
  freq_dom = (spec_fun(data_ts)$freq[which.max(spec_fun(data_ts)$spec)])/48 #frequency of the dominant peak
  cycle_dom = 1/freq_dom
  spec_dom = max(spec_fun(data_ts)$spec)
  spec_dom_norm =  spec_dom*(1/mean(spec_fun(data_ts)$spec))
  dfree = spec_fun(data_ts)$df
  lower_ci = (dfree*spec_dom)/(qchisq(c(0.975),dfree))
  spec_null_dom <- spec_null_fun(data_ts)$spec[which(abs((spec_null_fun(data_ts)$freq)/48-freq_dom)==min(abs((spec_null_fun(data_ts)$freq)/48-freq_dom)))] #the corresponding dominant peak from the null continuum
  sig = lower_ci > spec_null_dom
  cycle_category <- as.factor(ifelse(1/freq_dom > length(data_ts)/2,"No cyclicity",
                                     ifelse(1/freq_dom > 49,"Supra-annual",
                                            ifelse(1/freq_dom >= 47,"Annual","Sub-annual"))))
  sig_cycle <- as.factor(ifelse(cycle_category!="No cyclicity" & sig==TRUE,TRUE,FALSE))

  spectrum_output_ind <- data.frame(species = species_name,
                                    phenophase = pheno,
                                    id = individuals[i],
                                    freq_dom = freq_dom, #frequency of the dominant peak
                                    cycle_dom = cycle_dom,
                                    spec_dom = spec_dom, #spectrum of dominant peak
                                    spec_dom_norm =  spec_dom_norm, #normalised spectrum so that mean power is equal to 1 enabling comparisons of "Power" between individuals
                                    dfree = dfree, #degrees of freedom of the smoothed periodogram
                                    lower_ci = lower_ci, #lower CI of dominant peak of smoothed periodogram
                                    spec_null_dom = spec_null_dom, #the corresponding dominant peak from the null continuum
                                    sig = sig, #If the null spectrum is lower than the lower confidence interval, dominant peak can be considered "signficantly different" to null hypothesis
                                    cycle_category = cycle_category,
                                    # Is the dominant cycle both significant and represent a regular cycle?
                                    sig_cycle = sig_cycle)
  spectrum_output <- rbind(spectrum_output, spectrum_output_ind)
  }

  ##Plot periodograms for visual inspection
  spectrum_output$xloc <- 0.4
  spectrum_output$yloc <- max(fourier_df$spec_norm, na.rm=T)*0.75
  spectrum_output$yloc2 <- max(fourier_df$spec_norm, na.rm=T)*0.5
  spectrum_output$cycle_months <- ifelse(spectrum_output$sig %in% "TRUE",
                                         format(round(spectrum_output$cycle_dom/4,0)), "ns")

  Periodograms <- ggplot(fourier_df) +
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
    # annotate("text", label = c("id","2"), size = 4, x = 0.4, y = 60) +
    facet_wrap(~id,
               ncol = 1,
               strip.position = "left") +
    geom_text(data = spectrum_output,
              aes(x = xloc,
                  y = yloc,
                  label = cycle_category)) +
    geom_text(data = spectrum_output,
              aes(x = xloc,
                  y = yloc2,
                  label = paste("cycle = ",cycle_months, "months")))

  # Periodograms
  grid.arrange(raw_plots, Periodograms, widths = c(1,1))

  # #Summarise dominant peaks
  # fourier_df_sample_summary <- ddply(fourier_df, .(id), summarize,
  #                                  maxFreq = freq[which.max(spec)],
  #                                  maxSpec_norm = max(spec_norm),
  #                                  cycle = 1/maxFreq)
  #
  # #Dominant phenology cycle for sample
  # Median_cycle <- 1/median(spectrum_output$freq_dom)
  # Median_cycle



  # ###################################################################################
  # # out-of-phase detection comparing to basic cosine with same frequency
  # # phase difference not working will
  # # seems like based on first years?
  # ###################################################################################
  # #-----
  # # Find phase of significant dominant cycles
  # #-----
  # #Function to extract data from co-Fourier anaysis regarding phase difference in radians
  # #Simulate a "template" time series using the sample modal cycle peaking at the beginning of January 1937 to act as guide in co-Fourier analysis
  # w = 2*pi*(1/cycle_dom) #wavelength in radians
  # simulated_phase_ts <- ts(4*cos(w*1:720),start=c(first_year,1), end=c(last_year,48),freq=48) #1:360 -> 1:720 -> no interuption in cosinus for longer time series
  # ts_cross <- ts.intersect(simulated_phase_ts,data_ts)
  # phase_radians <- spec_cross_fun(ts_cross)$phase[which.max(spec_cross_fun(ts_cross)$spec[,1])] #the phase difference (in radians) at the dominant frequency
  # # different formats of timing
  # phase_degree <- deg(phase_radians)
  # phase_week <- phase_degree*48/360
  # phase_doy <- round((phase_week*7.6))
  # timing_peak <- as.Date(phase_doy, origin = "1938-01-01")
  # timing_peak <- format.Date(timing_peak, "%m-%d")
  #
  # # return data as data frame
  # return(data.frame(freq_dom = freq_dom,
  #                   cycle_dom = cycle_dom,
  #                   sig = sig,
  #                   cycle_category = cycle_category,
  #                   phase_radians = phase_radians,
  #                   phase_week = phase_week,
  #                   timing_peak = timing_peak
  # ))

  return(spectrum_output)
}



firsttest <- peak_detection(data,
                            species_name = c("Pericopsis elata"), #"Erythrophleum suaveolens",Irvingia grandifolia
                            pheno = "fruit") #flowers, leaf_turnover, leaf_dormancy, fruit, fruit_drop

# sectest <- data %>%
#   filter(phenophase == "flowers") %>%
#   filter(species_full %in% c("Antrocaryon nannanii",
#                              "Erythrophleum suaveolens",
#                              "Irvingia grandifolia",
#                              "Pericopsis elata",
#                              "Petersianthus macrocarpus")) %>% #,"Anonidium mannii"
#   group_by(species_full) %>%
#   do(peak_detection(.))
