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
spec_fun <- function(x) spectrum(x, #spans = c(3,3),
                                 plot=F,demean=T,detrend=T)
# spectrum function for null hypothesis spectrum
spec_null_fun <- function(x) spectrum(x,spans=c(55,57), # a large span from Bush
                                      plot=F,demean=T,detrend=T)
# Function for co-Fourier analysis of a simulated cosine curve against the empirical data
spec_cross_fun <- function(x) spec.pgram(x,#spans=
                                         plot=F, demean=T,detrend=T)

# peak_detection <- function(
#   data,
#   species_name,
#   pheno,
#   # species_name = "Anonidium mannii",
#   # pheno = "flowers",
#   perio_plot = TRUE
# ){

  #----------------------------------------------------------------------
  #----- subset data as timeseries
  #----------------------------------------------------------------------
  # # get out specific individual
  # if(missing(species_name)){
  #   data_subset <- data
  # } else {
  #   data_subset <- data %>%
  #     filter(species_full %in% species_name) %>%
  #     filter(phenophase %in% pheno)
  # }
data_subset <- data %>%
  filter(species_full %in% "Irvingia grandifolia") %>%
  filter(phenophase %in% "fruit")

  # convert date
  data_subset$date <- as.Date(paste(data_subset$year,
                                    round((data_subset$week*7.6)-7),sep="-"), "%Y-%j")

  # # average by date
  # data_subset <- data_subset %>%
  #   group_by(species_full, date) %>%
  #   dplyr::summarise(mean_value = mean(value),
  #                    scaled_value = ifelse(any(value > 0), 1, 0))

  # sort dataframe according to date
  data_subset <- data_subset %>%
    dplyr::arrange(id,date)
  # add a month column
  data_subset$month <- format.Date(data_subset$date, "%m")

  # # data as timeseries, species_level
  # first_year <- as.numeric(format.Date(data_subset$date[1], "%Y"))
  # last_year <- as.numeric(format.Date(last(data_subset$date), "%Y"))
  # data_ts <- ts(data = data_subset$scaled_value, start = first_year, frequency = 48)
  # str(data_ts)

  #-------------------------
  # overview plots
  #-----------------------------

  #Summarise phenology scores across sample for each month and year
  #(to give proportions for boxplots)
  data_df_sample_summary <- ddply(data_subset,.(year,month),summarize,
                                  propPhenophase = length(value[value > 0])/length(id))
  data_df_sample_summary$month <- mapvalues(data_df_sample_summary$month,from=c("01","02","03","04","05","06","07","08","09","10","11","12"), to=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
  data_df_sample_summary$month <- factor(data_df_sample_summary$month,ordered=T,levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

  #Write circular boxplot graphic
  Phenology_circular_boxplot <- ggplot(data_df_sample_summary,aes(x=factor(month),y=propPhenophase))+
    geom_boxplot(fill="orange",outlier.colour="lightgrey")+
    coord_polar(start=0)+
    theme_minimal(base_size=16)+
    theme(legend.position="none",axis.title.x=element_blank(),axis.ticks=element_blank())+
    ylab("")
  # raw data plot
  p_lin <- ggplot(data_subset) +
    geom_line(aes(x = date,
                  y = value)) +
    scale_x_date(date_breaks = "1 years",
                 date_labels = "%Y",
                 limits = as.Date(c('1937-01-01','1956-12-31'))) +
    theme_minimal() +
    theme(panel.grid.major.x = element_line(colour = "grey50", size = 0.3))
  #Periogorams for each individual
  periodogram <- ggplot(Fourier_df) +
    geom_line(aes(x = freq,
                  y = spec_norm)) +
    theme_minimal(base_size=14)+
    ylab("Power (standardised spectrum)")+
    xlab("Cycle frequency") +
    labs(title = unique(data_subset$species_full))
  # empty plot
  empty <- ggplot() + theme_void()
  #all
  # grid.arrange(periodogram,p_lin, heights = c(1,1))
  grid.arrange(periodogram, Phenology_circular_boxplot, p_lin, empty, ncol=2, widths = c(3,1))


  #-----------
  data_df_sample_summary<-ddply(data_subset,.(year,month),summarize,
                                propPhenophase=length(value[value>0])/length(id))

  data_df_sample_summary$month<-mapvalues(data_df_sample_summary$month,from=c("01","02","03","04","05","06","07","08","09","10","11","12"), to=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
  data_df_sample_summary$month<-factor(data_df_sample_summary$month,ordered=T,levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

  #Write circular boxplot graphic
  Phenology_circular_boxplot<-ggplot(data_df_sample_summary,aes(x=factor(month),y=propPhenophase))+
    geom_boxplot(fill="orange",outlier.colour="lightgrey")+
    coord_polar(start=0)+
    theme_minimal(base_size=16)+
    theme(legend.position="none",axis.title.x=element_blank(),axis.ticks=element_blank())+
    ylab("")

  ##Figure 1b - timeseries plot of raw data for each individual

  raw_plots<-ggplot(data_subset,aes(x=date,y=value,colour=id))+
    geom_line()+
    scale_y_continuous(labels=NULL)+
    ylab("Canopy score")+
    xlab("")+
    theme_light(base_size=14)+
    theme(#legend.position="none",
          strip.background = element_blank()
          , strip.text = element_blank()) +
    facet_wrap(~id,ncol=1)

  Phenology_circular_boxplot
  raw_plots

  #----------------------------------------------------------------------
  #----- fourier transform + periodogram
  #----------------------------------------------------------------------
  fourier_df <- data.frame()
  spectrum_output <- data.frame()
  individuals <- unique(data_subset$id)

  # time series at the individual level

  for (i in 1:length(individuals)){
  data_ind <- data_subset %>%
    filter(id %in% individuals[i])
  first_year <- as.numeric(format.Date(data_ind$date[1], "%Y"))
  last_year <- as.numeric(format.Date(last(data_ind$date), "%Y"))

  data_ts <- ts(data = data_ind$value, start = c(first_year,1), frequency = 48) #, end = c(last_year,48)

  # fourier transform
  fourier_df_ind <- data.frame(id = individuals[i],
                               freq = spec_fun(data_ts)$freq/48,
                               spec = spec_fun(data_ts)$spec,
                               spec_norm = spec_fun(data_ts)$spec*(1/mean(spec_fun(data_ts)$spec)))

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

  spectrum_output_ind <- data.frame(id = individuals[i],
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

  # #Periogorams for each individual on one plot - Figure 1c
  # Periodograms <- ggplot(data = fourier_df,aes(x = freq,
  #                                              y = spec_norm,
  #                                              group = id,
  #                                              colour = id)) +
  #   geom_line()+
  #   theme_minimal(base_size=14)+
  #   ylab("Power (standardised spectrum)")+
  #   xlab("Cycle frequency (cycles per month)")+
  #   ggtitle("Smoothed periodograms of phenology cycles by individual")

  Periodograms <- ggplot(fourier_df, aes(x = freq,
                                          y = spec_norm,
                                          colour = id))+
    geom_line()+
    scale_y_continuous(labels=NULL)+
    ylab("Power (standardised spectrum)")+
    xlab("Cycle frequency (cycles per month)")+
    ggtitle("Smoothed periodograms of phenology cycles by individual") +
    theme_light(base_size=14)+
    theme(#legend.position="none",
      strip.background = element_blank()
      , strip.text = element_blank()) +
    facet_wrap(~id, ncol=1)

  # Periodograms
  grid.arrange(raw_plots, Periodograms, widths = c(1,1))

  #Summarise dominant peaks
  fourier_df_sample_summary <- ddply(fourier_df, .(id), summarize,
                                   maxFreq = freq[which.max(spec)],
                                   maxSpec_norm = max(spec_norm),
                                   cycle = 1/maxFreq)

  #Dominant phenology cycle for sample
  Median_cycle <- 1/median(fourier_df_sample_summary$maxFreq)
  Median_cycle



  ###################################################################################
  # out-of-phase detection comparing to basic cosine with same frequency
  # phase difference not working will
  # seems like based on first years?
  ###################################################################################
  #-----
  # Find phase of significant dominant cycles
  #-----
  #Function to extract data from co-Fourier anaysis regarding phase difference in radians
  #Simulate a "template" time series using the sample modal cycle peaking at the beginning of January 1937 to act as guide in co-Fourier analysis
  w = 2*pi*(1/cycle_dom) #wavelength in radians
  simulated_phase_ts <- ts(4*cos(w*1:720),start=c(first_year,1), end=c(last_year,48),freq=48) #1:360 -> 1:720 -> no interuption in cosinus for longer time series
  ts_cross <- ts.intersect(simulated_phase_ts,data_ts)
  phase_radians <- spec_cross_fun(ts_cross)$phase[which.max(spec_cross_fun(ts_cross)$spec[,1])] #the phase difference (in radians) at the dominant frequency
  # different formats of timing
  phase_degree <- deg(phase_radians)
  phase_week <- phase_degree*48/360
  phase_doy <- round((phase_week*7.6))
  timing_peak <- as.Date(phase_doy, origin = "1938-01-01")
  timing_peak <- format.Date(timing_peak, "%m-%d")

  # return data as data frame
  return(data.frame(freq_dom = freq_dom,
                    cycle_dom = cycle_dom,
                    sig = sig,
                    cycle_category = cycle_category,
                    phase_radians = phase_radians,
                    phase_week = phase_week,
                    timing_peak = timing_peak
  ))
# }

#
#
# firsttest <- peak_detection(data,
#                             species_name = c("Erythrophleum suaveolens"), #"Anonidium mannii",
#                             pheno = "flowers")

sectest <- data %>%
  filter(phenophase == "leaf_turnover") %>%
  filter(species_full %in% c("Antrocaryon nannanii",
                             "Erythrophleum suaveolens",
                             "Irvingia grandifolia",
                             "Pericopsis elata",
                             "Petersianthus macrocarpus")) %>% #,"Anonidium mannii"
  group_by(species_full) %>%
  do(peak_detection(.))
