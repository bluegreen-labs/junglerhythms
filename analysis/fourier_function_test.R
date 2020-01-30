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
spec_fun <- function(x) spectrum(x,#spans = 3,
                                 plot=F,demean=T,detrend=T)
# spectrum function for null hypothesis spectrum
spec_null_fun <- function(x) spectrum(x,spans=c(55,57), # a large span from Bush
                                      plot=F,demean=T,detrend=T)
# Function for co-Fourier analysis of a simulated cosine curve against the empirical data
spec_cross_fun <- function(x) spec.pgram(x,#spans=
                                         plot=F, demean=T,detrend=T)

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
  # get out specific individual
  if(missing(species_name)){
    data_subset <- data
  } else {
    data_subset <- data %>%
      filter(species_full %in% species_name) %>%
      filter(phenophase %in% pheno)
  }

  # ids with to few years or events
  a <- data.frame(total_years = (tapply(data_subset$value, data_subset$id, length))/48,
                  total_events = tapply(data_subset$value, data_subset$id, sum))
  a$ID <- rownames(a)
  a <- a[!(a$total_years < 2 | a$total_events < 4),]
  ids_ok <- a$ID

  data_subset <- data_subset[(data_subset$id %in% ids_ok),]

  # convert date
  data_subset$date <- as.Date(paste(data_subset$year,
                                  round((data_subset$week*7.6)-7),sep="-"), "%Y-%j")
  # sort dataframe according to date
  data_subset <- data_subset %>%
    dplyr::arrange(date)



  # # data as timeseries, per individual
  # d_first_year <- data_subset %>%
  #   group_by(species_full, id) %>%
  #   dplyr::summarise(first_year = min(year))

  data_ls <- tapply(data_subset$value, data_subset$id, FUN = function(x)ts(x, start = 1938, frequency = 48))
  str(data_ls)
  IDs <- names(data_ls)

  #----------------------------------------------------------------------
  #----- fourier transform + periodogram
  #----------------------------------------------------------------------
  # fourier transform
  Fourier_df <- data.frame()
  for (i in 1:length(data_ls)){
    d <- data.frame(freq = spec_fun(data_ls[[i]])$freq,
                  spec = spec_fun(data_ls[[i]])$spec,
                  spec_norm = spec_fun(data_ls[[i]])$spec*(1/mean(spec_fun(data_ls[[i]])$spec)),
                  ID = as.factor(IDs[i]))
    Fourier_df<-rbind(Fourier_df,d)
  }


  # raw data plot
  p_lin <- ggplot(data_subset) +
    geom_line(aes(x = date,
                    y = value,
                    colour = id)) +
    scale_x_date(date_breaks = "1 years",
                   date_labels = "%Y",
                   limits = as.Date(c('1937-01-01','1956-12-31'))) +
    theme_minimal() +
    theme(panel.grid.major.x = element_line(colour = "grey50", size = 0.3)) +
    facet_wrap(~id,
              ncol=1)
  #Periogorams for each individual
  periodogram <- ggplot(data=Fourier_df) +
    geom_line(aes(x = freq,
                 y = spec_norm,
                 group = ID,
                 colour=ID)) +
    theme_minimal(base_size=14)+
    ylab("Power (standardised spectrum)")+
    xlab("Cycle frequency") +
    labs(title = unique(data_subset$species_full))

  grid.arrange(periodogram,p_lin, heights = c(1,1))

  #----------------------------------------------------------------------
  #----- Periodogram analysis with confidence intervals to find frequency and phase
  #----------------------------------------------------------------------
  # Function to extract key variables from spectrum outputs
  confidence_fun <- function(p){
    ts <- data_ls[[p]]
    d <- data.frame(ID=IDs[p])
    d$length <- length(ts)
    d$freq_dom <- (spec_fun(ts)$freq[which.max(spec_fun(ts)$spec)])/48 #frequency of the dominant peak
    d$cycle_dom <- 1/d$freq_dom
    d$spec_dom <- max(spec_fun(ts)$spec) #spectrum of dominant peak
    d$spec_dom_norm =  d$spec_dom*(1/mean(spec_fun(ts)$spec)) #normalised spectrum so that mean power is equal to 1 enabling comparisons of "Power" between individuals
    dfree <- spec_fun(ts)$df #degrees of freedom of the smoothed periodogram
    lower_ci <- (dfree*d$spec_dom)/(qchisq(c(0.975),dfree)) #lower CI of dominant peak of smoothed periodogram
    spec_null_dom <- spec_null_fun(ts)$spec[which(abs((spec_null_fun(ts)$freq)/48-d$freq_dom)==min(abs((spec_null_fun(ts)$freq)/48-d$freq_dom)))] #the corresponding dominant peak from the null continuum
    d$sig <- lower_ci>spec_null_dom #If the null spectrum is lower than the lower confidence interval, dominant peak can be considered "signficantly different" to null hypothesis
    data.frame(d)
  }

  # Run loop to apply spectrum functions to each individual time series in data_ls
  confidence_results <- ldply(1:length(data_ls),.fun = confidence_fun)

  # Categorise dominant cycles and screen results for positive false results where dominant peak is same as full length of the data - "No cyclicity"
  confidence_results$cycle_category<-as.factor(ifelse(1/confidence_results$freq_dom > confidence_results$length/2,"No cyclicity",
                                                    ifelse(1/confidence_results$freq_dom > 49,"Supra-annual",
                                                           ifelse(1/confidence_results$freq_dom >= 47,"Annual","Sub-annual"))))

  # Is the dominant cycle both significant and represent a regular cycle?
  confidence_results$sig_cycle <- as.factor(ifelse(confidence_results$cycle_category!="No cyclicity" & confidence_results$sig==TRUE,TRUE,FALSE))

  print(head(confidence_results))


  #-----
  # Find phase of significant dominant cycles and assess synchrony between individuals
  #-----

  #Function to extract data from co-Fourier anaysis regarding phase difference in radians
  #Simulate a "template" time series using the sample modal cycle peaking at the beginning of January 1937 to act as guide in co-Fourier analysis
  phase_fun <- function(p){
    w = 2*pi*(1/confidence_results$cycle_dom[p]) #wavelength in radians
    # print(w)
    simulated_phase_ts <- ts(4*cos(w*1:360),start=c(1938),freq=48)# during make ts, all started in 1937 # d_first_year$first_year[p]
    # print(simulated_phase_ts)
    ts_cross <- ts.intersect(simulated_phase_ts,data_ls[[p]])
    d<-data.frame(ID=IDs[p])
    # d$cycle_used <- w[p]
    d$phase_radians<-spec_cross_fun(ts_cross)$phase[which.max(spec_cross_fun(ts_cross)$spec[,1])] #the phase difference (in radians) at the dominant frequency
    d
  }


  phase_results <- ldply(1:length(data_ls),.fun = phase_fun)
  # different formats of timing
  phase_results$phase_degree <- deg(phase_results$phase_radians)
  phase_results$phase_week <- phase_results$phase_degree*48/360
  phase_results$phase_doy <- round((phase_results$phase_week*7.6))
  phase_results$timing_peak <- as.Date(phase_results$phase_doy, origin = "1937-01-01")
  phase_results$timing_peak <- format.Date(phase_results$timing_peak, "%m-%d")

  print(head(phase_results))

  ##Combine confidence and phase results
  results <- merge(a,confidence_results,"ID")
  results <- merge(results,phase_results,"ID")

  print(head(results))


  # return data as data frame
  return(results)
  # return(data.frame(NULL))
}

#
#
# firsttest <- peak_detection(data,
#                             species_name = c("Erythrophleum suaveolens"), #"Anonidium mannii",
#                             pheno = "flowers")

sectest <- data %>%
  filter(phenophase == "flowers") %>%
  filter(species_full %in% c("Irvingia grandifolia")) %>% #,"Anonidium mannii", "Erythrophleum suaveolens",
  group_by(species_full) %>%
  do(peak_detection(.))
