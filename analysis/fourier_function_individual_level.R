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

# # Main spectrum function from Bush
# spec_fun <- function(x) spectrum(x, #spans = c(3,3), # adding spans changes significance of dominant frequencies found!!!
#                                  plot=F, demean=T, detrend=T)
# # spectrum function for null hypothesis spectrum
# spec_null_fun <- function(x) spectrum(x, spans=c(55,57), # a large span from Bush c(55,57)
#                                       plot=F, demean=T, detrend=T)
# # Function for co-Fourier analysis of a simulated cosine curve against the empirical data
# spec_cross_fun <- function(x) spec.pgram(x,#spans=
#                                          plot=F, demean=T, detrend=T)


spans_lookup <- list(observations = c(96,144, 192, 240, 288, 336, 384, 432, 480, 528,576,624,672,720,768,816,864,912,960, 1008),
                     spans_smooth = list(NULL, NULL, NULL, NULL,
                                       NULL, NULL, NULL, NULL,
                                       3,3,3,3,
                                       3,3,3,
                                       c(3,3), c(3,3), c(3,3), c(3,3)),
                   # spans_smooth=list(NULL, NULL, NULL, NULL,
                   #                   3,3,3,3,
                   #                   c(3,3), c(3,3), c(3,3),
                   #                   c(3,5), c(3,5), c(3,5),
                   #                   c(5,5), c(5,5), c(5,5),
                   #                   c(5,7), c(5,7), c(5,7)),
                   spans_super_smooth = list(c(5,7), c(7,9), c(11,11), c(13,13),
                                           c(15,17), c(19,19), c(19,21), c(23,23),
                                           c(25,27), c(27,29), c(29,31),
                                           c(33,33), c(35,37), c(37,39),
                                           c(39,41), c(45,45), c(45,47),
                                           c(49,51), c(49,51), c(49,51)))

# spectrum function for normal smoother periodogram
spec_fun <- function(x) spectrum(x,spans = spans_lookup$spans_smooth[[which.min(abs(spans_lookup$observations - length(x)))]],
                               plot = F, demean = T, detrend = T)

# spectrum function for null hypothesis spectrum (super-smoothed spectrum of data -> smoothed periodogram bandwith to approx. 1)
spec_null_fun <- function(x) spectrum(x, spans = spans_lookup$spans_super_smooth[[which.min(abs(spans_lookup$observations - length(x)))]],
                                      plot = F, demean = T, detrend = T)

#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
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
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
# peak_detection <- function(
#   data,
#   species_name,
#   pheno,
#   # species_name = "Anonidium mannii",
#   # pheno = "flowers",
#   perio_plot = TRUE){
#
# # species_name = c("Irvingia grandifolia", "Erythrophleum suaveolens")
# # pheno = "flowers"
# # data_subset2 <- data %>%
# #   filter(species_full %in% species_name) %>%
# #   filter(phenophase %in% pheno)
#
#   spectrum_output <- data.frame()
#
#   for (j in 1:length(species_name)){


  #----------------------------------------------------------------------
  #----- subset data as timeseries
  #----------------------------------------------------------------------
  # # # # get out specific individual
  # # if(missing(species_name[j])){
  # #   data_subset <- data
  # # } else {
  #   data_subset <- data %>%
  #     filter(species_full %in% species_name[j]) %>%
  #     filter(phenophase %in% pheno)
  # # }
  species_name = "Albizia adianthifolia" #Erythrophleum suaveolens, Irvingia grandifolia
  pheno = "leaf_turnover"
  data_subset <- data %>%
    filter(species_full %in% species_name) %>%
    filter(phenophase %in% pheno)

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

  #---
  # grow dataset to full range
  #---
  individuals <- unique(data_subset$id)
  data_grow <- data.frame()
  for (i in 1:length(individuals)){
    data_ind_grow <- data_subset %>%
      filter(id %in% individuals[i])
    years <- sort(rep(min(data_ind_grow$year):max(data_ind_grow$year), 48))
    days <- rep(round((1:48 * 7.6) - 7),length(unique(years)))
    test_full <- data.frame(date = as.Date(paste(years, days, sep = "-"), "%Y-%j"))

    data_ind_grow <- merge(data_ind_grow, test_full, by = "date", all.y = TRUE)
    data_ind_grow$species_full <- unique(na.omit(data_ind_grow$species_full)) #fill species in empty years
    data_ind_grow$id <- unique(na.omit(data_ind_grow$id))
    data_ind_grow$year <- format.Date(data_ind_grow$date, "%Y")

    data_grow <- rbind(data_grow, data_ind_grow)
  }

  #---
  # fill in short missing consecutive years
  #---
  # group consecutive missing years
  data_id <- data_subset %>%
    filter(id %in% "353")
  a <- unique(data_id$year)

  data_grow_id <- data_grow %>%
    filter(id %in% "353")
  b <- as.numeric(unique(data_grow_id$year))
  missing_years <- setdiff(b,a)
  consec_years <- cumsum(c(1, abs(missing_years[-length(missing_years)] - missing_years[-1]) > 1))
  missing_years <- as.data.frame(cbind(missing_years,consec_years))
  colnames(missing_years)[1] <- 'year'
  test <- missing_years %>%
    group_by(consec_years) %>%
    dplyr::summarise(lgh_consec_years = length(consec_years))
  missing_years <- merge(missing_years, test, by = "consec_years", all.x = TRUE)
  data_ind_full <- merge(data_grow_id, missing_years, by = "year", all.x = TRUE)
  # if only 2 consecutive years of missing data -> fill with zeros. If longer, keep NA
  data_ind_full$value <- ifelse(is.na(data_ind_full$value) & data_ind_full$lgh_consec_years < 3 , 0, data_ind_full$value)

  #---
  # only keep longest consecutive timeline for fourier analysis
  #---
  ### repeat to only keep longest timeline for fourier analysis
  # remove rows with NA's in values -> individuals with 'no_data' in the archive
  timelines <- data_ind_full[!(is.na(data_ind_full$value)),]
  years_timeline <- as.numeric(unique(timelines$year))
  consec_timelines <- cumsum(c(1, abs(years_timeline[-length(years_timeline)] - years_timeline[-1]) > 1))
  years_timeline <- as.data.frame(cbind(years_timeline,consec_timelines))
  colnames(years_timeline)[1] <- 'year'
  test <- years_timeline %>%
    group_by(consec_timelines) %>%
    dplyr::summarise(lgh_consec_timelines = length(consec_timelines))
  years_timeline <- merge(years_timeline, test, by = "consec_timelines", all.x = TRUE)

  timelines <- merge(timelines, years_timeline, by = "year", all.x = TRUE)

  timelines$value <- ifelse(timelines$lgh_consec_timelines == max(timelines$lgh_consec_timelines) , timelines$value, NA)
  timelines <- timelines[!(is.na(timelines$value)),]

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

  # timeseries plot of raw data for each individual
  raw_plots <- ggplot(data_grow,
                    aes(x = date,
                        y = value))+
    geom_line()+
    scale_y_continuous(labels=NULL)+
    ylab("Observations") +
    xlab("") +
    ggtitle(paste(species_name[j], "-", pheno)) +
    theme_light(base_size = 14) +
    theme(legend.position = "none") +
    facet_wrap(~id,
               ncol = 1,
               strip.position = "left")

  # grid.arrange(raw_plots, Phenology_circular_boxplot, ncol=2, widths = c(2,1))

  #----------------------------------------------------------------------
  #----- fourier transform + periodogram
  #----------------------------------------------------------------------
  fourier_df <- data.frame()
  spectrum_output_species <- data.frame()

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
                               spec_norm = spec_fun(data_ts)$spec*(1/mean(spec_fun(data_ts)$spec)),
                               freq_null = spec_null_fun(data_ts)$freq/48,
                               spec_null = spec_null_fun(data_ts)$spec,
                               spec_norm_null = spec_null_fun(data_ts)$spec*(1/mean(spec_null_fun(data_ts)$spec)))

  fourier_df <- rbind(fourier_df,fourier_df_ind)


  # # extract key variables from spectrum outputs
  #----- Periodogram analysis with confidence intervals to find frequency and phase
  freq_dom = (spec_fun(data_ts)$freq[which.max(spec_fun(data_ts)$spec)])/48 # frequency of the dominant peak
  cycle_dom = 1/freq_dom
  spec_dom = max(spec_fun(data_ts)$spec) # spectrum of dominant peak
  spec_dom_norm =  spec_dom*(1/mean(spec_fun(data_ts)$spec)) # normalised spectrum so that mean power is equal to 1 enabling comparisons of "Power" between individuals
  dfree = spec_fun(data_ts)$df # degrees of freedom of the smoothed periodogram
  lower_ci = (dfree*spec_dom)/(qchisq(c(0.975),dfree)) # lower CI of dominant peak of smoothed periodogram
  spec_null_dom <- spec_null_fun(data_ts)$spec[which(abs((spec_null_fun(data_ts)$freq)/48-freq_dom)==min(abs((spec_null_fun(data_ts)$freq)/48-freq_dom)))] # the corresponding dominant peak from the null continuum
  sig = lower_ci > spec_null_dom # If the null spectrum is lower than the lower confidence interval, dominant peak can be considered "signficantly different" to null hypothesis
  cycle_category <- as.factor(ifelse(1/freq_dom > length(data_ts)/2,"No cyclicity",
                                     ifelse(1/freq_dom > 49,"Supra-annual",
                                            ifelse(1/freq_dom >= 47,"Annual","Sub-annual"))))
  # Is the dominant cycle both significant and represent a regular cycle?
  sig_cycle <- as.factor(ifelse(cycle_category!="No cyclicity" & sig==TRUE,TRUE,FALSE))

  spectrum_output_ind <- data.frame(species = species_name[j],
                                    phenophase = pheno,
                                    id = individuals[i],
                                    freq_dom = freq_dom,
                                    cycle_dom = cycle_dom,
                                    spec_dom = spec_dom,
                                    spec_dom_norm =  spec_dom_norm,
                                    dfree = dfree,
                                    lower_ci = lower_ci,
                                    spec_null_dom = spec_null_dom,
                                    sig = sig,
                                    cycle_category = cycle_category,
                                    sig_cycle = sig_cycle)
  spectrum_output_species <- rbind(spectrum_output_species, spectrum_output_ind)
  }


  ##Plot periodograms for visual inspection
  spectrum_output_species$xloc <- 0.4
  spectrum_output_species$yloc <- max(fourier_df$spec_norm, na.rm=T)*0.75
  spectrum_output_species$yloc2 <- max(fourier_df$spec_norm, na.rm=T)*0.5
  spectrum_output_species$cycle_months <- ifelse(spectrum_output_species$sig %in% "TRUE",
                                         format(round(spectrum_output_species$cycle_dom/4,0)), NA)

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
    facet_wrap(~id,
               ncol = 1,
               strip.position = "left") +
    geom_text(data = spectrum_output_species,
              aes(x = xloc,
                  y = yloc,
                  label = cycle_category)) +
    geom_text(data = spectrum_output_species,
              aes(x = xloc,
                  y = yloc2,
                  label = paste("cycle = ",cycle_months, "months")))

  # # Periodograms
  # if(perio_plot == TRUE){
  #   plot_name <- paste("~/Desktop/FFT/",species_name[j], "_", pheno,".png",sep = "")
  #   png(plot_name, width = 925, height = 700)
  #   grid.arrange(raw_plots, Periodograms, widths = c(1,1))
  #   dev.off()
  # } else {
  #   grid.arrange(raw_plots, Periodograms, widths = c(1,1))
  # }



  spectrum_output <- rbind(spectrum_output, spectrum_output_species)

#   }
#
#   return(spectrum_output)
# }



firsttest <- peak_detection(data,
                            species_name = c("Scorodophloeus zenkeri","Irvingia grandifolia"), #"Erythrophleum suaveolens",Irvingia grandifolia, Pericopsis elata
                            pheno = "flowers", #flowers, leaf_turnover, leaf_dormancy, fruit, fruit_drop
                            perio_plot = FALSE)


test <- firsttest %>%
  group_by(species) %>%
  dplyr::summarise(nr_id = length(id),
                   signif_perio = sum(sig=="TRUE"),
                   cycle = list(cycle_months[!is.na(cycle_months)]))

                     #count(which(sig==TRUE)))


# sectest <- data %>%
#   filter(phenophase == "flowers") %>%
#   filter(species_full %in% c("Antrocaryon nannanii",
#                              "Erythrophleum suaveolens",
#                              "Irvingia grandifolia",
#                              "Pericopsis elata",
#                              "Petersianthus macrocarpus")) %>% #,"Anonidium mannii"
#   group_by(species_full) %>%
#   do(peak_detection(.))

overview <- read.csv("data/SI_table2.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)

species_all <- overview$species_full

fulltest2 <- peak_detection(data,
                            species_name = species_all,
                            pheno = "leaf_turnover", #flowers, leaf_turnover, leaf_dormancy, fruit, fruit_drop
                            perio_plot = FALSE)

test2 <- fulltest2 %>%
  group_by(species) %>%
  dplyr::summarise(nr_id = length(id),
                   signif_perio = sum(sig=="TRUE"),
                   cycle = list(cycle_months[!is.na(cycle_months)]))

