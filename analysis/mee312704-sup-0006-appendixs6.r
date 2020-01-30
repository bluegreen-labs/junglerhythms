#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
#-------------------------------------------------------------------------------#

############################
##Supporting information: R code for Fourier analysis of phenology
# Bush et al. "Fourier analysis to detect phenological cycles using tropical field data and simulations."
# 25th October 2016
############################

# This code is
# designed to work for a single species sample and will need to be scaled-up for
# analysis of more than one species

# line 19 - Upload libraries
# line 35 - Prep data for analysis
# line 115 - Initial analysis of phenological pattern (circular boxplots from Fig 1a)
# line 153 - Intitial observation of the Periodogram (peridograms per individual from Fig 1c)
# line 208 - Periodogram analysis with confidence intervals (extracting data for Figures 1c-e, 2 and 3)
# line 332 - Power analyses of simulated data (Figure 4)

##########################
#### Upload libraries ####
##########################

library(plyr)
library(ggplot2)
library(reshape)
# library("RColorBrewer", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")
library(RColorBrewer)
library(gtools)
# library("grid", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")
library(grid)
library(circular)
# library("fields", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")
library(fields)
library(DescTools)
# library("grid", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")
library(truncnorm)

##############################
####Prep data for analysis####
##############################

### Data format introduction ###

# To work this code, empirical or simulated data should be arranged initially as a named list of
# individual time series (called "data_ls"). These should each be of class “time series” and
# have no gaps (any small gaps, for example up to 3 data points long, could be
# interpolated using a simple linear estimator). This will then also be arranged
# into data frame format (data_df) for some graphics / analyses.
#  e.g. str(data_ls)
#   List of 10
#  $ 1 : Time-Series [1:360] from 1986 to 2016: 3.4934 2.0508 0.0587 0 0 ...
#  $ 2 : Time-Series [1:360] from 1986 to 2016: 3.525 2.105 0.121 0 0 ...
#  $ 3 : Time-Series [1:360] from 1986 to 2016: 3.846 2.661 0.763 0 0 ...
#  $ 4 : Time-Series [1:360] from 1986 to 2016: 3.725 2.452 0.522 0 0 ...
#  $ 5 : Time-Series [1:360] from 1986 to 2016: 3.17 1.49 0 0 0 ...
#  $ 6 : Time-Series [1:360] from 1986 to 2016: 3.19 1.52 0 0 0 ...
#  $ 7 : Time-Series [1:360] from 1986 to 2016: 3.3 1.71 0 0 0 ...
#  $ 8 : Time-Series [1:360] from 1986 to 2016: 3.917 2.785 0.906 0 0 ...
#  $ 9 : Time-Series [1:360] from 1986 to 2016: 3.828 2.63 0.728 0 0 ...
#  $ 10: Time-Series [1:360] from 1986 to 2016: 3.03 1.25 0 0 0 ...

### Simulate data ###

#The following code for simulating data can be used as an example set to run
# all following code

# Set fundamental variables for the simulated time series sample. The paramaters
#below are identical for each  individual, but these can also be varied by individual if required
individuals<-10 #sample size
months<-rep(360, length.out=individuals) #length of timeseries (months)
sd<-rep(2,length.out=individuals) #level of SD around random noise term
wavelength<-rep(12,length.out=individuals) #Dominant wavelength (months)
p<-rep(0,length.out=individuals) #phase
R<-rep(4,length.out=individuals) #amplitude of sinewave

#Derive further parameters from those set above
f<-1/wavelength #frequency
w<-2*pi*f #wavelength in radians
z<-rep(NA,length.out=individuals) #empty object for random noise term
t<-rep(NA,length.out=individuals) #empty object for time signifier
ts<-rep(NA,length.out=individuals) #empty object for sine curve output
names<-paste(months,sd,wavelength,sep=",")

#Assemble the variables above and run loop to populate time series with data
conds<-data.frame(ID=c(1:individuals),name=names,t,months,sd,wavelength,f,w,R,p,z,ts)
conds<-lapply(split(conds,conds$ID),function(x) as.list(x))

for (i in 1:individuals){
  conds[[i]]$t<-1:conds[[i]]$months
  conds[[i]]$z<-rnorm(conds[[i]]$t,mean=0,sd=conds[[i]]$sd)
  conds[[i]]$ts<-conds[[i]]$R*cos(conds[[i]]$w*conds[[i]]$t+conds[[i]]$p)+conds[[i]]$z
  conds[[i]]$ts[conds[[i]]$ts<0]<-0
  conds[[i]]$ts<-(conds[[i]]$ts/max(conds[[i]]$ts))*4
}

sim_list<-llply(conds,.fun=function(x) ts(x$ts,start=c(1986,1),freq=12)) #Output object in list form

#Format list and dataframe outputs for further analysis (arrive here if inputting own empirical data)
#make list of timeseries "data_ls" (can input own empirical data here as list of timeseries)

data_ls <- sim_list #substitute "sim_list" for named list of empirical data
str(data_ls)

#make dataframe of timeseries "data_df"
individuals <- length(data_ls)
IDs <- names(data_ls)

data_df<-data.frame()
for(i in 1:individuals){
  d=data.frame(Date=seq(from=strptime(paste(start(data_ls[[i]])[1],start(data_ls[[i]])[2],01),"%Y %m %d"),to=strptime(paste(end(data_ls[[i]])[1],end(data_ls[[i]])[2],01),"%Y %m %d"), by="month"))
  d$Year=format(seq(from=strptime(paste(start(data_ls[[i]])[1],start(data_ls[[i]])[2],01),"%Y %m %d"),to=strptime(paste(end(data_ls[[i]])[1],end(data_ls[[i]])[2],01),"%Y %m %d"), by="month"),"%Y")
  d$Month=format(seq(from=strptime(paste(start(data_ls[[i]])[1],start(data_ls[[i]])[2],01),"%Y %m %d"),to=strptime(paste(end(data_ls[[i]])[1],end(data_ls[[i]])[2],01),"%Y %m %d"), by="month"),"%m")
  d$TS=data_ls[[i]][1:length(data_ls[[i]])]
  d$ID= factor(IDs[i])
  data_df<-rbind(data_df,d)
}

########################################################################
#### Initial analysis of phenological pattern at the species level #####
########################################################################

###To make raw data plots Figure 1a and 1b (main text)

##Figure 1a - circular calendar boxplots of phenology activity

#Summarise phenology scores across sample for each month and year
#(to give proportions for boxplots)
data_df_sample_summary<-ddply(data_df,.(Year,Month),summarize,
                              propPhenophase=length(TS[TS>0])/length(IDs))

data_df_sample_summary$Month<-mapvalues(data_df_sample_summary$Month,from=c("01","02","03","04","05","06","07","08","09","10","11","12"), to=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
data_df_sample_summary$Month<-factor(data_df_sample_summary$Month,ordered=T,levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

#Write circular boxplot graphic
Phenology_circular_boxplot<-ggplot(data_df_sample_summary,aes(x=factor(Month),y=propPhenophase))+
  geom_boxplot(fill="orange",outlier.colour="lightgrey")+
  coord_polar(start=0)+
  theme_minimal(base_size=16)+
  theme(legend.position="none",axis.title.x=element_blank(),axis.ticks=element_blank())+
  ylab("")

##Figure 1b - timeseries plot of raw data for each individual

raw_plots<-ggplot(data_df,aes(x=Date,y=TS,colour=ID))+
  geom_line()+
  scale_y_continuous(labels=NULL)+
  ylab("Canopy score")+
  xlab("")+
  theme_light(base_size=14)+
  theme(legend.position="none",strip.background = element_blank(), strip.text = element_blank())+
  facet_wrap(~ID,ncol=1)

Phenology_circular_boxplot
raw_plots

####################################################
#### Intitial observation of the Periodogram #######
###################################################

###To display smoothed periodograms - Figure 1c (main text)

##Parameters ##

# Spans for the Daniel kernel to smooth the spectrum The following
# list allows the spec_fun function to find appropriate spans for the Daniel
# kernel to successively apply to the spectrum to give a smoothed periodogram
# with bandwidth similar to 0.1 (span size is linked to the length of the
# oringal timeseries data, we give options here for data from 24 to 360 months
# long).
spans_lookup<-list(months=c(24,48,72,96,120,144,168,192,216,240,264,288,312,336,360),
                   spans_smooth=list(3,3,3,3,c(3,3),c(3,5),c(3,5),c(5,5),c(5,5),c(5,7),c(5,7),c(7,7),c(7,9),c(7,9),c(7,9)))

##Functions ##

#Spectrum fucntion
spec_fun<-function(x) spectrum(x,spans=spans_lookup$spans_smooth[[which.min(abs(spans_lookup$months-length(x))) ]],
                               plot=F,demean=T,detrend=T) #spectrum function for normal smoother periodogram


##Run fourier analysis using function "spectrum" on each list object (individual
##time series) to give Fourier outputs for each individual

Fourier_df<-data.frame()
for (i in 1:length(data_ls)){
  d<-data.frame(freq=spec_fun(data_ls[[i]])$freq/12)
  d$spec=spec_fun(data_ls[[i]])$spec
  d$spec_norm=spec_fun(data_ls[[i]])$spec*(1/mean(spec_fun(data_ls[[i]])$spec))
  d$ID=as.factor(IDs[i])
  Fourier_df<-rbind(Fourier_df,d)
}

##Plot periodograms for visual inspection

#Periogorams for each individual on one plot - Figure 1c
Periodograms<-ggplot(data=Fourier_df,aes(x=freq,y=spec_norm,group=ID,colour=ID)) +
  geom_line()+
  theme_minimal(base_size=14)+
  ylab("Power (standardised spectrum)")+
  xlab("Cycle frequency (cycles per month")+
  ggtitle("Smoothed periodograms of phenology cycles by individual")

Periodograms

#Summarise dominant peaks
Fourier_df_sample_summary<-ddply(Fourier_df, .(ID), summarize,
                                 maxFreq=freq[which.max(spec)],
                                 maxSpec_norm=max(spec_norm))

#Dominant phenology cycle for sample
Median_cycle<-1/median(Fourier_df_sample_summary$maxFreq)
Median_cycle

#####################################################################################
#### Periodogram analysis with confidence intervals to find frequency and phase #####
#####################################################################################

##Find frequency of dominant cycles and test against null hypothesis using 95%
##confidence interval derived from #Bloomfield et al. (2004)

##Parameters ##

#The following list allows the spec_fun function to find appropriate spans for
#the Daniel kernel to successively apply to the spectrum to give a smoothed
#periodogram with bandwidth similar to 0.1, and a super-smoothed spectrum (the
#null continuum) with bandwidth similar to 1. Span size is linked to the length
#of the oringal timeseries data, we give options here for data from 24 to 360
#months long.

spans_lookup<-list(months=c(24,48,72,96,120,144,168,192,216,240,264,288,312,336,360),
                   spans_smooth=list(3,3,3,3,c(3,3),c(3,5),c(3,5),c(5,5),c(5,5),c(5,7),c(5,7),c(7,7),c(7,9),c(7,9),c(7,9)),
                   spans_super_smooth=list(c(5,7),c(11,11),c(15,17),c(19,21),c(25,27),c(29,31),c(35,37),c(39,41),c(43,45),c(47,49),c(55,57),c(59,61),c(65,67),c(73,75),c(75,79)))

##Functions ##

#Function to calculate smoothed spectrum - > same function as before
spec_fun<-function(x) spectrum(x,spans=spans_lookup$spans_smooth[[which.min(abs(spans_lookup$months-length(x)))]],
                               plot=F,demean=T,detrend=T) #spectrum function for normal smoother periodogram

#Function to calculate the null continuum (super-smoothed spectrum of data) for null hypothesis test
spec_null_fun<-function(x) spectrum(x,spans=spans_lookup$spans_super_smooth[[which.min(abs(spans_lookup$months-length(x)))]],
                                    plot=F,demean=T,detrend=T) #spectrum function for null hypothesis spectrum

#Function to extract key variables from spectrum outputs
confidence_fun <- function(p){
  ts <- data_ls[[p]]
  d <- data.frame(ID=IDs[p])
  d$length <- length(ts)
  d$freq_dom <- (spec_fun(ts)$freq[which.max(spec_fun(ts)$spec)])/12 #frequency of the dominant peak
  d$cycle_dom <- 1/d$freq_dom
  d$spec_dom <- max(spec_fun(ts)$spec) #spectrum of dominant peak
  d$spec_dom_norm =  d$spec_dom*(1/mean(spec_fun(ts)$spec)) #normalised spectrum so that mean power is equal to 1 enabling comparisons of "Power" between individuals
  df <- spec_fun(ts)$df #degrees of freedom of the smoothed periodogram
  lower_ci <- (df*d$spec_dom)/(qchisq(c(0.975),df)) #lower CI of dominant peak of smoothed periodogram
  spec_null_dom <- spec_null_fun(ts)$spec[which(abs((spec_null_fun(ts)$freq)/12-d$freq_dom)==min(abs((spec_null_fun(ts)$freq)/12-d$freq_dom)))] #the corresponding dominant peak from the null continuum
  d$sig <- lower_ci>spec_null_dom #If the null spectrum is lower than the lower confidence interval, dominant peak can be considered "signficantly different" to null hypothesis
  #d$bw_smooth<-spec_fun(ts)$bandwidth #bandwidth of smoothed periodogram
  #d$bw_super_smooth<-spec_null_fun(ts)$bandwidth #bandwidth of super-smoothed periodogram
  data.frame(d)
}

## Run loop to apply spectrum functions to each individual time series in data_ls
confidence_results <- ldply(1:length(data_ls),.fun=confidence_fun)

#Categorise dominant cycles and screen results for positive false results where dominant peak is same as full length of the data - "No cyclicity"
confidence_results$cycle_category<-as.factor(ifelse(1/confidence_results$freq_dom > confidence_results$length/2,"No cyclicity",
                                                    ifelse(1/confidence_results$freq_dom > 13,"Supra-annual",
                                                           ifelse(1/confidence_results$freq_dom >= 11,"Annual","Sub-annual"))))

#Is the dominant cycle both significant and represent a regular cycle?
confidence_results$sig_cycle <- as.factor(ifelse(confidence_results$cycle_category!="No cyclicity" & confidence_results$sig==TRUE,TRUE,FALSE))

#Summary of Fourier outputs for sample
summary(confidence_results)

### Find phase of significant dominant cycles and assess synchrony between individuals

##Parameters ##

#Calculate the modal dominant cycle length of the sample
Dominant_cycle <- as.numeric(names(sort(-table(round(confidence_results$cycle_dom[confidence_results$sig_cycle==TRUE],digits=0))))[1])

#Simulate a "template" time series using the sample modal cycle peaking at the beginning of January 1986 to act as guide in co-Fourier analysis
w=2*pi*(1/Dominant_cycle) #wavelength in radians
simulated_phase_ts <- ts(4*cos(w*1:360),start=c(1986,2),freq=12)

##Functions ##

# Function for co-Fourier analysis of a simulated cosine curve against the empirical data
spec_cross_fun <- function(x)spec.pgram(x,spans=spans_lookup$spans_smooth[[which.min(abs(spans_lookup$months-length(x)))]],
                                        demean=T,detrend=T,plot=F)

#Function to extract data from co-Fourier anaysis regarding phase difference in radians
phase_fun<-function(p){
  ts_cross<-ts.intersect(simulated_phase_ts,data_ls[[p]])
  d<-data.frame(ID=IDs[p])
  d$sample_dominant_cycle<-Dominant_cycle
  #d$sim_freq_dom<- 1/12*spec_cross_fun(ts_cross)$freq [which.max(spec_cross_fun(ts_cross)$spec[,1])] #frequency of the dominant peak for the simulated data
  d$phase_radians<-spec_cross_fun(ts_cross)$phase[which.max(spec_cross_fun(ts_cross)$spec[,1])] #the phase difference (in radians) at the dominant frequency
  #d$coh<-spec_cross_fun(ts_cross)$coh[which.max(spec_cross_fun(ts_cross)$spec[,1])]
  d
}

## Run to apply the co-fourier spectrum function to each individual time series
# in combination with simulated guide time series for sample modal cycle to identify relative phase difference for each
phase_results<-ldply(1:length(data_ls),.fun=phase_fun)

##Combine confidence and phase results
results<-merge(confidence_results,phase_results,"ID")

#Convert phase (radians) to peak month (months since January 1st) if dominant cycle of individual data matches sample dominant cycle
results$peak_monthsSinceJan1st1986<- ifelse(round(results$cycle_dom)!=results$sample_dominant_cycle,NA,
                                            ifelse(results$sig_cycle==FALSE,NA,
                                                   ifelse(results$phase>0, results$phase/ (2*pi/results$sample_dominant_cycle),
                                                          (results$phase+2*pi)/ (2*pi/results$sample_dominant_cycle))))


#Convert months since Jan1st 1986 to the first peak calendar month
months_lookup<-data.frame(monthsSinceJan1st1986 = 0:11,CalendarMonth=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

results<-merge(results,adply(results,.margins=1,.expand=TRUE,.fun=summarise,
                             peak_firstCalendarMonth=months_lookup$CalendarMonth[which(months_lookup$monthsSinceJan1st1986== ifelse(round(peak_monthsSinceJan1st1986)>=12,round(peak_monthsSinceJan1st1986)-12,round(peak_monthsSinceJan1st1986)))])[,c(1,13)],
               "ID",all.x=T)

summary(results)

#Can also summarise these data for the whole sample:
summary_df<-data.frame(
  cycleLength_mean=round(mean(results$cycle_dom[results$sig_cycle==TRUE],na.rm=T),digits=1),
  cycleLength_sd=round(sd(results$cycle_dom[results$sig_cycle==TRUE],na.rm=T),digits=2),
  cycleLength_median=round(median(results$cycle_dom[results$sig_cycle==TRUE],na.rm=T),digits=1),
  cycleLength_mode=as.numeric(names(sort(-table(round(results[results$sig_cycle==TRUE,]$cycle_dom,digits=0))))[1]),
  phase_mean=round(mean.circular(results$phase[which(results$cycle_dom==results$sample_dominant_cycle&results$sig_cycle==T)],na.rm=T),digits=2),
  phase_sd=round(sd.circular(results$phase[which(results$cycle_dom==results$sample_dominant_cycle&results$sig_cycle==T)],na.rm=T),digits=2)) #level of synchronicity between individuals
summary_df$peak_monthsSinceJan1st1986_mean<-ifelse(summary_df$phase_mean<0,(summary_df$phase_mean+2*pi)/(2*pi/summary_df$cycleLength_median),(summary_df$phase_mean)/(2*pi/summary_df$cycleLength_median))
summary_df$peak_monthsSinceJan1st1986_sd<-summary_df$phase_sd/(2*pi/summary_df$cycleLength_median)
summary_df$peak_firstCalendarMonth_mean<-months_lookup$CalendarMonth[which(months_lookup$monthsSinceJan1st1986== ifelse(round(summary_df$peak_monthsSinceJan1st1986_mean)>=12,round(summary_df$peak_monthsSinceJan1st1986_mean)-12,round(summary_df$peak_monthsSinceJan1st1986_mean)))]

summary_df

#########################################
#### Power analyses - simulated data ####
#########################################

###Power anaysis of time series length and noise to detect periodicity
###(simulated data) - FIGURE 4 (code for testing against both plausable null
###hypotheses as described in supporting information S2)

#Functions

#Smoothed spectrum
spec_fun<-    function(x) spectrum(x,spans=spans_comb[[which.min(abs(months_comb-length(x)))]] ,plot=F,demean=T,detrend=T) #spectrum function for normal smoother periodogram
#Super-smoothed spectrum for null continuum null hypothesis comparison
spec_null_fun<-function(x) spectrum(x,spans=spans_super_smooth_comb[[which.min(abs(months_comb-length(x)))]] ,plot=F,demean=T,detrend=T) #spectrum function for null hypothesis spectrum

##Spans
#The following lists give appropriate successive spans for the corresponding
#length of data wto give a smoothed periodogram with bandwidth in magnitude of
#0.1 (spans_comb) and a null continuum smoother periodogram with bandwidth in
#magnitude of 1 (spans_super_smooth_comb)
months_comb<-c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30)*12
spans_comb<-list(3,3,3,3,c(3,3),c(3,5),c(3,5),c(5,5),c(5,5),c(5,7),c(5,7),c(7,7),c(7,9),c(7,9),c(7,9))
spans_super_smooth_comb<-list(c(5,7),c(11,11),c(15,17),c(19,21),c(25,27),c(29,31),c(35,37),c(39,41),c(43,45),c(47,49),c(55,57),c(59,61),c(65,67),c(73,75),c(75,79))#spans_smooth_fun<-function(x) ifelse(odd(round(sqrt(sqrt(x))))==TRUE,round(sqrt(sqrt(x))),round(sqrt(sqrt(x)))+1)

#Extracting key variables from spectrum outputs
test_fun<-function(ts) {
  start<-c(sample(start(ts)[1]:(end(ts)[1]-window_trial),1),sample(1:12,1))
  end<- c( ifelse(start[2]>1,start[1]+window_trial,start[1]+(window_trial-1)), ifelse(start[2]>1,start[2]-1,12))
  ts<-window(ts,start=start,end=end) #cut window from simulated data
  d<-list(sd=sd_trial) #record regularity of peak (SD of truncated normal distribution)
  d$replacement=replacement_trial #record event detectability - percentage positive scores replaced with zeros
  d$window<-window_trial #length of window
  d$months<-length(ts) #check length of timeseries against window length
  d$freq_max<-spec_fun(ts)$freq[which.max(spec_fun(ts)$spec)]/12 #frequency of dominant peak
  d$spec<-max(spec_fun(ts)$spec) #spectrum of dominant peak
  d$df<-spec_fun(ts)$df #degrees of freedom for CI calculations
  d$lower_ci<-(d$df*d$spec)/(qchisq(c(0.975),d$df)) #lower CI
  d$null<-spec_null_fun(ts)$freq/12
  d$null<-spec_null_fun(ts)$spec[which.min(abs(d$null-d$freq_max))]   #Corresponding peak in null continuum
  d$sig_null<-d$lower_ci>d$null #Significant compared to null continuum?
  d$white<-sum(spec_fun(ts)$spec)/length(spec_fun(ts)$spec) #Corresponding peak in white spectrum (for comparison as null hypothesis)
  d$sig_white<-d$lower_ci>d$white #Significant compared to white noise?
  d$dom_cycle<-d$freq_max>1/13.5&d$freq_max<1/10.5 #Does dominant peak fall within expected frequency interval?
  d$cyclicity<-as.factor(ifelse(1/d$freq_max>d$months/2,"No cyclicity","Cyclicity")) #Screen for false positive results when dominant peak is same as full length of time series "No cyclicity"
  d<-data.frame(d)
}

# Set fundamental variables for the simulated time series sample. The variables
#below are identical for each  individual, but these can also be varied by individual
peak_mean<-6
peak_sd<-c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6)
replacement<-c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8)
windows<-c(5,10,15)

years<-35
wavelength<-12
trials<-10000

peak_scores<-c(0.5,1,2,2,2,3,3,3,4,4,4,4) #Similar to distribution of positive scores in flowering data from Lope long-term study
mid_peak_scores<-c(0.1,0.1,0.1,0.1,0.1,0.5,0.5,0.5,0.5,1,1,1,1,1.5,1.5,2,3,4)

#Loop to simulate data and analyse using Fourier to extract information about
#dominant cycles under different noise and data length scenarios
results_trials<-data.frame()
for (t in 1:trials){
  ts<-vector()
  sd_trial<-sample(peak_sd,1)
  replacement_trial<-sample(replacement,1)
  window_trial<-sample(windows,1)
  for(y in (1:years)){
    peak_month<-round(rtruncnorm(1,a=2,b=wavelength-1,mean=peak_mean,sd=sd_trial),0)
    ts_year<-rep(0,length.out=wavelength)
    ts_year[[peak_month]] <-sample(peak_scores,1)
    ts_year[[peak_month+1]] <-sample(mid_peak_scores,1)
    ts_year[[peak_month-1]] <-sample(mid_peak_scores,1)
    ts<-c(ts,ts_year)
  }
  ts[ts>0][sample(length(ts[ts>0]),round(replacement_trial*length(ts[ts>0]),0),replace=F)]<-  0
  ts<-ts(ts,start=c(1986,1),freq=12)
  results_trials<-rbind(results_trials,test_fun(ts))
}

summary(results_trials)

#Determine if estimates of dominant peak represent a real cycle, are significant and fall within expected frequency interval
#Compared to null continuum
results_trials$null_dom<-ifelse(results_trials$cyclicity=="Cyclicity"&results_trials$sig_null==TRUE&results_trials$dom_cycle==T,1,0)
results_trials$null_other<-ifelse(results_trials$cyclicity=="Cyclicity"&results_trials$sig_null==TRUE&results_trials$dom_cycle==F,1,0)
#Compared to white noise
results_trials$white_dom<-ifelse(results_trials$cyclicity=="Cyclicity"&results_trials$sig_white==TRUE&results_trials$dom_cycle==T,1,0)
results_trials$white_other<-ifelse(results_trials$cyclicity=="Cyclicity"&results_trials$sig_white==TRUE&results_trials$dom_cycle==F,1,0)

#Inspect for each window length 5,10 and 15
results_trials1<-results_trials[which(results_trials$window==5),]

#Summarise for each level of noise
results_props<-ddply(results_trials1,.(sd,replacement),summarize,
                     prop_null_dom=sum(null_dom)/length(null_dom),
                     prop_null_other=sum(null_other)/length(null_other),
                     prop_white_dom=sum(white_dom)/length(white_dom),
                     prop_white_other=sum(white_other)/length(white_other))

#Creat matrix for plotting
matrix_data<-rbind(data.frame(results_props[,c(1:2)],proportion=results_props[,3],hypothesis="null",peak="dominant"),
                   data.frame(results_props[,c(1:2)],proportion=results_props[,4],hypothesis="null",peak="other"),
                   data.frame(results_props[,c(1:2)],proportion=results_props[,5],hypothesis="white",peak="dominant"),
                   data.frame(results_props[,c(1:2)],proportion=results_props[,6],hypothesis="white",peak="other"))

matrix_pp<-list(null_dominant<-subset(matrix_data,hypothesis=="null"&peak=="dominant"),null_other<-subset(matrix_data,hypothesis=="null"&peak=="other"),white_dominant<-subset(matrix_data,hypothesis=="white"&peak=="dominant"),white_other<-subset(matrix_data,hypothesis=="white"&peak=="other"))
zz<-llply(matrix_pp,.fun=function(x) matrix(x$proportion,ncol=length(unique(x$sd)),nrow=length(unique(x$replacement)),byrow=F))
xx <- unique(matrix_data$replacement)
yy <- unique(matrix_data$sd)

set.seed(1)
par(mfrow=c(1,1))

colfunc<-colorRampPalette(c("#ffffd9","#41b6c4","#253494","#081d58"))

image(xx,yy,z=zz[[1]],col=colfunc(100),zlim=range(unlist(zz),image.plot(xx,yy,z=zz[[1]],col=colfunc(100),zlim=range(unlist(zz)))))
image(xx,yy,z=zz[[2]],col=colfunc(100),zlim=range(unlist(zz),image.plot(xx,yy,z=zz[[2]],col=colfunc(100),zlim=range(unlist(zz)))))
image(xx,yy,z=zz[[3]],col=colfunc(100),zlim=range(unlist(zz),image.plot(xx,yy,z=zz[[3]],col=colfunc(100),zlim=range(unlist(zz)))))
image(xx,yy,z=zz[[4]],col=colfunc(100),zlim=range(unlist(zz),image.plot(xx,yy,z=zz[[4]],col=colfunc(100),zlim=range(unlist(zz)))))
