library(seewave)

# read climate data
climate <- read.csv("data/yangambi_km5_monthly_kasongo.csv")
climate <- climate %>%
  filter(year > 1936,
         year < 1959) %>%
  # mutate(col = ifelse(precip < drought_mm,"dry","wet"))
  mutate(col = ifelse(month %in% c(1,2,6,7,12), "dry", "wet"))
climate$date <- paste(climate$year,climate$month,"15",sep = "-")
climate$date <- as.Date(climate$date, "%Y-%m-%d")
climate$yr_month <- format(as.Date(climate$date), "%Y-%m")

climate.avg <- read.csv("data/ClimData_monthly_avg.csv",header = TRUE,sep = ",")
climate.avg = climate.avg[,(names(climate.avg) %in% c("Month","insol_JR","tmax_JR"))]
colnames(climate.avg)[1] <- "month"
climate <- merge(climate, climate.avg, by = "month", all.x = TRUE)

# climate time series
climate <- climate %>%
  dplyr::arrange(date)
first_year <- as.numeric(climate$year[1])
climate_insol_ts <- ts(data = climate$insol_JR, start=c(1938, 1), end=c(1957, 12), frequency = 12) # frequency = 48)
str(climate_insol_ts)
#----------------------------------------------------------------------
#----------------------------------------------------------------------

data_subset <- data %>%
    filter(species_full %in% "Albizia adianthifolia") %>% #Pericopsis elata #Scorodophloeus zenkeri
    filter(phenophase %in% "leaf_dormancy")

# convert date
data_subset$date <- as.Date(paste(data_subset$year,
                                  round((data_subset$week*7.6)-7),sep="-"), "%Y-%j")
data_subset$yr_month <- format(data_subset$date, "%Y-%m")
data_subset$yr <- format(data_subset$date, "%Y")


# average by date
data_subset <- data_subset %>%
  group_by(species_full, yr_month) %>% # date
  dplyr::summarise(mean_value = mean(value),
                   scaled_value = ifelse(any(value > 0), 1, 0), # any observation: value 1
                   sum_value = sum(value),
                   nr_indiv = length(value),
                   yr = min(yr))

# sort dataframe according to date
data_subset <- data_subset %>%
  dplyr::arrange(yr_month)
# # add a month column
# data_subset$month <- format.Date(data_subset$yr_month, "%m")
# data_subset$year <- format.Date(data_subset$yr_month, "%Y")

# data as timeseries, species_level
# first_year <- as.numeric(format.Date(data_subset$date[1], "%Y"))
first_year <- as.numeric(data_subset$yr[1])
data_ts <- ts(data = data_subset$mean_value, start = first_year, frequency = 12) # frequency = 48)
str(data_ts)

#-------------------------------------------------
# coherence <- ccoh(data_ts, climate_insol_ts, channel = c(1,1), wl = 24, ovlp = 1)
coherence <- coh(data_ts, climate_insol_ts, f = 12) #, xlim = c(0,12/1000))
lags      <- coherence[,1]
vals      <- coherence[,2]
val.max   <- vals[which.max(vals)]
lag.max   <- lags[which.max(vals)]

freq_dom = (spec_fun(data_ts)$freq[which.max(spec_fun(data_ts)$spec)])/48 #frequency of the dominant peak
cycle_dom = 1/freq_dom


##########################
# find cff max lag = cross correlation
##########################
data_sp <- data %>%
  filter(species_full %in% "Pericopsis elata",
         phenophase == "leaf_dormancy")
# convert date to year-month --> to be able to match with monthly climate data
data_sp$date <- as.Date(
  paste(data_sp$year,
        round(data_sp$week*7.6),sep="-"), "%Y-%j")
data_sp$yr_month <- format(as.Date(data_sp$date), "%Y-%m")
# phenological data to monthly resolution
data_sp <- data_sp %>%
  group_by(yr_month) %>%
  dplyr::summarise(mean_value = mean(value),
                   sum_value = sum(value),
                   nr_indiv = length(value))
# remove if data = NA
data_sp <- data_sp[!(is.na(data_sp$yr_month)),]
# merge with climate data
data_sp <- merge(data_sp, climate, by = "yr_month", all.x = TRUE)
# cross-correlation + confidence interval
corr.precip <- ccf(data_sp$precip, data_sp$mean_value,lag = 6, pl = FALSE)
plot(corr.precip)
##########################


##########################
# fourier transform
##########################
spec_fun <- function(x) spectrum(x,#spans = 10,
                                 plot=T,demean=T,detrend=T) #spectrum function for normal smoother periodogra

ftest <- data_sp %>%
  # group_by(species_full) %>%
  do(freq_list = matrix(spec_fun(.$mean_value)$spec)) %>%
  unnest(freq_list)

ftestclimate <- data_sp %>%
  # group_by(species_full) %>%
  do(freq_list = matrix(spec_fun(.$mean_value)$spec)) %>%
  unnest(freq_list)
