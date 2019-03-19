# site years
library(tidyverse)

precip <- readxl::read_xlsx("~/Downloads/Données_climatiques_INERA_Yangambi_Km5_1931_2010_monthly.xlsx", sheet = 1)

tmin <- readxl::read_xlsx("~/Downloads/Données_climatiques_INERA_Yangambi_Km5_1931_2010_monthly.xlsx", sheet = 6)

tmax <- readxl::read_xlsx("~/Downloads/Données_climatiques_INERA_Yangambi_Km5_1931_2010_monthly.xlsx", sheet = 5)

tmean <- (tmax + tmin) / 2

sun_hours <- readxl::read_xlsx("~/Downloads/Données_climatiques_INERA_Yangambi_Km5_1931_2010_monthly.xlsx", sheet = 4)

tmean <- tmean %>% gather(key = "month", value = "value", -Année) %>%
 rename(year = "Année") %>%
  mutate(date = as.Date(sprintf("%s-%s-15",year,month),"%Y-%m-%d"),
         measurement = "t_mean")

precip <- precip %>% gather(key = "month", value = "value", -Année) %>%
  rename(year = "Année") %>%
  mutate(date = as.Date(sprintf("%s-%s-15",year,month),"%Y-%m-%d"),
         measurement = "precip")

sun_hours <- sun_hours %>% gather(key = "month", value = "value", -Année) %>%
  rename(year = "Année") %>%
  mutate(date = as.Date(sprintf("%s-%s-15",year,month),"%Y-%m-%d"),
         measurement = "sun_hours")

climate_data <- rbind(tmean,precip, sun_hours)
climate_data$month <- as.numeric(climate_data$month)

saveRDS(climate_data, "data/yangambi_monthly_climate.rds")

climate <- climate_data %>% group_by(measurement, month) %>%
  summarize(value = mean(value))

p <- ggplot() + geom_line(data = climate,
                          aes(x = month,
                y = value,
                group = measurement,
                colour = measurement))
print(p)
