library(tidyverse)
library(ggplot2)
library(MODISTools)
library(ggthemes)

# load sites
sites <- read.table("data-raw/yangambi_sites.csv",
                    header = TRUE,
                    sep = ",")

# download and/or read data
if(file.exists("data-raw/yangambi_mod13q1_evi_values.csv")){
  VI <- mt_batch_subset(sites,
                product = "MOD13Q1",
                start = "2001-01-01",
                end = "2018-12-31",
                band = "250m_16_days_EVI",
                internal = TRUE)

  write.table(VI, "data-raw/yangambi_mod13q1_evi_values.csv",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = ",")
} else {
  VI <- read.table("data-raw/yangambi_mod13q1_evi_values.csv",
                   header = TRUE,
                   sep = ",")
}

if(file.exists("data-raw/yangambi_mod13q1_qa_values.csv")){
  QA <- mt_batch_subset(sites,
                      product = "MOD13Q1",
                      band = "250m_16_days_pixel_reliability",
                      start = "2001-01-01",
                      end = "2018-12-31",
                      internal = TRUE)

  write.table(QA, "data-raw/yangambi_mod13q1_qa_values.csv",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = ",")

} else {
  VI <- read.table("data-raw/yangambi_mod13q1_evi_values.csv",
                   header = TRUE,
                   sep = ",")
}

# screen for QA
VI$value[QA$value > 1] <- NA

# create mean values by DOY
VI_s <- VI %>%
  mutate(value = value * as.numeric(scale),
         doy = as.numeric(format(as.Date(calendar_date), "%j")),
         site = toupper(site)) %>%
  group_by(site, doy) %>%
  summarize(EVI = mean(value, na.rm = TRUE))

# plot EVI by site
p <- ggplot(VI_s, aes(doy, EVI)) +
  geom_point(col = "green") +
  geom_smooth(method = "loess", span = 0.3) +
  facet_wrap(~site)

print(p)

