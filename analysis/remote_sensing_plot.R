library(tidyverse)
library(ggplot2)
library(MODISTools)
library(ggthemes)
library(gridExtra)

# load sites
sites <- read.table("data-raw/yangambi_sites.csv",
                    header = TRUE,
                    sep = ",")

# download and/or read data
if(!file.exists("data-raw/yangambi_mod13q1_evi_values.csv")){
  message("downloading VI data")
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
  message("reading VI data")
  VI <- read.table("data-raw/yangambi_mod13q1_evi_values.csv",
                   header = TRUE,
                   sep = ",")
}

if(!file.exists("data-raw/yangambi_mod13q1_nir_values.csv")){
  message("downloading NIR data")
  NIR <- mt_batch_subset(sites,
                        product = "MOD13Q1",
                        start = "2001-01-01",
                        end = "2018-12-31",
                        band = "250m_16_days_NIR_reflectance",
                        internal = TRUE)

  write.table(NIR, "data-raw/yangambi_mod13q1_nir_values.csv",
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE,
              sep = ",")
} else {
  message("reading NIR data")
  NIR <- read.table("data-raw/yangambi_mod13q1_nir_values.csv",
                   header = TRUE,
                   sep = ",")
}

if(!file.exists("data-raw/yangambi_mod13q1_qa_values.csv")){
  message("downloading QA data")
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
  message("reading QA data")
  QA <- read.table("data-raw/yangambi_mod13q1_qa_values.csv",
                   header = TRUE,
                   sep = ",")
}

if(!file.exists("data-raw/yangambi_mod17a2h_gpp_values.csv")){
  message("downloading GPP data")
  GPP <- mt_batch_subset(sites,
                        product = "MOD17A2H",
                        band = "Gpp_500m",
                        start = "2001-01-01",
                        end = "2018-12-31",
                        internal = TRUE)

  write.table(GPP, "data-raw/yangambi_mod17a2h_gpp_values.csv",
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE,
              sep = ",")

} else {
  message("reading GPP data")
  GPP <- read.table("data-raw/yangambi_mod17a2h_gpp_values.csv",
                   header = TRUE,
                   sep = ",")
}

if(!file.exists("data-raw/yangambi_mcd15a3h_fapar_values.csv")){
  message("downloading FAPAR data")
  FAPAR <- mt_batch_subset(sites,
                         product = "MCD15A3H",
                         band = "Fpar_500m",
                         start = "2001-01-01",
                         end = "2018-12-31",
                         internal = TRUE)

  write.table(FAPAR, "data-raw/yangambi_mcd15a3h_fapar_values.csv",
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE,
              sep = ",")

} else {
  message("reading FAPAR data")
  FAPAR <- read.table("data-raw/yangambi_mcd15a3h_fapar_values.csv",
                    header = TRUE,
                    sep = ",")
}


# screen for QA
VI$value[QA$value > 1] <- NA

# create mean values by DOY
FAPAR_s <- FAPAR %>%
  mutate(value = value * as.numeric(scale),
         doy = as.numeric(format(as.Date(calendar_date), "%j")),
         site = toupper(site)) %>%
  group_by(site, doy) %>%
  summarize(FAPAR = mean(value, na.rm = TRUE))

# create mean values by DOY
VI_s <- VI %>%
  mutate(value = value * as.numeric(scale),
         doy = as.numeric(format(as.Date(calendar_date), "%j")),
         site = toupper(site)) %>%
  group_by(site, doy) %>%
  summarize(EVI = mean(value, na.rm = TRUE))

# create mean values by DOY
GPP_s <- GPP %>%
  mutate(value = value * 0.0125,
         doy = as.numeric(format(as.Date(calendar_date), "%j")),
         site = toupper(site)) %>%
  group_by(site, doy) %>%
  summarize(GPP = mean(value, na.rm = TRUE))


# plot EVI by site
p1 <- ggplot(FAPAR_s, aes(doy, FAPAR)) +
  geom_point(col = "black") +
  geom_smooth(method = "loess", span = 0.3, col = "black") +
  ylab("MCD15A3H FAPAR (%)") +
  xlab("") +
  theme_bw()

# plot EVI by site
p2 <- ggplot(VI_s, aes(doy, EVI)) +
  geom_point(col = "black") +
  geom_smooth(method = "loess", span = 0.3, col = "black") +
  ylab("MCD13Q1 EVI") +
  xlab("DOY") +
  theme_bw()

# plot EVI by site
p3 <- ggplot(GPP_s, aes(doy, GPP)) +
  geom_point(col = "black") +
  geom_smooth(method = "loess", span = 0.3, col = "black") +
  ylab(bquote("MOD17A2H GPP (g C/" ~ m^2 ~ "d)")) +
  xlab("") +
  theme_bw()

gridExtra::grid.arrange(p1, p2, p3)

# with sofun
gridExtra::grid.arrange(p1, p2, p3, nrow = 3)

