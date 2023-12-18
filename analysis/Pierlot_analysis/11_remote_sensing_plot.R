library(tidyverse)
library(ggplot2)
#library(MODISTools)
library(ggthemes)
library(showtext)
font_add_google(
  "Lato",
  regular.wt = 300,
  bold.wt = 700)

# load sites
sites <- read.table("data-raw/yangambi_sites.csv",
                    header = TRUE,
                    sep = ",")

# download and/or read data
if(!file.exists("data-raw/yangambi_mod13q1_evi_values.csv")){
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

if(!file.exists("data-raw/yangambi_mod13q1_qa_values.csv")){
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
  QA <- read.table("data-raw/yangambi_mod13q1_qa_values.csv",
                   header = TRUE,
                   sep = ",")
}

# screen for QA
VI$value[QA$value > 1] <- NA

# create mean values by DOY
VI_s <- VI %>%
  dplyr::mutate(value = value * as.numeric(scale),
                doy = as.numeric(format(as.Date(calendar_date), "%j")),
                site = toupper(site)) %>%
  group_by(site, doy) %>%
  dplyr::summarise(EVI = mean(value, na.rm = TRUE),
                   EVI_sd = sd(value, na.rm = TRUE))
VI_conf <- VI_s %>%
  group_by(doy) %>%
  dplyr::summarise(EVImax = max(EVI, na.rm = TRUE) + max(EVI_sd),
                   EVImin = min(EVI, na.rm = TRUE) - max(EVI_sd),
                   EVImean = mean(EVI, na.rm = TRUE))

# plot EVI by site
p_modis <- ggplot(VI_s) +
  annotate("rect", xmin = 0, xmax = 61, ymin = 0.25, ymax = 0.77, alpha = .1) + # jan - febr
  annotate("rect", xmin = 152.5, xmax = 213.5, ymin = 0.25, ymax = 0.77, alpha = .1) + # jun - jul
  annotate("rect", xmin = 335.5, xmax = 365, ymin = 0.25, ymax = 0.77, alpha = .1) + # dec
  annotate("text", x = 58, y = 0.72, label = "EVI-MOD13Q1\n2000 – 2018") +
  # geom_point(aes(doy, EVI, shape = site), col = "grey40") +
  geom_ribbon(data = VI_conf, aes(x = doy, ymin = EVImin, ymax = EVImax), fill="grey60", alpha=0.5) +
  geom_point(data = VI_conf, aes(doy, EVImean), col = "grey30") +
  geom_smooth(aes(doy, EVI), span = 0.3, se = FALSE, col = "grey30", size = 1.2) +
  #geom_line(aes(doy, EVI + EVI_sd)) +
  #geom_line(aes(doy, EVI - EVI_sd)) +
  labs(#title = "MOD13Q1",
    #subtitle = "mean ...",
    x = "",
    y = "EVI") + # - MOD13Q1
  scale_x_continuous(limits = c(0,365),
                     breaks = seq(0,365,30.5),
                     labels = month.abb) +
  # scale_y_continuous(limits = c(0.4,0.6),
  #                    breaks = c(0.4,0.5,0.6),
  #                    labels = scales::number_format(accuracy = 0.1)) +
  theme_minimal() +
  theme(text = element_text(family = "Lato", color = "#22211d"),
        panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        # axis.text.x=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 2,size = 10), # vjust to center the label
        axis.title.x=element_blank(),
        axis.title.y = element_text(vjust = 3),
        legend.position = "none",
        plot.margin=unit(c(0,0,0,0.5),"cm")
  )

# print(p_modis)

ggsave("manuscript/leaf_phenology/figures/Pierlot_fig6.png", p_modis,
       device = "png", width = 4.5, height = 3)
