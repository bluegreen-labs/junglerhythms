#----- reset your R session. -----------------------------------------
rm(list=ls())
#----- load required libraries ---------------------------------------
library(plyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(gtools)
library(grid)
library(gridExtra)
library(tseries)
library(showtext)
font_add_google(
  "Lato",
  regular.wt = 300,
  bold.wt = 700)
#----- source required functions --------------------------------------
source("R/climate_ccf_function.R")
#----------------------------------------------------------------------
# Suppress summarise info
# because `summarise()` ungrouping output (override with `.groups` argument)
# comes up a lot in the consol, which is annoying
options(dplyr.summarise.inform = FALSE)


#----------------------------------------------------------------------
#--------   Phenology data --------------------------------------------
#----------------------------------------------------------------------
data <- readRDS("data/jungle_rhythms_data_Pierlot_species-96.rds")
species_list <- unique(data$species_full)
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#--- Climate data -----------------------------------------------------
#----------------------------------------------------------------------
climate <- read.csv("data/yangambi_km5_monthly_kasongo.csv")
climate$date <- paste(climate$year,climate$month,"15",sep = "-")
climate$date <- as.Date(climate$date, "%Y-%m-%d")
climate$date_monthly <- format(as.Date(climate$date), "%Y-%m")

climate.avg <- read.csv("data/ClimData_monthly_avg.csv",header = TRUE,sep = ",")
climate.avg = climate.avg[,(names(climate.avg) %in% c("Month","insol_JR","tmax_JR"))]
colnames(climate.avg)[1] <- "month"
climate <- merge(climate, climate.avg, by = "month", all.x = TRUE)
climate <- climate %>%
  dplyr::arrange(date)


#----------------------------------------------------------------------
#-------- get some metadata on the species ----------------------------
#----------------------------------------------------------------------
overview <- read.csv("data/Pierlot_summ_species_pheno_characteristics.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)
# get extra info on species from 'overview' file
# in summary figure, only species with events are included
overview <- overview %>%
  dplyr::select(species_full,
                phenophase,
                deciduousness,
                total_nr_events) # species with too few events (< 5) not included for cross correlation analysis

#--- leaf turnover ---------------------------------------------------
# species-specific timelines, average across individuals
timelines_sp_turn <- data %>%
  filter(phenophase == "leaf_turnover") %>%
  group_by(species_full, date) %>%
  dplyr::summarise(mean_value = mean(value, na.rm = TRUE))
timelines_sp_turn$mean_value[is.nan(timelines_sp_turn$mean_value)] <- NA
timelines_sp_turn$scaled_value <- ifelse(timelines_sp_turn$mean_value > 0, 1,
                                         ifelse(timelines_sp_turn$mean_value == 0, 0, NA))
timelines_sp_turn$phenophase <- "leaf_turnover"
timelines_sp_turn <- merge(timelines_sp_turn, overview, by = c("species_full","phenophase"), all.x = TRUE)

#--- leaf dormancy ----------------------------------------------------
timelines_sp_dorm <- data %>%
  filter(phenophase == "leaf_dormancy") %>%
  group_by(species_full, date) %>%
  dplyr::summarise(mean_value = mean(value, na.rm = TRUE))
timelines_sp_dorm$mean_value[is.nan(timelines_sp_dorm$mean_value)] <- NA
timelines_sp_dorm$scaled_value <- ifelse(timelines_sp_dorm$mean_value > 0, 1,
                                         ifelse(timelines_sp_dorm$mean_value == 0, 0, NA))
timelines_sp_dorm$phenophase <- "leaf_dormancy"
timelines_sp_dorm <- merge(timelines_sp_dorm, overview, by = c("species_full","phenophase"), all.x = TRUE)
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#--- cross correlations climate - phenology timeseries ----------------
#--- phenology timeseries with too few events (< 5) are ---------------
#---not analysed within the function (set to NA) ----------------------
#----------------------------------------------------------------------
# precip
crosscorr_precip_turn <- climate_ccf(data = timelines_sp_turn,
                                     climate = climate,
                                     climate.variable = "precip",
                                     species_name = species_list,
                                     pheno = "leaf_turnover")
crosscorr_precip_dorm <- climate_ccf(data = timelines_sp_dorm,
                                     climate = climate,
                                     climate.variable = "precip",
                                     species_name = species_list,
                                     pheno = "leaf_dormancy")
crosscorr_precip <- rbind(crosscorr_precip_turn, crosscorr_precip_dorm)
# sun
crosscorr_sun_turn <- climate_ccf(data = timelines_sp_turn,
                                  climate = climate,
                                  climate.variable = "sun",
                                  species_name = species_list,
                                  pheno = "leaf_turnover")
crosscorr_sun_dorm <- climate_ccf(data = timelines_sp_dorm,
                                  climate = climate,
                                  climate.variable = "sun",
                                  species_name = species_list,
                                  pheno = "leaf_dormancy")
crosscorr_sun <- rbind(crosscorr_sun_turn, crosscorr_sun_dorm)

# temp
crosscorr_temp_turn <- climate_ccf(data = timelines_sp_turn,
                                   climate = climate,
                                   climate.variable = "temp",
                                   species_name = species_list,
                                   pheno = "leaf_turnover")
crosscorr_temp_dorm <- climate_ccf(data = timelines_sp_dorm,
                                   climate = climate,
                                   climate.variable = "temp",
                                   species_name = species_list,
                                   pheno = "leaf_dormancy")
crosscorr_temp <- rbind(crosscorr_temp_turn, crosscorr_temp_dorm)

# merge
# output <- merge(overview[1], dormancy_precip, by = "species_full", all.x = TRUE)
output <- merge(crosscorr_precip, crosscorr_sun, by = c("species_full","phenophase"), all.x = TRUE)
output <- merge(output, crosscorr_temp, by = c("species_full","phenophase"), all.x = TRUE)

#----------------------------------------------------------------------
#--- analysis showed 2 specific groups based on cross-correlaitons ----
#--- group 1: species with a negative correlation to precipitation
#             and a positive in-phase correlation to sun hours and
#             maximum temperature;
#--- group 2: species with a negative correlation to sun hours and
#             maximum temperature
#----------------------------------------------------------------------

# set.seed(1)
# test <- output %>%
#   filter(phenophase == "leaf_turnover") %>%
#   select(!phenophase)
# test$corr_precip_timing <- as.numeric(test$corr_precip_timing)
# test$corr_sun_timing <- as.numeric(test$corr_sun_timing)
# test$corr_temp_timing <- as.numeric(test$corr_temp_timing)
# # test$species_full <- ifelse(is.na(test$corr_precip) & is.na(test$corr_sun) & is.na(test$corr_temp), NA, test$species_full)
# # test <- test %>%
# #   filter(!is.na(species_full))
#
#
# test$corr_precip <- ifelse(is.na(test$corr_precip), 0, test$corr_precip)
# test$corr_sun <- ifelse(is.na(test$corr_sun), 0, test$corr_sun)
# test$corr_temp <- ifelse(is.na(test$corr_temp), 0, test$corr_temp)
# test$corr_precip_timing <- ifelse(is.na(test$corr_precip_timing), 10, test$corr_precip_timing)
# test$corr_sun_timing <- ifelse(is.na(test$corr_sun_timing), 10, test$corr_sun_timing)
# test$corr_temp_timing <- ifelse(is.na(test$corr_temp_timing), 10, test$corr_temp_timing)
#
# rownames(test) <- test$species_full
# test <- test %>%
#   select(!species_full)
# # test <- na.omit(test)
# # Determine number of clusters
# wss <- (nrow(test)-1)*sum(apply(test,2,var))
# for (i in 2:15) wss[i] <- sum(kmeans(test,
#                                      centers=i)$withinss)
# plot(1:15, wss, type="b", xlab="Number of Clusters",
#      ylab="Within groups sum of squares")
#
#
# fit <- kmeans(test, 3)
# # get cluster means
# aggregate(test, by=list(fit$cluster), FUN=mean)
# # append cluster assignment
# test <- data.frame(test, fit$cluster)

# is.na didn't work in the big ifelse statement
output$corr_precip <- ifelse(is.na(output$corr_precip), 0, output$corr_precip)
output$corr_sun <- ifelse(is.na(output$corr_sun), 0, output$corr_sun)
output$corr_temp <- ifelse(is.na(output$corr_temp), 0, output$corr_temp)


output$groups <- ifelse(output$corr_precip < 0 &  output$corr_sun > 0 &    output$corr_temp > 0, "group1",
                        ifelse(output$corr_precip == 0 & output$corr_sun > 0 &    output$corr_temp > 0, "group1",
                               ifelse(output$corr_precip == 0 & output$corr_sun > 0 &    output$corr_temp == 0, "group1",
                                      ifelse(output$corr_precip == 0 & output$corr_sun == 0 &   output$corr_temp > 0, "group1",
                                             ifelse(output$corr_precip < 0 &  output$corr_sun == 0 &   output$corr_temp > 0, "group1",
                                                    ifelse(output$corr_precip < 0 &  output$corr_sun > 0 &    output$corr_temp == 0, "group1",
                                                           ifelse(output$corr_precip < 0 &  output$corr_sun == 0 &   output$corr_temp == 0, "group1",

                                                                  ifelse(output$corr_precip > 0 &  output$corr_sun < 0 &    output$corr_temp < 0, "group2",
                                                                         ifelse(output$corr_precip == 0 & output$corr_sun < 0 &    output$corr_temp < 0, "group2",
                                                                                ifelse(output$corr_precip == 0 & output$corr_sun < 0 &    output$corr_temp == 0, "group2",
                                                                                       ifelse(output$corr_precip == 0 & output$corr_sun == 0 &   output$corr_temp < 0, "group2",
                                                                                              ifelse(output$corr_precip > 0 &  output$corr_sun == 0 &   output$corr_temp < 0, "group2",
                                                                                                     ifelse(output$corr_precip > 0 &  output$corr_sun < 0 &    output$corr_temp == 0, "group2",
                                                                                                            ifelse(output$corr_precip > 0 &  output$corr_sun == 0 &   output$corr_temp == 0, "group2", NA))))))))))))))

# back to NA
output$corr_precip <- ifelse(output$corr_precip == 0, NA, output$corr_precip)
output$corr_sun <- ifelse(output$corr_sun == 0, NA, output$corr_sun)
output$corr_temp <- ifelse(output$corr_temp == 0, NA, output$corr_temp)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# write.table(output, "data/Pierlot_timeseries_crosscorrelations.csv",
#             quote = FALSE,
#             col.names = TRUE,
#             row.names = FALSE,
#             sep = ",")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
output <- merge(overview, output, by = c("species_full","phenophase"), all.x = TRUE)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# manuscript figure 3
# only species that have events of dormancy or turnover
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# dormancy
df_plot <- output %>%
  filter(total_nr_events > 0)

df_plot$precip_phase <- ifelse(df_plot$total_nr_events <5, "few-events",
                               ifelse(df_plot$corr_precip_timing == "0" & df_plot$corr_precip < 0, "no-lag-neg",
                                      ifelse(df_plot$corr_precip_timing == "0" & df_plot$corr_precip > 0, "no-lag-pos",
                                             ifelse(df_plot$corr_precip_timing %in% c("-1","-2","-3") & df_plot$corr_precip < 0, "lag-neg",
                                                    ifelse(df_plot$corr_precip_timing %in% c("-1","-2","-3") & df_plot$corr_precip > 0, "lag-pos",
                                                           NA)))))
df_plot$precip_phase <- ifelse(is.na(df_plot$precip_phase) & df_plot$total_nr_events >= 5, "h-no-corr", df_plot$precip_phase)


df_plot$insol_phase <- ifelse(df_plot$total_nr_events <5, "few-events",
                              ifelse(df_plot$corr_sun_timing == "0" & df_plot$corr_sun < 0, "no-lag-neg",
                                     ifelse(df_plot$corr_sun_timing == "0" & df_plot$corr_sun > 0, "no-lag-pos",
                                            ifelse(df_plot$corr_sun_timing %in% c("-1","-2","-3") & df_plot$corr_sun < 0, "lag-neg",
                                                   ifelse(df_plot$corr_sun_timing %in% c("-1","-2","-3") & df_plot$corr_sun > 0, "lag-pos",
                                                          NA)))))
df_plot$insol_phase <- ifelse(is.na(df_plot$insol_phase) & df_plot$total_nr_events >= 5, "h-no-corr", df_plot$insol_phase)

df_plot$tmax_phase <- ifelse(df_plot$total_nr_events <5, "few-events",
                             ifelse(df_plot$corr_temp_timing == "0" & df_plot$corr_temp < 0, "no-lag-neg",
                                    ifelse(df_plot$corr_temp_timing == "0" & df_plot$corr_temp > 0, "no-lag-pos",
                                           ifelse(df_plot$corr_temp_timing %in% c("-1","-2","-3") & df_plot$corr_temp < 0, "lag-neg",
                                                  ifelse(df_plot$corr_temp_timing %in% c("-1","-2","-3") & df_plot$corr_temp > 0, "lag-pos",
                                                         NA)))))
df_plot$tmax_phase <- ifelse(is.na(df_plot$tmax_phase) & df_plot$total_nr_events >= 5, "h-no-corr", df_plot$tmax_phase)


#-----------------------------
# evergreen - dormancy
#-----------------------------
df_ever_dorm <- df_plot %>%
  filter(phenophase == "leaf_dormancy",
         grepl("evergreen",deciduousness))

# for each climate variable
# get % of species with cross-correlations in different phases and signs
counts_ever_dorm <- length(df_ever_dorm$species_full)
ever_dorm_prec <- df_ever_dorm %>%
  group_by(precip_phase) %>%
  dplyr::summarise(value = length(corr_precip)/counts_ever_dorm*100) %>%
  dplyr::rename(relation = precip_phase) %>%
  add_column(variable = "precipitation")
ever_dorm_sun <- df_ever_dorm %>%
  group_by(insol_phase) %>%
  dplyr::summarise(value = length(corr_sun)/counts_ever_dorm*100) %>%
  dplyr::rename(relation = insol_phase) %>%
  add_column(variable = "sun hours")
ever_dorm_temp <- df_ever_dorm %>%
  group_by(tmax_phase) %>%
  dplyr::summarise(value = length(corr_temp)/counts_ever_dorm*100) %>%
  dplyr::rename(relation = tmax_phase) %>%
  add_column(variable = "tmax")
# rbind all
ever_dorm_summary <- rbind(ever_dorm_prec, ever_dorm_sun, ever_dorm_temp)
# % of species with 'few-events' and 'no-corr' will be added
# 1/2 at both the negative and positive correlation ends
# small hack to divide by 2 with seperate labels for relation
ever_dorm_summary$value <- ifelse(ever_dorm_summary$relation %in% c("few-events","h-no-corr"),
                                  ever_dorm_summary$value / 2, ever_dorm_summary$value)
hack <- ever_dorm_summary %>%
  filter(relation %in% "few-events") %>%
  mutate(relation = "few-events-neg")
hack2 <- ever_dorm_summary %>%
  filter(relation %in% "h-no-corr") %>%
  mutate(relation = "h-no-corr-neg")
ever_dorm_summary <- rbind(ever_dorm_summary, hack, hack2)
# ever_dorm_summary <- ever_dorm_summary %>%
#   add_column(class = "1evergeen-dormancy")
ever_dorm_summary <- ever_dorm_summary %>%
  add_column(class = "1evergeen",
             pheno = "senescence")
#-----------------------------
# evergreen - turnover
#-----------------------------
df_ever_turn <- df_plot %>%
  filter(phenophase == "leaf_turnover",
         grepl("evergreen",deciduousness))

# for each climate variable
# get % of species with cross-correlations in different phases and signs
counts_ever_turn <- length(df_ever_turn$species_full)
ever_turn_prec <- df_ever_turn %>%
  group_by(precip_phase) %>%
  dplyr::summarise(value = length(corr_precip)/counts_ever_turn*100) %>%
  dplyr::rename(relation = precip_phase) %>%
  add_column(variable = "precipitation")
ever_turn_sun <- df_ever_turn %>%
  group_by(insol_phase) %>%
  dplyr::summarise(value = length(corr_sun)/counts_ever_turn*100) %>%
  dplyr::rename(relation = insol_phase) %>%
  add_column(variable = "sun hours")
ever_turn_temp <- df_ever_turn %>%
  group_by(tmax_phase) %>%
  dplyr::summarise(value = length(corr_temp)/counts_ever_turn*100) %>%
  dplyr::rename(relation = tmax_phase) %>%
  add_column(variable = "tmax")
# rbind all
ever_turn_summary <- rbind(ever_turn_prec, ever_turn_sun, ever_turn_temp)
# % of species with 'few-events' and 'no-corr' will be added
# 1/2 at both the negative and positive correlation ends
# small hack to divide by 2 with seperate labels for relation
ever_turn_summary$value <- ifelse(ever_turn_summary$relation %in% c("few-events","h-no-corr"),
                                  ever_turn_summary$value / 2, ever_turn_summary$value)
hack <- ever_turn_summary %>%
  filter(relation %in% "few-events") %>%
  mutate(relation = "few-events-neg")
hack2 <- ever_turn_summary %>%
  filter(relation %in% "h-no-corr") %>%
  mutate(relation = "h-no-corr-neg")
ever_turn_summary <- rbind(ever_turn_summary, hack, hack2)
# ever_turn_summary <- ever_turn_summary %>%
#   add_column(class = "1evergeen-turnover")
ever_turn_summary <- ever_turn_summary %>%
  add_column(class = "1evergeen",
             pheno = "turnover")


#-----------------------------
# deciduous - dormancy
#-----------------------------
df_dec_dorm <- df_plot %>%
  filter(phenophase == "leaf_dormancy",
         grepl("deciduous",deciduousness))

# for each climate variable
# get % of species with cross-correlations in different phases and signs
counts_dec_dorm <- length(df_dec_dorm$species_full)
dec_dorm_prec <- df_dec_dorm %>%
  group_by(precip_phase) %>%
  dplyr::summarise(value = length(corr_precip)/counts_dec_dorm*100) %>%
  dplyr::rename(relation = precip_phase) %>%
  add_column(variable = "precipitation")
dec_dorm_sun <- df_dec_dorm %>%
  group_by(insol_phase) %>%
  dplyr::summarise(value = length(corr_sun)/counts_dec_dorm*100) %>%
  dplyr::rename(relation = insol_phase) %>%
  add_column(variable = "sun hours")
dec_dorm_temp <- df_dec_dorm %>%
  group_by(tmax_phase) %>%
  dplyr::summarise(value = length(corr_temp)/counts_dec_dorm*100) %>%
  dplyr::rename(relation = tmax_phase) %>%
  add_column(variable = "tmax")
# rbind all
dec_dorm_summary <- rbind(dec_dorm_prec, dec_dorm_sun, dec_dorm_temp)
# % of species with 'few-events' and 'no-corr' will be added
# 1/2 at both the negative and positive correlation ends
# small hack to divide by 2 with seperate labels for relation
dec_dorm_summary$value <- ifelse(dec_dorm_summary$relation %in% c("few-events","h-no-corr"),
                                 dec_dorm_summary$value / 2, dec_dorm_summary$value)
hack <- dec_dorm_summary %>%
  filter(relation %in% "few-events") %>%
  mutate(relation = "few-events-neg")
hack2 <- dec_dorm_summary %>%
  filter(relation %in% "h-no-corr") %>%
  mutate(relation = "h-no-corr-neg")
dec_dorm_summary <- rbind(dec_dorm_summary, hack, hack2)
# dec_dorm_summary <- dec_dorm_summary %>%
#   add_column(class = "2deciduous-dormancy")
dec_dorm_summary <- dec_dorm_summary %>%
  add_column(class = "deciduous",
             pheno = "senescence")


#-----------------------------
# deciduous - turnover
#-----------------------------
df_dec_turn <- df_plot %>%
  filter(phenophase == "leaf_turnover",
         grepl("deciduous",deciduousness))

# for each climate variable
# get % of species with cross-correlations in different phases and signs
counts_dec_turn <- length(df_dec_turn$species_full)
dec_turn_prec <- df_dec_turn %>%
  group_by(precip_phase) %>%
  dplyr::summarise(value = length(corr_precip)/counts_dec_turn*100) %>%
  dplyr::rename(relation = precip_phase) %>%
  add_column(variable = "precipitation")
dec_turn_sun <- df_dec_turn %>%
  group_by(insol_phase) %>%
  dplyr::summarise(value = length(corr_sun)/counts_dec_turn*100) %>%
  dplyr::rename(relation = insol_phase) %>%
  add_column(variable = "sun hours")
dec_turn_temp <- df_dec_turn %>%
  group_by(tmax_phase) %>%
  dplyr::summarise(value = length(corr_temp)/counts_dec_turn*100) %>%
  dplyr::rename(relation = tmax_phase) %>%
  add_column(variable = "tmax")
# rbind all
dec_turn_summary <- rbind(dec_turn_prec, dec_turn_sun, dec_turn_temp)
# % of species with 'few-events' and 'no-corr' will be added
# 1/2 at both the negative and positive correlation ends
# small hack to divide by 2 with seperate labels for relation
dec_turn_summary$value <- ifelse(dec_turn_summary$relation %in% c("few-events","h-no-corr"),
                                 dec_turn_summary$value / 2, dec_turn_summary$value)
hack <- dec_turn_summary %>%
  filter(relation %in% "few-events") %>%
  mutate(relation = "few-events-neg")
hack2 <- dec_turn_summary %>%
  filter(relation %in% "h-no-corr") %>%
  mutate(relation = "h-no-corr-neg")
dec_turn_summary <- rbind(dec_turn_summary, hack, hack2)
# dec_turn_summary <- dec_turn_summary %>%
#   add_column(class = "2deciduous-turnover")
dec_turn_summary <- dec_turn_summary %>%
  add_column(class = "deciduous",
             pheno = "turnover")


corr_summary <- rbind(ever_dorm_summary, ever_turn_summary, dec_dorm_summary, dec_turn_summary)



##################################################################
new_labels <- c("1evergeen" = "", "deciduous" = "")

p_all <- ggplot(corr_summary,
                aes(x = variable,
                    y = ifelse(relation %in% c("no-lag-neg","lag-neg","few-events-neg","h-no-corr-neg"), -value, value),
                    fill = relation)) +
  geom_col() +
  scale_y_continuous(limits = c(-73,68),
                     breaks = c(-50,-25,0,25,50),
                     sec.axis = dup_axis(name = "",
                                         breaks = c(-25,25),
                                         labels = c("negative correlations", "positive correlations"))) +
  scale_x_discrete(limits = rev(levels(as.factor(corr_summary$variable)))) +
  coord_flip() +
  scale_fill_manual(values = c("no-lag-pos" = "#018571",
                               "no-lag-neg" = "#018571",
                               "lag-pos" = "#80cdc1",
                               "lag-neg" = "#80cdc1",
                               "few-events" = "grey90",
                               "few-events-neg" = "grey90",
                               "h-no-corr" = "#dfc27d",
                               "h-no-corr-neg" = "#dfc27d"),
                    breaks = c("no-lag-neg","lag-neg","h-no-corr-neg","few-events-neg"),
                    labels = c(" in-phase   "," lag   "," no correlation   "," too few events   ")) +
  # scale_fill_manual(values = c("grey90", "grey90","#dfc27d", "#dfc27d", "#80cdc1", "#80cdc1", "#018571","#018571"),
  #                   breaks = c("no-lag-neg","lag-neg","h-no-corr-neg","few-events-neg"),
  #                   labels = c(" in-phase   "," lag   "," no correlation   "," too few events   ")
  #                   ) +
  geom_hline(yintercept = 0, color =c("grey60"), size = 0.5) +
  labs(y = "% species with neg. or pos. correlations",
       x = "")  +
  theme(text = element_text(family = "Lato", color = "#22211d"),
        panel.grid.major.x = element_line(colour = "grey60"),#, size = 0.3),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(size = 11),
        strip.background = element_rect(fill="white"),
        strip.placement = "outside",
        # axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(),
        axis.title.y = element_text(),
        axis.title.x = element_text(),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm")
  ) +
  facet_wrap( ~ class + pheno,
              ncol = 1,
              labeller = labeller(class = new_labels),
              strip.position = "left")
# p_all

ann_text <- data.frame(x = 3,
                       y = 53,
                       relation = "few-events-neg",
                       class = c("1evergeen","1evergeen","deciduous","deciduous"),
                       pheno = c("senescence", "turnover","senescence", "turnover"),
                       notes = c(paste("n = ",counts_ever_dorm),
                                 paste("n = ",counts_ever_turn),
                                 paste("n = ",counts_dec_dorm),
                                 paste("n = ",counts_dec_turn)))
p1 <- p_all +
  geom_text(data = ann_text,
            aes(x, y, label = notes),
            size = 3,
            hjust = 0,
            family = "Lato", color = "#22211d")
p1
p2 <- grid.arrange(arrangeGrob(p1,
                               right = textGrob("deciduous                              evergreen",
                                                gp=gpar(family = "Lato", color = "#22211d",fontsize=11),
                                                rot = 90,
                                                hjust = 0.4,
                                                vjust = -57)))


ggsave("manuscript/leaf_phenology/figures/Pierlot_fig3.png", p2,
       device = "png", width = 6.5, height = 5.3)

# # stationarity and autocorr precip
# adf.test(climate$precip) # p value lower then 0.05 = stationary
# kpss.test(climate$precip) # p value higher then 0.05 = stationary
# autocorr <- acf(climate$precip)


