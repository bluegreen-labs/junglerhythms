#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
# graphics.off()
#----- load required packages --------------------------------------------------#
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(stats)
library(agricolae)
#-------------------------------------------------------------------------------#

#-----------------------------------------------------------------------
# read in summary data (file made in summary_statistics_elizabeth.r)
#-----------------------------------------------------------------------
overview <- read.csv("data/species_meta_data_phase2.csv",
                    header = TRUE,
                    sep = ",",
                    stringsAsFactors = FALSE)

# # # first remove if too few events to calculate index (minimum 5 events within a species)
# overview$sd_intrasp_onset_leaf_dormancy_weeks <- ifelse(overview$site_years_with_leaf_dormancy < 5, NA, overview$sd_intrasp_onset_leaf_dormancy_weeks)
# overview$sd_intrasp_onset_leaf_turnover_weeks <- ifelse(overview$site_years_with_leaf_turnover < 5, NA, overview$sd_intrasp_onset_leaf_turnover_weeks)

# # first remove if too few events to calculate index (minimum 2 events per individual)
overview$mean_synchrony_individuals_onset_leaf_dormancy_weeks <- ifelse(overview$mean_nr_events_within_individuals_leaf_dormancy < 2, NA, overview$mean_synchrony_individuals_onset_leaf_dormancy_weeks)
overview$mean_synchrony_individuals_onset_leaf_turnover_weeks <- ifelse(overview$mean_synchrony_individuals_onset_leaf_turnover_weeks < 2, NA, overview$mean_synchrony_individuals_onset_leaf_turnover_weeks)



# to sort easily for making the eventual tables -. uniform label deciduous or evergreen
overview$dec_label <- ifelse(overview$deciduousness %in% c("deciduous","deciduous*","deciduous**","deciduous* (?)","deciduous** (?)"), "dec",
                             ifelse(overview$deciduousness %in% c("evergreen","evergreen*","evergreen**","evergreen* (?)","evergreen** (?)"), "ever", "NA"))



#---------------------------
# synchrony at the individual level
overview_ind <- read.csv("data/synchrony_individuals_phase2.csv",
                         header = TRUE,
                         sep = ",",
                         stringsAsFactors = FALSE)
# first remove if too few events to calculate index (minimum 3 events per individual)
overview_ind$onset_sd_weeks <- ifelse(overview_ind$nr_events_onset < 3, NA, overview_ind$onset_sd_weeks)




#-----------------------------------------------------------------------
# data deciduousness
#-----------------------------------------------------------------------
# dormancy
df_ld <- overview %>%
  dplyr::select(species_full,
                dec_label,
                ratio_site_years_with_leaf_dormancy,
                mean_duration_leaf_dormancy_weeks,
                sd_intrasp_onset_leaf_dormancy_weeks,
                median_intrasp_onset_leaf_dormancy_weeks)
df_ld$phenophase <- "leaf_dormancy"
df_ld <- df_ld %>%
  dplyr::rename("ratio_site_years" = ratio_site_years_with_leaf_dormancy,
                "mean_duration" = mean_duration_leaf_dormancy_weeks,
                "sd_intrasp_onset" = sd_intrasp_onset_leaf_dormancy_weeks,
                "timing" = median_intrasp_onset_leaf_dormancy_weeks)

# turnover
df_lt <- overview %>%
  dplyr::select(species_full,
                dec_label,
                ratio_site_years_with_leaf_turnover,
                mean_duration_leaf_turnover_weeks,
                sd_intrasp_onset_leaf_turnover_weeks,
                median_intrasp_onset_leaf_turnover_weeks)
df_lt$phenophase <- "leaf_turnover"
df_lt <- df_lt %>%
  dplyr::rename("ratio_site_years" = ratio_site_years_with_leaf_turnover,
                "mean_duration" = mean_duration_leaf_turnover_weeks,
                "sd_intrasp_onset" = sd_intrasp_onset_leaf_turnover_weeks,
                "timing" = median_intrasp_onset_leaf_turnover_weeks)

# rbind dormancy and turnover
df_all <- rbind(df_ld,df_lt)

df_all$timing2 <- df_all$timing/52*48
df_all$season <- ifelse(df_all$timing2 <=8, "LDS",
                        ifelse(df_all$timing2 > 8 & df_all$timing2 <= 20, "SWS",
                               ifelse(df_all$timing2 > 20 & df_all$timing2 <=28, "SDS",
                                      ifelse(df_all$timing2 > 28 & df_all$timing2 <=44, "LWS",
                                             ifelse(df_all$timing2 > 44, "LDS",NA)))))

df_all$LDSrest <- ifelse(df_all$timing2 <=8, "LDS",
                         ifelse(df_all$timing2 > 8 & df_all$timing2 <= 20, "rest",
                                ifelse(df_all$timing2 > 20 & df_all$timing2 <=28, "rest",
                                       ifelse(df_all$timing2 > 28 & df_all$timing2 <=44, "rest",
                                              ifelse(df_all$timing2 > 44, "LDS",NA)))))
df_all$halfyearrest <- ifelse(df_all$timing2 <=8, "halfyear",
                         ifelse(df_all$timing2 > 8 & df_all$timing2 <= 20, "rest",
                                ifelse(df_all$timing2 > 20 & df_all$timing2 <=28, "rest",
                                       ifelse(df_all$timing2 > 28 & df_all$timing2 <=44, "halfyear",
                                              ifelse(df_all$timing2 > 44, "halfyear",NA)))))
# #-----------------------------------------------------------------------
# # some statistics
# #-----------------------------------------------------------------------
# all deciduous compared to all evergreen
a <- aov(df_all$ratio_site_years ~ df_all$dec_label)
summary(a)
a <- aov(df_all$mean_duration ~ df_all$dec_label)
summary(a)
a <- aov(df_all$mean_duration ~ df_all$phenophase)
summary(a)
test <- df_all %>%
  group_by(phenophase) %>%
  dplyr::summarise(mean_value = mean(mean_duration, na.rm=TRUE),
                   sd_value = sd(mean_duration, na.rm=TRUE))
test


a <- aov(df_all$sd_intrasp_onset ~ df_all$dec_label)
summary(a)
boxplot(df_all$sd_intrasp_onset ~ df_all$dec_label)

tx <- with(df_all, interaction(dec_label, LDSrest))
a <- aov(sd_intrasp_onset ~ tx, data = df_all)
summary(a)
posthoc <- HSD.test(a, "tx", group=TRUE)
posthoc

test <- df_all %>%
  group_by(dec_label, LDSrest) %>%
  dplyr::summarise(mean_value = mean(sd_intrasp_onset, na.rm=TRUE),
                   sd_value = sd(sd_intrasp_onset, na.rm=TRUE))
test

# within deciduous, campare dormancy v turnover
only_dec <- df_all %>%
  filter(dec_label == "ever")
a <- aov(only_dec$ratio_site_years ~ only_dec$phenophase)
summary(a)
a <- aov(only_dec$mean_duration ~ only_dec$phenophase)
summary(a)
a <- aov(only_dec$sd_intrasp_onset ~ only_dec$phenophase)
summary(a)

tx <- with(only_dec, interaction(phenophase, LDSrest))
a <- aov(sd_intrasp_onset ~ tx, data = only_dec)
summary(a)
posthoc <- HSD.test(a, "tx", group=TRUE)
posthoc




test <- only_dec %>%
  filter(!is.na(sd_intrasp_onset)) %>%
  group_by(LDSrest) %>%
  dplyr::summarise(mean_value = mean(sd_intrasp_onset),
                   sd_value = sd(sd_intrasp_onset),
                   nr = length(sd_intrasp_onset))
test

mean(only_dec$sd_intrasp_onset, na.rm=TRUE)
sd(only_dec$sd_intrasp_onset, na.rm=TRUE)








#-----------------------------------------------------------------------
# # merge intraspecific with intra-annual dissynchrony to see if their is a statistical difference
# test <- df_all[,(names(df_all) %in% c("species_full",
#                                             "dec_label",
#                                             "sd_intrasp_onset",
#                                             "phenophase"))]
# test <- test %>%
#   rename("onset_sd_weeks" = sd_intrasp_onset)
# test$onset <- "intraspecific"
# test2 <- overview_ind[,(names(overview_ind) %in% c("species_full",
#                                              "dec_label",
#                                              "onset_sd_weeks",
#                                              "phenophase"))]
# test2$onset <- "intraindividual"
# test3 <- rbind(test,test2)
# a <- aov(test3$onset_sd_weeks ~ test3$onset)
# summary(a)
# mean(test$onset_sd_weeks, na.rm = TRUE)
# sd(test$onset_sd_weeks, na.rm = TRUE)
# mean(test2$onset_sd_weeks, na.rm = TRUE)
# sd(test2$onset_sd_weeks, na.rm = TRUE)
#-----------------------------------------------------------------------
overview_ind <- overview_ind %>%
  group_by(species_full, phenophase) %>%
  dplyr::summarize(mean_interann = mean(onset_sd_weeks))
overview_sp <- df_all %>%
  dplyr::select(species_full,
                dec_label,
                phenophase,
                season,
                LDSrest,
                halfyearrest)
overview_ind <- merge(overview_ind, overview_sp, by = c("species_full", "phenophase"), all.x = TRUE)

a <- aov(overview_ind$mean_interann ~ overview_ind$dec_label)
summary(a)
boxplot(overview_ind$mean_interann ~ overview_ind$dec_label)

tx <- with(overview_ind, interaction(phenophase, dec_label))
a <- aov(mean_interann ~ tx, data = overview_ind)
summary(a)
posthoc <- HSD.test(a, "tx", group=TRUE)
posthoc


ind_dorm <- overview_ind %>%
  filter(phenophase == "leaf_dormancy") #leaf_dormancy, leaf_turnover
a <- aov(ind_dorm$mean_interann ~ ind_dorm$dec_label)
summary(a)
boxplot(ind_dorm$mean_interann ~ ind_dorm$dec_label)

ind_dec <- overview_ind %>%
  filter(dec_label == "ever")
a <- aov(ind_dec$mean_interann ~ ind_dec$LDSrest)
summary(a)
boxplot(ind_dec$mean_interann ~ ind_dec$LDSrest)

test <- ind_dec %>%
  group_by(LDSrest) %>%
  dplyr::summarise(mean_value = mean(mean_interann, na.rm=TRUE),
                   sd_value = sd(mean_interann, na.rm=TRUE))
test

mean(ind_dec$mean_interann, na.rm=TRUE)
sd(ind_dec$mean_interann, na.rm=TRUE)


tx <- with(ind_dec, interaction(phenophase, LDSrest))
a <- aov(mean_interann ~ tx, data = ind_dec)
summary(a)
posthoc <- HSD.test(a, "tx", group=TRUE)
posthoc
#-----------------------------------------------------------------------
# plots dec_label
#-----------------------------------------------------------------------
p1_dec <- ggplot(data = df_all) +
  geom_boxplot(aes(x = dec_label,
                   y = ratio_site_years,
                   color = phenophase),
               size = 1) +
  scale_colour_manual(values = c("black","darkgreen")) +
  labs(y = "",
       x = "") +
  scale_y_continuous(name = "% event years",
                     limits = c(0,1),
                     breaks = seq(0,1,0.25)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 4),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0.5),"cm")
  )

p2_dec <- ggplot(data = df_all) +
  geom_boxplot(aes(x = dec_label,
                   y = mean_duration,
                   color = phenophase),
               size = 1) +
  scale_colour_manual(values = c("black","darkgreen")) +
  labs(y = "event duration (w)",
       x = "") +
  scale_y_continuous(limits = c(0,15),
                     breaks = seq(0,15,5)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 4),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0.5),"cm")
  )

p3_dec <- ggplot(data = df_all) +
  geom_boxplot(aes(x = dec_label,
                   y = sd_intrasp_onset,
                   color = phenophase),
               size = 1) +
  scale_colour_manual(values = c("black","darkgreen")) +
  labs(y = "intra-specific synchr. (w)",
       x = "") +
  scale_y_continuous(limits = c(0,22),
                     breaks = seq(0,20,5)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 4),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0.5),"cm")
  )

p4_dec <-
  ggplot(data = overview_ind) +
  geom_boxplot(aes(x = dec_label,
                   y = onset_sd_weeks,
                   color = phenophase),
               size = 1) +
  scale_colour_manual(values = c("black","darkgreen")) +
  labs(y = "inter-annual synchr. (w)",
       x = "") +
  scale_y_continuous(limits = c(0,20),
                     breaks = seq(0,20,5)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 4),
        legend.position = "none",
        plot.margin = unit(c(0,0,1,0.5),"cm")
  )

#-----------------------------------------------------------------------
# combine plots
#-----------------------------------------------------------------------
# make the widths of the plots with different y-axis equal
# https://stackoverflow.com/questions/30402930/align-x-axes-of-box-plot-and-line-plot-using-ggplot
p1_dec <- ggplot_gtable(ggplot_build(p1_dec))
p2_dec <- ggplot_gtable(ggplot_build(p2_dec))
p1_dec$widths <-p2_dec$widths
p_all <- grid.arrange(p1_dec, p2_dec, p3_dec, p4_dec, heights = c(1,1,1,1.2))

# pdf("~/Desktop/figure2_summary_plots.pdf",2.5,8)
# plot(p_all)
# dev.off()



# #-----------------------------------------------------------------------
# # data ecology
# #-----------------------------------------------------------------------
# dfecol <- overview %>%
#   filter(ecology == "shade" | ecology == "sun")
# # dormancy
# dfecol_ld = dfecol[,(names(dfecol) %in% c("ecology",
#                                        "ratio_site_years_with_leaf_dormancy",
#                                        "mean_duration_leaf_dormancy_weeks",
#                                        "sd_intrasp_onset_leaf_dormancy_weeks",
#                                        "mean_synchrony_individuals_onset_leaf_dormancy_weeks",
#                                        "mean_distance_onset_leaf_dormancy_weeks"))]
# dfecol_ld$phenophase <- "leaf_dormancy"
# dfecol_ld <- dfecol_ld %>%
#   rename("ratio_site_years" = ratio_site_years_with_leaf_dormancy,
#          "mean_duration" = mean_duration_leaf_dormancy_weeks,
#          "sd_intrasp_onset" = sd_intrasp_onset_leaf_dormancy_weeks,
#          "mean_synchrony_individuals_onset" = mean_synchrony_individuals_onset_leaf_dormancy_weeks,
#          "mean_distance_onset" = mean_distance_onset_leaf_dormancy_weeks)
#
# # turnover
# dfecol_lt = dfecol[,(names(dfecol) %in% c("ecology",
#                                        "ratio_site_years_with_leaf_turnover",
#                                        "mean_duration_leaf_turnover_weeks",
#                                        "sd_intrasp_onset_leaf_turnover_weeks",
#                                        "mean_synchrony_individuals_onset_leaf_turnover_weeks",
#                                        "mean_distance_onset_leaf_turnover_weeks"))]
# dfecol_lt$phenophase <- "leaf_turnover"
# dfecol_lt <- dfecol_lt %>%
#   rename("ratio_site_years" = ratio_site_years_with_leaf_turnover,
#          "mean_duration" = mean_duration_leaf_turnover_weeks,
#          "sd_intrasp_onset" = sd_intrasp_onset_leaf_turnover_weeks,
#          "mean_synchrony_individuals_onset" = mean_synchrony_individuals_onset_leaf_turnover_weeks,
#          "mean_distance_onset" = mean_distance_onset_leaf_turnover_weeks)
#
# # rbind dormancy and turnover
# dfecol_all <- rbind(dfecol_ld,dfecol_lt)

# #-----------------------------------------------------------------------
# # plots ecology
# #-----------------------------------------------------------------------
# p1_ecol <- ggplot(data = dfecol_all) +
#   geom_boxplot(aes(x = ecology,
#                    y = ratio_site_years,
#                    color = phenophase),
#                size = 1) +
#   scale_colour_manual(values = c("black","darkgreen")) +
#   labs(y = "",
#        x = "") +
#   scale_y_continuous(limits = c(0,1),
#                      breaks = seq(0,1,0.5)) +
#   theme_minimal() +
#   theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
#         panel.grid.minor.x =  element_blank(),
#         panel.grid.minor.y = element_blank(),
#         panel.background = element_blank(),
#         plot.background = element_rect(fill = 'white', colour = 'white'),
#         strip.text = element_text(hjust = 0),
#         axis.line.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.title.y = element_blank(),
#         legend.position = "none",
#         plot.margin = unit(c(0,1,0,0),"cm")
#   )
#
# p2_ecol <- ggplot(data = dfecol_all) +
#   geom_boxplot(aes(x = ecology,
#                    y = mean_duration,
#                    color = phenophase),
#                size = 1) +
#   scale_colour_manual(values = c("black","darkgreen")) +
#   labs(y = "",
#        x = "") +
#   scale_y_continuous(limits = c(0,15),
#                      breaks = seq(0,15,5)) +
#   theme_minimal() +
#   theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
#         panel.grid.minor.x =  element_blank(),
#         panel.grid.minor.y = element_blank(),
#         panel.background = element_blank(),
#         plot.background = element_rect(fill = 'white', colour = 'white'),
#         strip.text = element_text(hjust = 0),
#         axis.line.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.title.y = element_blank(),
#         legend.position = "none",
#         plot.margin = unit(c(0,1,0,0),"cm")
#   )
#
# p3_ecol <- ggplot(data = dfecol_all) +
#   geom_boxplot(aes(x = ecology,
#                    y = sd_intrasp_onset,
#                    color = phenophase),
#                size = 1) +
#   scale_colour_manual(values = c("black","darkgreen")) +
#   labs(y = "",
#        x = "") +
#   scale_y_continuous(limits = c(0,22),
#                      breaks = seq(0,20,5)) +
#   theme_minimal() +
#   theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
#         panel.grid.minor.x =  element_blank(),
#         panel.grid.minor.y = element_blank(),
#         panel.background = element_blank(),
#         plot.background = element_rect(fill = 'white', colour = 'white'),
#         strip.text = element_text(hjust = 0),
#         axis.line.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.title.y = element_blank(),
#         legend.position = "none",
#         plot.margin = unit(c(0,1,0,0),"cm")
#   )
#
# p4_ecol <- ggplot(data = dfecol_all) +
#   geom_boxplot(aes(x = ecology,
#                    y = mean_synchrony_individuals_onset,
#                    color = phenophase),
#                size = 1) +
#   geom_vline(xintercept = 0) +
#
#
#   scale_colour_manual(values = c("black","darkgreen")) +
#   labs(y = "",
#        x = "") +
#   scale_y_continuous(limits = c(0,20),
#                      breaks = seq(0,20,5)) +
#   theme_minimal() +
#   theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
#         panel.grid.minor.x =  element_blank(),
#         panel.grid.minor.y = element_blank(),
#         panel.background = element_blank(),
#         plot.background = element_rect(fill = 'white', colour = 'white'),
#         strip.text = element_text(hjust = 0),
#         axis.line.x = element_blank(),
#         axis.text.x = element_text(size = 12, colour = "black"),
#         axis.title.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.title.y = element_blank(),
#         legend.position = "none",
#         plot.margin = unit(c(0,1,1,0),"cm")
#   )
#
#-----------------------------------------------------------------------
# combine plots
#-----------------------------------------------------------------------
# # make the widths of the plots with different y-axis equal
# # https://stackoverflow.com/questions/30402930/align-x-axes-of-box-plot-and-line-plot-using-ggplot
# p1_dec <- ggplot_gtable(ggplot_build(p1_dec))
# p2_dec <- ggplot_gtable(ggplot_build(p2_dec))
# p1_dec$widths <-p2_dec$widths

# p1_ecol <- ggplot_gtable(ggplot_build(p1_ecol))
# p2_ecol <- ggplot_gtable(ggplot_build(p2_ecol))
# p1_ecol$widths <-p2_ecol$widths

# grid.arrange(p1_dec,p1_ecol, p2_dec, p2_ecol, p3_dec, p3_ecol, p4_dec, p4_ecol, ncol = 2, heights = c(1,1,1,1.2))





#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

# p5 <- ggplot(data = df_all) +
#   geom_boxplot(aes(x = dec_label,
#                    y = mean_distance_onset,
#                    color = phenophase),
#                size = 1) +
#   scale_colour_manual(values = c("black","darkgreen")) +
#   labs(y = "mean distance of onset phenophase to all species in the community (weeks)",
#        x = "") +
#   theme_minimal() +
#   theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
#         panel.grid.minor.x =  element_blank(),
#         panel.grid.minor.y = element_blank(),
#         panel.background = element_blank(),
#         plot.background = element_rect(fill = 'white', colour = 'white'),
#         strip.text = element_text(hjust = 0),
#         axis.line.x = element_blank(),
#         axis.text.x = element_text(size = 12, colour = "black"),
#         axis.title.x = element_blank(),
#         legend.position = "none",
#         plot.margin = unit(c(0,0,1,1),"cm")
#   )




# par(mfrow=c(4,4),
#     mar=c(6,5,2,2))
# # ratio eventyears:siteyears
# boxplot(df$ratio_site_years_with_leaf_dormancy ~ df$dec_label,ylab = "ratio years with canopy dormancy")
# boxplot(df$ratio_site_years_with_leaf_turnover ~ df$dec_label,ylab = "ratio years with canopy turnover")
# boxplot(dfecol$ratio_site_years_with_leaf_dormancy ~ dfecol$ecology,ylab = "ratio years with canopy dormancy")
# boxplot(dfecol$ratio_site_years_with_leaf_turnover ~ dfecol$ecology,ylab = "ratio years with canopy turnover")
#
# # duration event
# boxplot(df$mean_duration_leaf_dormancy_weeks ~ df$dec_label,ylab = "mean duration canopy dormancy (weeks)")
# boxplot(df$mean_duration_leaf_turnover_weeks ~ df$dec_label,ylab = "mean duration canopy turnover (weeks)")
# boxplot(dfecol$mean_duration_leaf_dormancy_weeks ~ dfecol$ecology,ylab = "mean duration canopy dormancy (weeks)")
# boxplot(dfecol$mean_duration_leaf_turnover_weeks ~ dfecol$ecology,ylab = "mean duration canopy turnover (weeks)")
#
# # intraspecific synchrony index
# boxplot(df$sd_intrasp_onset_leaf_dormancy_weeks ~ df$dec_label,ylab = "intraspecific synchrony index leaf dormancy (weeks)")
# boxplot(df$sd_intrasp_onset_leaf_turnover_weeks ~ df$dec_label,ylab = "intraspecific synchrony index leaf turnover (weeks)")
# boxplot(dfecol$sd_intrasp_onset_leaf_dormancy_weeks ~ dfecol$ecology,ylab = "intraspecific synchrony index leaf dormancy (weeks)")
# boxplot(dfecol$sd_intrasp_onset_leaf_turnover_weeks ~ dfecol$ecology,ylab = "intraspecific synchrony index leaf turnover (weeks)")
#
# # intra-individual synchrony index
# boxplot(df$mean_synchrony_individuals_onset_leaf_dormancy_weeks ~ df$dec_label,ylab = "intra-individual synchrony index leaf dormancy (weeks)")
# boxplot(df$mean_synchrony_individuals_onset_leaf_turnover_weeks ~ df$dec_label,ylab = "intra-individual synchrony index leaf turnover (weeks)")
# boxplot(dfecol$mean_synchrony_individuals_onset_leaf_dormancy_weeks ~ dfecol$ecology,ylab = "intra-individual synchrony index leaf dormancy (weeks)")
# boxplot(dfecol$mean_synchrony_individuals_onset_leaf_turnover_weeks ~ dfecol$ecology,ylab = "intra-individual synchrony index leaf turnover (weeks)")

