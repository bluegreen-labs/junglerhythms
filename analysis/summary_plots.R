#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
# graphics.off()
#----- load required packages --------------------------------------------------#
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(stats)
#-------------------------------------------------------------------------------#

#-----------------------------------------------------------------------
# read in summary data (file made in summary_statistics_elizabeth.r)
#-----------------------------------------------------------------------
overview <- read.csv("data/species_meta_data.csv",
                    header = TRUE,
                    sep = ",",
                    stringsAsFactors = FALSE)

# first remove if too few events to calculate index (minimum 5 events within a species)
overview$sd_intrasp_onset_leaf_dormancy_weeks <- ifelse(overview$site_years_with_leaf_dormancy < 5, NA, overview$sd_intrasp_onset_leaf_dormancy_weeks)
overview$sd_intrasp_onset_leaf_turnover_weeks <- ifelse(overview$site_years_with_leaf_turnover < 5, NA, overview$sd_intrasp_onset_leaf_turnover_weeks)
# first remove if too few events to calculate index (minimum 2 events per individual)
overview$mean_synchrony_individuals_onset_leaf_dormancy_weeks <- ifelse(overview$mean_nr_events_within_individuals_leaf_dormancy < 2, NA, overview$mean_synchrony_individuals_onset_leaf_dormancy_weeks)
overview$mean_synchrony_individuals_onset_leaf_turnover_weeks <- ifelse(overview$mean_synchrony_individuals_onset_leaf_turnover_weeks < 2, NA, overview$mean_synchrony_individuals_onset_leaf_turnover_weeks)


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



#-----------------------------------------------------------------------
# data deciduousness
#-----------------------------------------------------------------------
dfdec <- overview %>%
  filter(deciduousness == "deciduous" | deciduousness == "evergreen")
# dormancy
dfdec_ld = dfdec[,(names(dfdec) %in% c("deciduousness",
                                       "ratio_site_years_with_leaf_dormancy",
                                       "mean_duration_leaf_dormancy_weeks",
                                       "sd_intrasp_onset_leaf_dormancy_weeks",
                                       "mean_synchrony_individuals_onset_leaf_dormancy_weeks",
                                       "mean_distance_onset_leaf_dormancy_weeks"))]
dfdec_ld$phenophase <- "leaf_dormancy"
dfdec_ld <- dfdec_ld %>%
  rename("ratio_site_years" = ratio_site_years_with_leaf_dormancy,
         "mean_duration" = mean_duration_leaf_dormancy_weeks,
         "sd_intrasp_onset" = sd_intrasp_onset_leaf_dormancy_weeks,
         "mean_synchrony_individuals_onset" = mean_synchrony_individuals_onset_leaf_dormancy_weeks,
         "mean_distance_onset" = mean_distance_onset_leaf_dormancy_weeks)

# turnover
dfdec_lt = dfdec[,(names(dfdec) %in% c("deciduousness",
                                       "ratio_site_years_with_leaf_turnover",
                                       "mean_duration_leaf_turnover_weeks",
                                       "sd_intrasp_onset_leaf_turnover_weeks",
                                       "mean_synchrony_individuals_onset_leaf_turnover_weeks",
                                       "mean_distance_onset_leaf_turnover_weeks"))]
dfdec_lt$phenophase <- "leaf_turnover"
dfdec_lt <- dfdec_lt %>%
  rename("ratio_site_years" = ratio_site_years_with_leaf_turnover,
         "mean_duration" = mean_duration_leaf_turnover_weeks,
         "sd_intrasp_onset" = sd_intrasp_onset_leaf_turnover_weeks,
         "mean_synchrony_individuals_onset" = mean_synchrony_individuals_onset_leaf_turnover_weeks,
         "mean_distance_onset" = mean_distance_onset_leaf_turnover_weeks)

# rbind dormancy and turnover
dfdec_all <- rbind(dfdec_ld,dfdec_lt)


#-----------------------------------------------------------------------
# plots deciduousness
#-----------------------------------------------------------------------
p1_dec <- ggplot(data = dfdec_all) +
  geom_boxplot(aes(x = deciduousness,
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

p2_dec <- ggplot(data = dfdec_all) +
  geom_boxplot(aes(x = deciduousness,
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

p3_dec <- ggplot(data = dfdec_all) +
  geom_boxplot(aes(x = deciduousness,
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

p4_dec <- ggplot(data = dfdec_all) +
  geom_boxplot(aes(x = deciduousness,
                   y = mean_synchrony_individuals_onset,
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

# pdf("~/Desktop/summary_plots.pdf",2.5,8)
# plot(p_all)
# dev.off()

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

# p5 <- ggplot(data = dfdec_all) +
#   geom_boxplot(aes(x = deciduousness,
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
# boxplot(dfdec$ratio_site_years_with_leaf_dormancy ~ dfdec$deciduousness,ylab = "ratio years with canopy dormancy")
# boxplot(dfdec$ratio_site_years_with_leaf_turnover ~ dfdec$deciduousness,ylab = "ratio years with canopy turnover")
# boxplot(dfecol$ratio_site_years_with_leaf_dormancy ~ dfecol$ecology,ylab = "ratio years with canopy dormancy")
# boxplot(dfecol$ratio_site_years_with_leaf_turnover ~ dfecol$ecology,ylab = "ratio years with canopy turnover")
#
# # duration event
# boxplot(dfdec$mean_duration_leaf_dormancy_weeks ~ dfdec$deciduousness,ylab = "mean duration canopy dormancy (weeks)")
# boxplot(dfdec$mean_duration_leaf_turnover_weeks ~ dfdec$deciduousness,ylab = "mean duration canopy turnover (weeks)")
# boxplot(dfecol$mean_duration_leaf_dormancy_weeks ~ dfecol$ecology,ylab = "mean duration canopy dormancy (weeks)")
# boxplot(dfecol$mean_duration_leaf_turnover_weeks ~ dfecol$ecology,ylab = "mean duration canopy turnover (weeks)")
#
# # intraspecific synchrony index
# boxplot(dfdec$sd_intrasp_onset_leaf_dormancy_weeks ~ dfdec$deciduousness,ylab = "intraspecific synchrony index leaf dormancy (weeks)")
# boxplot(dfdec$sd_intrasp_onset_leaf_turnover_weeks ~ dfdec$deciduousness,ylab = "intraspecific synchrony index leaf turnover (weeks)")
# boxplot(dfecol$sd_intrasp_onset_leaf_dormancy_weeks ~ dfecol$ecology,ylab = "intraspecific synchrony index leaf dormancy (weeks)")
# boxplot(dfecol$sd_intrasp_onset_leaf_turnover_weeks ~ dfecol$ecology,ylab = "intraspecific synchrony index leaf turnover (weeks)")
#
# # intra-individual synchrony index
# boxplot(dfdec$mean_synchrony_individuals_onset_leaf_dormancy_weeks ~ dfdec$deciduousness,ylab = "intra-individual synchrony index leaf dormancy (weeks)")
# boxplot(dfdec$mean_synchrony_individuals_onset_leaf_turnover_weeks ~ dfdec$deciduousness,ylab = "intra-individual synchrony index leaf turnover (weeks)")
# boxplot(dfecol$mean_synchrony_individuals_onset_leaf_dormancy_weeks ~ dfecol$ecology,ylab = "intra-individual synchrony index leaf dormancy (weeks)")
# boxplot(dfecol$mean_synchrony_individuals_onset_leaf_turnover_weeks ~ dfecol$ecology,ylab = "intra-individual synchrony index leaf turnover (weeks)")

