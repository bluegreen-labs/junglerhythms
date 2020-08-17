

#-----------------------------------------------------------------------
#------------ Leaf flushing  -------------------------------------------
#-----------------------------------------------------------------------
data_flush <- data %>%
  filter(phenophase == "leaf_dormancy") %>%
  filter(basal_area_site > 0)
data_flush$flushing_date <- paste(data_flush$year, data_flush$week, sep = "-")

data_flush_ones <- data_flush[which(data_flush$value != 0),]
flushing_timing <- data_flush_ones %>%
  group_by(species_full, id) %>%
  do(event_length(.))

flushing_timing$flushing_date <- paste(flushing_timing$year_end, flushing_timing$week_end +1, sep = "-")

flushing_timing <- flushing_timing[,(names(flushing_timing) %in% c("species_full",
                                                                   "id",
                                                                   "flushing_date"))]
flushing_timing$flushing_value <- 1

data_flush <- merge(data_flush, flushing_timing, by = c("species_full","id","flushing_date"), all.x = TRUE)
data_flush$flushing_value <- ifelse(is.na(data_flush$flushing_value),"0", data_flush$flushing_value)
data_flush$flushing_value <- ifelse(is.na(data_flush$value),NA, data_flush$flushing_value)
data_flush$flushing_value <- as.numeric(data_flush$flushing_value)

data_flush <- data_flush %>%
  group_by(species_full, week) %>%
  summarise(mean_week = mean(flushing_value, na.rm = TRUE),
            count_week = sum(flushing_value),
            total_week = length(flushing_value))

# filter out species with too few total siteyears to have a meaningfull average year
data_flush <- data_flush %>%
  filter(total_week >= minimum_siteyears)

data_flush$mean_week <- ifelse(data_flush$mean_week < minimum_event_frequency, 0, data_flush$mean_week)

data_flush_site <- inner_join(data_flush, census_site, by = c("species_full"))
data_flush_plot <- inner_join(data_flush, census_plot, by = c("species_full"))


final_flush_site <- data_flush_site %>%
  group_by(week) %>%
  summarise(ss = sum(mean_week*basal_area_site, na.rm = TRUE)/total_basal_area_site *100) #,

final_flush_plot <- data_flush_plot %>%
  group_by(Plot, week) %>%
  summarise(ss = sum(mean_week*basal_area_plot, na.rm = TRUE))
final_flush_plot <- inner_join(final_flush_plot, total_basal_area_plot, by = c("Plot"))
final_flush_plot <- final_flush_plot %>%
  mutate(ss = ss/total_basal_area_plot*100)
#-----------------------------------------------------------------------

p_flushing <- ggplot() +
  geom_point(data = final_flush_plot,
             aes(week, ss, shape = Plot),
             col="grey40") +
  geom_smooth(data = final_flush_site,
              aes(week, ss), span = 0.2, se = FALSE, col = "black", size = 1.2) +
  geom_point(data = final_flush_site,
             aes(week, ss),
             col="black",
             size=2) +
  scale_x_continuous(limits = c(1,49),
                     breaks = seq(1,48,4),
                     labels = month.abb) +
  # scale_y_continuous(limits = c(0,2.5),
  #                    breaks = seq(0,2,0.5),
  #                    labels = scales::number_format(accuracy = 0.1)) +
  annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 1, alpha = .1) + # jan - febr
  annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 1, alpha = .1) + # jun - jul
  annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 1, alpha = .1) + # dec
  labs(y = "freq. onset canopy flushing",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        # axis.text.x = element_text(angle = 90, hjust = 1,vjust = 1.5,size = 10), # vjust to center the label
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 3),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0.5),"cm")
  )
print(p_flushing)




#-----------------------------------------------------------------------
#------------ Onset Leaf turnover  -------------------------------------------
#-----------------------------------------------------------------------
data_onset_turnover <- data %>%
  filter(phenophase == "leaf_turnover") %>%
  filter(basal_area_site > 0)
data_onset_turnover$turnover_date <- paste(data_onset_turnover$year, data_onset_turnover$week, sep = "-")

data_onset_turnover_ones <- data_onset_turnover[which(data_onset_turnover$value != 0),]
turnover_timing <- data_onset_turnover_ones %>%
  group_by(species_full, id) %>%
  do(event_length(.))

turnover_timing$turnover_date <- paste(turnover_timing$year_start, turnover_timing$week_start, sep = "-")

turnover_timing <- turnover_timing[,(names(turnover_timing) %in% c("species_full",
                                                                   "id",
                                                                   "turnover_date"))]
turnover_timing$turnover_value <- 1

data_onset_turnover <- merge(data_onset_turnover, turnover_timing, by = c("species_full","id","turnover_date"), all.x = TRUE)
data_onset_turnover$turnover_value <- ifelse(is.na(data_onset_turnover$turnover_value),"0", data_onset_turnover$turnover_value)
data_onset_turnover$turnover_value <- ifelse(is.na(data_onset_turnover$value),NA, data_onset_turnover$turnover_value)
data_onset_turnover$turnover_value <- as.numeric(data_onset_turnover$turnover_value)

data_onset_turnover <- data_onset_turnover %>%
  group_by(species_full, week) %>%
  summarise(mean_week = mean(turnover_value, na.rm = TRUE),
            count_week = sum(turnover_value),
            total_week = length(turnover_value))

# filter out species with too few total siteyears to have a meaningfull average year
data_onset_turnover <- data_onset_turnover %>%
  filter(total_week >= minimum_siteyears)

data_onset_turnover$mean_week <- ifelse(data_onset_turnover$mean_week < minimum_event_frequency, 0, data_onset_turnover$mean_week)

data_onset_turnover_site <- inner_join(data_onset_turnover, census_site, by = c("species_full"))
data_onset_turnover_plot <- inner_join(data_onset_turnover, census_plot, by = c("species_full"))


final_onset_turnover_site <- data_onset_turnover_site %>%
  group_by(week) %>%
  summarise(ss = sum(mean_week*basal_area_site, na.rm = TRUE)/total_basal_area_site *100) #,

final_onset_turnover_plot <- data_onset_turnover_plot %>%
  group_by(Plot, week) %>%
  summarise(ss = sum(mean_week*basal_area_plot, na.rm = TRUE))
final_onset_turnover_plot <- inner_join(final_onset_turnover_plot, total_basal_area_plot, by = c("Plot"))
final_onset_turnover_plot <- final_onset_turnover_plot %>%
  mutate(ss = ss/total_basal_area_plot*100)
#-----------------------------------------------------------------------

p_onset_turnover <- ggplot() +
  geom_point(data = final_onset_turnover_plot,
             aes(week, ss, shape = Plot),
             col="grey40") +
  geom_smooth(data = final_onset_turnover_site,
              aes(week, ss), span = 0.2, se = FALSE, col = "black", size = 1.2) +
  geom_point(data = final_onset_turnover_site,
             aes(week, ss),
             col="black",
             size=2) +
  scale_x_continuous(limits = c(1,49),
                     breaks = seq(1,48,4),
                     labels = month.abb) +
  # scale_y_continuous(limits = c(0,2.5),
  #                    breaks = seq(0,2,0.5),
  #                    labels = scales::number_format(accuracy = 0.1)) +
  annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 1, alpha = .1) + # jan - febr
  annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 1, alpha = .1) + # jun - jul
  annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 1, alpha = .1) + # dec
  labs(y = "freq. onset canopy turnover",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        # axis.text.x = element_text(angle = 90, hjust = 1,vjust = 1.5,size = 10), # vjust to center the label
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 3),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0.5),"cm")
  )
print(p_onset_turnover)

#-----------------------------------------------------------------------
#------------ Onset Leaf dormancy  -------------------------------------------
#-----------------------------------------------------------------------
data_onset_dormancy <- data %>%
  filter(phenophase == "leaf_dormancy") %>%
  filter(basal_area_site > 0)
data_onset_dormancy$dormancy_date <- paste(data_onset_dormancy$year, data_onset_dormancy$week, sep = "-")

data_onset_dormancy_ones <- data_onset_dormancy[which(data_onset_dormancy$value != 0),]
dormancy_timing <- data_onset_dormancy_ones %>%
  group_by(species_full, id) %>%
  do(event_length(.))

dormancy_timing$dormancy_date <- paste(dormancy_timing$year_start, dormancy_timing$week_start, sep = "-")

dormancy_timing <- dormancy_timing[,(names(dormancy_timing) %in% c("species_full",
                                                                   "id",
                                                                   "dormancy_date"))]
dormancy_timing$dormancy_value <- 1

data_onset_dormancy <- merge(data_onset_dormancy, dormancy_timing, by = c("species_full","id","dormancy_date"), all.x = TRUE)
data_onset_dormancy$dormancy_value <- ifelse(is.na(data_onset_dormancy$dormancy_value),"0", data_onset_dormancy$dormancy_value)
data_onset_dormancy$dormancy_value <- ifelse(is.na(data_onset_dormancy$value),NA, data_onset_dormancy$dormancy_value)
data_onset_dormancy$dormancy_value <- as.numeric(data_onset_dormancy$dormancy_value)

data_onset_dormancy <- data_onset_dormancy %>%
  group_by(species_full, week) %>%
  summarise(mean_week = mean(dormancy_value, na.rm = TRUE),
            count_week = sum(dormancy_value),
            total_week = length(dormancy_value))

# filter out species with too few total siteyears to have a meaningfull average year
data_onset_dormancy <- data_onset_dormancy %>%
  filter(total_week >= minimum_siteyears)

data_onset_dormancy$mean_week <- ifelse(data_onset_dormancy$mean_week < minimum_event_frequency, 0, data_onset_dormancy$mean_week)

data_onset_dormancy_site <- inner_join(data_onset_dormancy, census_site, by = c("species_full"))
data_onset_dormancy_plot <- inner_join(data_onset_dormancy, census_plot, by = c("species_full"))


final_onset_dormancy_site <- data_onset_dormancy_site %>%
  group_by(week) %>%
  summarise(ss = sum(mean_week*basal_area_site, na.rm = TRUE)/total_basal_area_site *100) #,

final_onset_dormancy_plot <- data_onset_dormancy_plot %>%
  group_by(Plot, week) %>%
  summarise(ss = sum(mean_week*basal_area_plot, na.rm = TRUE))
final_onset_dormancy_plot <- inner_join(final_onset_dormancy_plot, total_basal_area_plot, by = c("Plot"))
final_onset_dormancy_plot <- final_onset_dormancy_plot %>%
  mutate(ss = ss/total_basal_area_plot*100)
#-----------------------------------------------------------------------

p_onset_dormancy <- ggplot() +
  geom_point(data = final_onset_dormancy_plot,
             aes(week, ss, shape = Plot),
             col="grey40") +
  geom_smooth(data = final_onset_dormancy_site,
              aes(week, ss), span = 0.2, se = FALSE, col = "black", size = 1.2) +
  geom_point(data = final_onset_dormancy_site,
             aes(week, ss),
             col="black",
             size=2) +
  scale_x_continuous(limits = c(1,49),
                     breaks = seq(1,48,4),
                     labels = month.abb) +
  # scale_y_continuous(limits = c(0,2.5),
  #                    breaks = seq(0,2,0.5),
  #                    labels = scales::number_format(accuracy = 0.1)) +
  annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 1, alpha = .1) + # jan - febr
  annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 1, alpha = .1) + # jun - jul
  annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 1, alpha = .1) + # dec
  labs(y = "freq. onset canopy dormancy",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        # axis.text.x = element_text(angle = 90, hjust = 1,vjust = 1.5,size = 10), # vjust to center the label
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 3),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0.5),"cm")
  )
print(p_onset_dormancy)


#-----------------------------------------------------------------------
#------------ plots -------  -------------------------------------------
#-----------------------------------------------------------------------
colnames(final_onset_dormancy_site)[2] <- "ss_onset_dorm"
colnames(final_onset_turnover_site)[2] <- "ss_onset_turn"
colnames(final_flush_site)[2] <- "ss_flush"

combined_onset <- merge(final_onset_dormancy_site, final_onset_turnover_site, by = c("week"), all.x = TRUE)
combined_onset <- merge(combined_onset, final_flush_site, by = c("week"), all.x = TRUE)
combined_onset$ss_onset_senescence <- combined_onset$ss_onset_dorm + combined_onset$ss_onset_turn
combined_onset$ss_onset_flushing <- combined_onset$ss_flush + combined_onset$ss_onset_turn

p_combined_onset <- ggplot() +
  annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 1, alpha = .1) + # jan - febr
  annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 1, alpha = .1) + # jun - jul
  annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 1, alpha = .1) + # dec
  # # both
  # geom_smooth(data = combined_onset,
  #             aes(week, ss_onset_senescence, colour = "line4"), span = 0.2, se = FALSE, size = 1.2, linetype = "twodash",
  #             show.legend = TRUE) +
  # geom_smooth(data = combined_onset,
  #             aes(week, ss_onset_flushing, colour = "line5"), span = 0.2, se = FALSE, size = 1.2, linetype = "twodash",
  #             show.legend = TRUE) +

  # turnover
  geom_smooth(data = combined_onset,
              aes(week, ss_onset_turn, colour = "line2"), span = 0.2, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(data = combined_onset,
             aes(week, ss_onset_turn),
             col="#018571",
             shape = 1,
             stroke = 1.3) +
  # dormancy
  geom_smooth(data = combined_onset,
              aes(week, ss_onset_dorm, colour = "line1"), span = 0.2, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(data = combined_onset,
             aes(week, ss_onset_dorm),
             col="#a6611a",
             size=2) +
  # flushing
  geom_smooth(data = combined_onset,
              aes(week, ss_flush, colour = "line3"), span = 0.2, se = FALSE, size = 1.2,
              show.legend = TRUE) +
  geom_point(data = combined_onset,
             aes(week, ss_flush),
             col="purple",
             size=2) +


  scale_colour_manual(values = c("line1" = "#a6611a", "line2" = "#018571", "line3" = "purple", "line4" = "grey30", "line5" = "red"),
                      labels = c(" dormacy"," turnover"," flushing"," senescence", " total flushing")) +
  scale_x_continuous(limits = c(1,49),
                     breaks = seq(1,48,4),
                     labels = month.abb) +
  labs(y = "freq. canopy phenophase",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        # axis.text.x = element_text(angle = 90, hjust = 1,vjust = 1.5,size = 10), # vjust to center the label
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 3),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        plot.margin = unit(c(0,0,0,0.5),"cm")
  )

# p_combined_onset <- ggplot() +
#   annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 1, alpha = .1) + # jan - febr
#   annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 1, alpha = .1) + # jun - jul
#   annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 1, alpha = .1) + # dec
#   # both
#   geom_line(data = combined_onset,
#               aes(week, ss_onset_senescence, colour = "line4"), size = 1.2, linetype = "twodash",
#               show.legend = TRUE) +
#   geom_line(data = combined_onset,
#               aes(week, ss_onset_flushing, colour = "line5"), size = 1.2, linetype = "twodash",
#               show.legend = TRUE) +
#
#   # # turnover
#   # geom_line(data = combined_onset,
#   #             aes(week, ss_onset_turn, colour = "line2"), size = 1.2,
#   #             show.legend = TRUE) +
#   # geom_point(data = combined_onset,
#   #            aes(week, ss_onset_turn),
#   #            col="#018571",
#   #            shape = 1,
#   #            stroke = 1.3) +
#   # # dormancy
#   # geom_line(data = combined_onset,
#   #             aes(week, ss_onset_dorm, colour = "line1"), size = 1.2,
#   #             show.legend = TRUE) +
#   # geom_point(data = combined_onset,
#   #            aes(week, ss_onset_dorm),
#   #            col="#a6611a",
#   #            size=2) +
#   # # flushing
#   # geom_line(data = combined_onset,
#   #             aes(week, ss_flush, colour = "line3"), size = 1.2,
#   #             show.legend = TRUE) +
#   # geom_point(data = combined_onset,
#   #            aes(week, ss_flush),
#   #            col="purple",
#   #            size=2) +
#
#
#   scale_colour_manual(values = c("line1" = "#a6611a", "line2" = "#018571", "line3" = "purple", "line4" = "grey30", "line5" = "red"),
#                       labels = c(" dormacy"," turnover"," flushing"," senescence", " total flushing")) +
#   scale_x_continuous(limits = c(1,49),
#                      breaks = seq(1,48,4),
#                      labels = month.abb) +
#   labs(y = "freq. canopy phenophase",
#        x = "") +
#   theme_minimal() +
#   theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
#         panel.grid.minor.x =  element_blank(),
#         panel.grid.minor.y = element_blank(),
#         panel.background = element_blank(),
#         plot.background = element_rect(fill = 'white', colour = 'white'),
#         strip.text = element_text(hjust = 0),
#         axis.line.x = element_blank(),
#         axis.text.x = element_blank(),
#         # axis.text.x = element_text(angle = 90, hjust = 1,vjust = 1.5,size = 10), # vjust to center the label
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(vjust = 3),
#         legend.position = "top",
#         legend.title = element_blank(),
#         legend.text = element_text(size = 11),
#         plot.margin = unit(c(0,0,0,0.5),"cm")
#   )

p_combined_onset


