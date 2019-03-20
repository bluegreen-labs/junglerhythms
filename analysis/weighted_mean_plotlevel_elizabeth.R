# join data for weighted mean
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(stats)
library(Hmisc)
#----------------------------------------------------------------------
#--------   Phenology data - species correction Meise   ---------------
#----------------------------------------------------------------------
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
# df <- df[which(df$value != 0),]
df$join_id <- paste0("R",df$image,"-",df$image_row)

metadata <- read.csv("data/phenology_archives_species_long_format_20190319.csv",
                     header = TRUE, sep = ",")
metadata$join_id <- paste(metadata$image,metadata$row, sep = "-")

# test merge the two tables based upon the unique join ID
data <- merge(df, metadata, by = c("join_id"), all.x = TRUE)
data$species_full <- paste(data$genus_Meise, data$species_Meise)

# remove column id.x and rename id.y to id (--> in id.y, empty ids are renamed to EK1, EK2, etc...)
data = data[,!(names(data) %in% "id.x")]
data <- data %>%
  rename("id" = id.y)
data$id <- as.character(data$id)
data$species_full <- as.character(data$species_full)
#----------------------------------------------------------------------


# read in census data and convert basal area
census <- read.csv2("data/YGB_ForestPlotsNET.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)
census$basal_area = pi*(census$C1DBH4/2000)^2
# calculate sum basal area across all species across all mixed plots
censusMix <- census %>%
  filter(grepl("MIX",Plot)) %>%
  group_by(Plot, Species) %>%
  dplyr::summarise(BA = sum(basal_area))

censusMix <- na.omit(censusMix, cols="basal_area")
censusMix <- censusMix %>% rename("species_full" = Species)


traits <- read.table("data/Dataset_traits_African_trees.csv",
                        header = TRUE,
                        sep = ",")
traits$layer <- ifelse(traits$Dmax_lit_cm > 39, "canopy", "understory")



data <- inner_join(data, censusMix, by = "species_full")
data <- inner_join(data, traits, by = "species_full")

climate <- read.table("~/Dropbox/Phenology_JR/Manuscript/Phenology_Leaf_draft/ClimData_test.csv",
                      header = TRUE,
                      sep = ",")


#-----------------------------------------------------------------------
#------------ Leaf turnover  -------------------------------------------
#-----------------------------------------------------------------------
data_LT <- data %>%
  filter(phenophase == "leaf_turnover") %>%
  # filter(layer == "canopy") %>%
  # filter(Ecology == "shade") %>%
  na.omit()

# data_LT$value <- 1

data_LT <- data_LT %>%
  group_by(Plot,species_full, week) %>%
  summarise(mean_week = mean(value, na.rm = TRUE))

data_LT <- inner_join(data_LT, censusMix, by = c("Plot","species_full"))

### reference website for confidence intervals
### https://stats.stackexchange.com/questions/224029/calculating-weighted-cis-and-interpretation
final_LT <- data_LT %>%
  group_by(Plot,week) %>%
  summarise(ss = sum(mean_week*BA, na.rm = TRUE)/(sum(BA))*100,
            ss2 = wtd.mean(mean_week, BA,na.rm = TRUE)*100,
            weighted_sd = sqrt(wtd.var(mean_week, BA,na.rm = TRUE))*100,
            weighted_CI = ss2 + 1.96*weighted_sd)

final_LT_all <- data_LT %>%
  group_by(week) %>%
  summarise(ss = sum(mean_week*BA, na.rm = TRUE)/(sum(BA))*100,
            ss2 = wtd.mean(mean_week, BA,na.rm = TRUE)*100,
            weighted_sd = sqrt(wtd.var(mean_week, BA,na.rm = TRUE))*100,
            weighted_CI = ss2 + 1.96*weighted_sd)
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#------------ Leaf turnover  -------------------------------------------
#-----------------------------------------------------------------------
data_LD <- data %>%
  filter(phenophase == "leaf_dormancy") %>%
  # filter(layer == "canopy") %>%
  # filter(Ecology == "shade") %>%
  na.omit()

# data_LD$value <- 1

data_LD <- data_LD %>%
  group_by(Plot,species_full, week) %>%
  summarise(mean_week = mean(value, na.rm = TRUE))

data_LD <- inner_join(data_LD, censusMix, by = c("Plot","species_full"))

### reference website for confidence intervals
### https://stats.stackexchange.com/questions/224029/calculating-weighted-cis-and-interpretation
final_LD <- data_LD %>%
  group_by(Plot,week) %>%
  summarise(ss = sum(mean_week*BA, na.rm = TRUE)/(sum(BA))*100,
            ss2 = wtd.mean(mean_week, BA,na.rm = TRUE)*100,
            weighted_sd = sqrt(wtd.var(mean_week, BA,na.rm = TRUE))*100,
            weighted_CI = ss2 + 1.96*weighted_sd)

final_LD_all <- data_LD %>%
  group_by(week) %>%
  summarise(ss = sum(mean_week*BA, na.rm = TRUE)/(sum(BA))*100,
            ss2 = wtd.mean(mean_week, BA,na.rm = TRUE)*100,
            weighted_sd = sqrt(wtd.var(mean_week, BA,na.rm = TRUE))*100,
            weighted_CI = ss2 + 1.96*weighted_sd)
#-----------------------------------------------------------------------

# data_peak1 <- data_w %>%
#   filter(week == c(1:8,41:48)) %>%
#   filter(mean_week > 0)
# final_peak1 <- data_peak1 %>%
#   group_by(species_full) %>%
#   summarise(ss = wtd.mean(mean_week, BA,na.rm = TRUE)*100)

# sp_driven <- data_LT %>%
#   filter(mean_week > 0)
# sp_driven$wtd <- sp_driven$mean_week * sp_driven$BA
# sp_driven$month <- ifelse(sp_driven$week %in% c(1:4),"Jan",
#                           ifelse(sp_driven$week %in% c(5:8),"Feb",
#                           ifelse(sp_driven$week %in% c(9:12),"Mar",
#                           ifelse(sp_driven$week %in% c(11:16),"Apr",
#                           ifelse(sp_driven$week %in% c(17:20),"May",
#                           ifelse(sp_driven$week %in% c(21:24),"Jun",
#                           ifelse(sp_driven$week %in% c(25:28),"Jul",
#                           ifelse(sp_driven$week %in% c(29:32),"Aug",
#                           ifelse(sp_driven$week %in% c(33:36),"Sep",
#                           ifelse(sp_driven$week %in% c(37:40),"Oct",
#                           ifelse(sp_driven$week %in% c(41:44),"Nov",
#                           ifelse(sp_driven$week %in% c(45:48),"Dec",
#                           NA))))))))))))
# # write to file
# write.table(sp_driven, "data/species_contribution_standlevel_turnover.csv",
#             quote = FALSE,
#             col.names = TRUE,
#             row.names = FALSE,
#             sep = ",")

####------------------------------------------------------------------------
####-----  This is Koens original code    ----------------------------------
####------------------------------------------------------------------------
# data_w <- data_w %>%
#   group_by(species_full, week) %>%
#   summarise(mean_week = mean(value, na.rm = TRUE) * unique(BAperc))
# final <- data_w %>%
#   group_by(week) %>%
#   summarise(ss = sum(mean_week, na.rm = TRUE))
####------------------------------------------------------------------------
par(mfrow=c(2,1))
plot(final_LT_all$ss,type="l",lwd=2,lty=1,ylim=c(0,2.6),col="red",ylab="frequency event",xlab="week (48-week year)",main="Leaf turnover")
# plot(1, type="n", xlab="", ylab="", xlim=c(0, 48), ylim=c(0, 5))
plots <- c('MIX2', 'MIX3', 'MIX4','MIX5','MIX6')
for (i in 1:length(plots)) {
  d <- subset(final_LT, Plot==plots[i])
  lines(x=d$week, y=d$ss,col="grey47")
}
rect(0,-0.5,8,15,col = rgb(0.5,0.5,0.5,1/4),border=NA) # Jan-Febr
rect(21,-0.5,28,15,col = rgb(0.5,0.5,0.5,1/4),border=NA) # Jun-Jul
rect(45,-0.5,48,15,col = rgb(0.5,0.5,0.5,1/4),border=NA) # Dec
# lines(final_LT$weighted_CI,type="l",lty=2)
points(final_LT_all$ss,pch=16,col="red")


b <- barplot(climate$prec,ylim=c(0,400),ylab="Precip & PAR",
             names.arg=c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
# lines(x=b,climate$insol, col="red")
lines(x=b,climate$PAR_Hauser, col="red",lty=2,lwd=2)


p_dormancy <- ggplot(final_LD) +
  # geom_line(aes(week, ss, shape = Plot), col="grey") +
  geom_point(aes(week, ss, shape = Plot), col="grey40") +
  # geom_smooth(aes(week, ss), span = 0.2, se = FALSE, col = "red", size = 1.2) +
  geom_line(data = final_LD_all, aes(week, ss), col="red", size=1.2) +
  geom_point(data = final_LD_all, aes(week, ss), col="red", size=2) +
  scale_x_continuous(limits = c(1,49),
                     breaks = seq(1,48,4),
                     labels = month.abb) +
  annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 2.1, alpha = .2) + # jan - febr
  annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 2.1, alpha = .2) + # jun - jul
  annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 2.1, alpha = .2) + # dec
  labs(y = "freq. canopy dormancy",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        # axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none",
        plot.margin=unit(c(0,0,0,0.2),"cm")
  )


p_turnover <- ggplot(final_LT) +
  # geom_line(aes(week, ss, shape = Plot), col="grey") +
  geom_point(aes(week, ss, shape = Plot), col="grey40") +
  # geom_smooth(aes(week, ss), span = 0.2, se = FALSE, col = "red", size = 1.2) +
  geom_line(data = final_LT_all, aes(week, ss), col="red", size=1.2) +
  geom_point(data = final_LT_all, aes(week, ss), col="red", size=2) +
  scale_x_continuous(limits = c(1,49),
                     breaks = seq(1,48,4),
                     labels = month.abb) +
  annotate("rect", xmin = 1, xmax = 9, ymin = 0, ymax = 2.5, alpha = .2) + # jan - febr
  annotate("rect", xmin = 21, xmax = 29, ymin = 0, ymax = 2.5, alpha = .2) + # jun - jul
  annotate("rect", xmin = 45, xmax = 49, ymin = 0, ymax = 2.5, alpha = .2) + # dec
  labs(y = "freq. canopy turnover",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        # axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none",
        plot.margin=unit(c(0,0,0.2,0.2),"cm")
  )


p_precip <- ggplot(climate) +
  geom_col(aes(x = Month,
               y = prec),
           col = "grey70",
           fill = "grey70") +
  scale_x_continuous(limits = c(0.5,12.5),
                     breaks = seq(1,12,1),
                     labels = month.abb) +
  labs(y = "precip. (mm)",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="none",
        plot.margin=unit(c(0,0,0,0.2),"cm")
  )

p_par <- ggplot(climate) +
  geom_line(aes(x = Month,
                y = PAR_Hauser),
            size = 1.2) +
  scale_x_continuous(limits = c(0.5,12.5),
                     breaks = seq(1,12,1),
                     labels = month.abb) +
  labs(y = "PAR",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        # axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none",
        plot.margin=unit(c(0,0,0,0.2),"cm")
  )
## for secondary y-axis
# geom_line(aes(x = Month,
#               y = PAR_Hauser/1.5),
#           size = 1.2) +
#   scale_y_continuous(sec.axis = sec_axis(~.*1.5, name = "PAR")) +


grid.arrange(p_modis, p_dormancy, p_turnover, p_par,p_precip, heights = c(3,3,3.3,1,2))
# grid.arrange(p_dormancy,p_turnover,p_climate, heights = c(1,1,1))
