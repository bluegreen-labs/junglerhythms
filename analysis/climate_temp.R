#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
# graphics.off()
#----- load required packages --------------------------------------------------#
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(grid)
library(tseries)

#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------

# read climate data
climate <- read.csv("data/yangambi_km5_monthly_kasongo.csv")
climate <- climate %>%
  filter(year >= 1931,
         year < 1959)
climate$date <- paste(climate$year,climate$month,"15",sep = "-")
climate$date <- as.Date(climate$date, "%Y-%m-%d")
climate$yr_month <- format(as.Date(climate$date), "%Y-%m")

# # stationarity and autocorr precip
# adf.test(climate$precip) # p value lower then 0.05 = stationary
# kpss.test(climate$precip) # p value higher then 0.05 = stationary
# autocorr <- acf(climate$precip)
#
# climate <- climate %>%
#   filter(year >= 1951,
#          year < 2012)
# # stationarity and autocorr tmax
# adf.test(climate$tmax) # p value lower then 0.05 = stationary
# kpss.test(climate$tmax) # p value higher then 0.05 = stationary
# # Check for autocorrelation in a timeseries
# autocorr <- acf(climate$tmax)


climate.avg <- read.csv("data/ClimData_monthly_avg.csv",header = TRUE,sep = ",")
climate.avg = climate.avg[,(names(climate.avg) %in% c("Month","tmax_JR"))]
colnames(climate.avg)[1] <- "month"
climate <- merge(climate, climate.avg, by = "month", all.x = TRUE)



# p_climate <- ggplot() +
#   geom_col(data = climate, aes(x = date,
#                                y = tmax)) +
#   labs(y = "precip. (mm/month)",
#        x = "Year") +
#   scale_x_date(date_breaks = "1 year",
#                date_labels = "%Y",
#                limits = as.Date(c('1951-01-01','1960-01-01'))) +
#   theme_minimal() +
#   theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
#         panel.grid.minor.x =  element_blank(),
#         panel.grid.minor.y = element_blank(),
#         panel.background = element_blank(),
#         plot.background = element_rect(fill = 'white', colour = 'white'),
#         strip.text = element_text(hjust = 0),
#         axis.line.x = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1),
#         legend.position="none"
#   )



#----------------------------------------------------------------------
#--------   Phenology data - species correction Meise   ---------------
#----------------------------------------------------------------------
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
# df <- df[which(df$value != 0),]
df$join_id <- paste0("R",df$image,"-",df$image_row)

metadata <- read.csv("data/phenology_archives_species_long_format_20190626.csv",
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
#----------------------------------------------------------------------

# #----------------------------------------------------------------------
# #-------- read in census data  ----------------------------------------
# #-------- get species-specific basal area at plot and site level ------
# #----------------------------------------------------------------------
overview <- read.csv("data/species_meta_data.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)
species_list <- overview$species_full


#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------- DORMANCY x tmax hours --------------------------------------------
#----------------------------------------------------------------------------------------

output_dormancy <- data.frame(matrix(0,nrow=length(species_list),ncol=17))
colnames(output_dormancy)[1] <- "species_full"
colnames(output_dormancy)[2:14] <- c("t-6","t-5","t-4","t-3","t-2","t-1","t","t+1","t+2","t+3","t+4","t+5","t+6")
colnames(output_dormancy)[15] <- "ci"

ci <- 0.95

for(i in 1:length(species_list)){
  data_sp <- data %>%
    filter(species_full %in% species_list[i], #"Staudtia kamerunensis",
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
  # data_sp <- data_sp %>%
  #   filter(year >= 1951)

  # cross-correlation + confidence interval
  if(max(data_sp$mean_value)>0){
    corr.tmax_JR <- ccf(data_sp$tmax_JR, data_sp$mean_value,lag = 6, pl = FALSE)
    # plot(corr.tmax_JR, main = paste(species_list[i], " - dormancy x tmax_JR"))
    ci_value <- qnorm((1 + ci)/2)/sqrt(corr.tmax_JR$n.used)
  } else {
    corr.tmax_JR$acf <- NA
    ci_value <- NA
  }

  output_dormancy[i,1] <- species_list[i]
  output_dormancy[i,2:14] <- as.numeric(corr.tmax_JR$acf)
  output_dormancy[i,15] <- as.numeric(ci_value)
}

# direct correlation at time t
colnames(output_dormancy)[16] <- "corr"
colnames(output_dormancy)[17] <- "corr.timing"
output_dormancy$corr <- ifelse(abs(output_dormancy[,8]) > output_dormancy$ci, output_dormancy[,8], # 0
                               ifelse(abs(output_dormancy[,7]) > output_dormancy$ci, output_dormancy[,7], # t-1
                                      ifelse(abs(output_dormancy[,6]) > output_dormancy$ci, output_dormancy[,6], # t-2
                                             ifelse(abs(output_dormancy[,5]) > output_dormancy$ci, output_dormancy[,5], # t-3
                                                    ifelse(abs(output_dormancy[,4]) > output_dormancy$ci, output_dormancy[,4], # t-4
                                                           ifelse(abs(output_dormancy[,3]) > output_dormancy$ci, output_dormancy[,3], # t-5
                                                                                                     NA))))))
output_dormancy$corr.timing <- ifelse(abs(output_dormancy[,8]) > output_dormancy$ci, "0", # 0
                                      ifelse(abs(output_dormancy[,7]) > output_dormancy$ci, "-1", # t-1
                                             ifelse(abs(output_dormancy[,6]) > output_dormancy$ci, "-2", # t-2
                                                    ifelse(abs(output_dormancy[,5]) > output_dormancy$ci, "-3", # t-3
                                                           ifelse(abs(output_dormancy[,4]) > output_dormancy$ci, "-4", # t-4
                                                                  ifelse(abs(output_dormancy[,3]) > output_dormancy$ci, "-5", # t-5
                                                                                                            NA))))))


dormancy_tmax_JR <- output_dormancy[,(names(output_dormancy) %in% c("species_full",
                                                                     "corr",
                                                                     "corr.timing"))]
colnames(dormancy_tmax_JR) <- c("species_full",
                                 "corr_dormancy_tmax_JR",
                                 "corr_dormancy_tmax_JR_timing")


#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------- turnover x tmax hours --------------------------------------------
#----------------------------------------------------------------------------------------

output_turnover <- data.frame(matrix(0,nrow=length(species_list),ncol=17))
colnames(output_turnover)[1] <- "species_full"
colnames(output_turnover)[2:14] <- c("t-6","t-5","t-4","t-3","t-2","t-1","t","t+1","t+2","t+3","t+4","t+5","t+6")
colnames(output_turnover)[15] <- "ci"

ci <- 0.95

for(i in 1:length(species_list)){
  data_sp <- data %>%
    filter(species_full %in% species_list[i], #"Staudtia kamerunensis",
           phenophase == "leaf_turnover")
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
  # data_sp <- data_sp %>%
  #   filter(year >= 1951)

  # cross-correlation + confidence interval
  if(max(data_sp$mean_value)>0){
    corr.tmax_JR <- ccf(data_sp$tmax_JR, data_sp$mean_value,lag = 6, pl = FALSE)
    # plot(corr.tmax_JR, main = paste(species_list[i], " - turnover x tmax_JR"))
    ci_value <- qnorm((1 + ci)/2)/sqrt(corr.tmax_JR$n.used)
  } else {
    corr.tmax_JR$acf <- NA
    ci_value <- NA
  }

  output_turnover[i,1] <- species_list[i]
  output_turnover[i,2:14] <- as.numeric(corr.tmax_JR$acf)
  output_turnover[i,15] <- as.numeric(ci_value)
}

# direct correlation at time t
colnames(output_turnover)[16] <- "corr"
colnames(output_turnover)[17] <- "corr.timing"
output_turnover$corr <- ifelse(abs(output_turnover[,8]) > output_turnover$ci, output_turnover[,8], # 0
                               ifelse(abs(output_turnover[,7]) > output_turnover$ci, output_turnover[,7], # t-1
                                      ifelse(abs(output_turnover[,6]) > output_turnover$ci, output_turnover[,6], # t-2
                                             ifelse(abs(output_turnover[,5]) > output_turnover$ci, output_turnover[,5], # t-3
                                                    ifelse(abs(output_turnover[,4]) > output_turnover$ci, output_turnover[,4], # t-4
                                                           ifelse(abs(output_turnover[,3]) > output_turnover$ci, output_turnover[,3], # t-5
                                                                                                     NA))))))
output_turnover$corr.timing <- ifelse(abs(output_turnover[,8]) > output_turnover$ci, "0", # 0
                                      ifelse(abs(output_turnover[,7]) > output_turnover$ci, "-1", # t-1
                                             ifelse(abs(output_turnover[,6]) > output_turnover$ci, "-2", # t-2
                                                    ifelse(abs(output_turnover[,5]) > output_turnover$ci, "-3", # t-3
                                                           ifelse(abs(output_turnover[,4]) > output_turnover$ci, "-4", # t-4
                                                                  ifelse(abs(output_turnover[,3]) > output_turnover$ci, "-5", # t-5
                                                                                                            NA))))))


turnover_tmax_JR <- output_turnover[,(names(output_turnover) %in% c("species_full",
                                                                     "corr",
                                                                     "corr.timing"))]
colnames(turnover_tmax_JR) <- c("species_full",
                                 "corr_turnover_tmax_JR",
                                 "corr_turnover_tmax_JR_timing")

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------- Merge  and write table, will be loaded in summary_table.R ------------------
#----------------------------------------------------------------------------------------

output <- merge(overview[1], dormancy_tmax_JR, by = "species_full", all.x = TRUE)
output <- merge(output, turnover_tmax_JR, by = "species_full", all.x = TRUE)

write.table(output, "data/timeseries_tmax_correlations.csv",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# get an idea of how many species show correlations, at what average timing etc...
# in groups of deciduous and evergreen
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------


overview <- read.csv("data/species_meta_data.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)
# clear
overview$deciduousness <- ifelse(overview$species_full %in% "Pericopsis elata", "deciduous*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Trichilia welwitschii", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Copaifera mildbraedii", "deciduous*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Tridesmostemon omphalocarpoides", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Omphalocarpum lecomteanum", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Fernandoa adolfi-friderici", "deciduous*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Tabernaemontana crassa", "evergreen*",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Trichilia tessmannii", "deciduous*",overview$deciduousness)
# not so sure, limited data
overview$deciduousness <- ifelse(overview$species_full %in% "Trichilia gilletii", "evergreen* (?)",overview$deciduousness)
# not so sure, unclear phenological data
overview$deciduousness <- ifelse(overview$species_full %in% "Radlkofera calodendron", "evergreen* (?)",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Gilletiodendron mildbraedii", "evergreen* (?)",overview$deciduousness)

## two stars, in literature found as evergreen or (sometimes) deciduous
## selected a class based on the actual data
overview$deciduousness <- ifelse(overview$species_full %in% "Celtis mildbraedii", "evergreen**",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Combretum lokele", "deciduous**",overview$deciduousness)

overview$deciduousness <- ifelse(overview$species_full %in% "Homalium africanum", "evergreen**",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Quassia silvestris", "evergreen**",overview$deciduousness)

# not so sure, unclear phenological data
overview$deciduousness <- ifelse(overview$species_full %in% "Homalium longistylum", "evergreen** (?)",overview$deciduousness)
overview$deciduousness <- ifelse(overview$species_full %in% "Irvingia gabonensis", "deciduous** (?)",overview$deciduousness)



# only keep the columns you want to work with
overview = overview[,(names(overview) %in% c("species_full",
                                             "deciduousness",
                                             "site_years_with_leaf_dormancy",
                                             "site_years_with_leaf_turnover"
))]

df <- merge(overview, output, by = "species_full", all.x = TRUE)

#-------------
df_dec_dorm <- df %>%
  filter(grepl("deciduous",deciduousness)) %>%
  filter(site_years_with_leaf_dormancy > 0)

df_dec_dorm$phase <- ifelse(df_dec_dorm$corr_dormancy_tmax_JR_timing == "0", "in-phase",
                            ifelse(df_dec_dorm$corr_dormancy_tmax_JR_timing %in% c("-1","-2","-3","-4","-5"), "lag",
                                   # ifelse(df_dec_dorm$corr_dormancy_tmax_JR_timing %in% c("1","2","3"), "lead",
                                   NA))

df_dec_dorm$relation <- ifelse(df_dec_dorm$corr_dormancy_tmax_JR > 0, "pos",
                               ifelse(df_dec_dorm$corr_dormancy_tmax_JR < 0, "neg", NA))
# tapply(df_dec_dorm$corr_dormancy_tmax_JR, list(df_dec_dorm$phase,df_dec_dorm$relation), length)
counts_dec_dorm <- length(df_dec_dorm$species_full)

df_dec_dorm <- df_dec_dorm[!(is.na(df_dec_dorm$phase)),]



p_dec_dorm_tmax <- ggplot(data = df_dec_dorm,
                     aes(x = phase,
                         group = relation,
                         fill = relation)) +
  geom_histogram(aes(y=..count../counts_dec_dorm),#sum(..count..)),
                 stat = "count",
                 position = position_dodge2(preserve = "single")) +
  scale_y_continuous(limits = c(0,0.46),
                     breaks = seq(0,0.46,0.1)) +
  labs(y = "% corr with tmax",
       x = "",
       title = "") + #"Dormancy") +
  scale_fill_manual(values = c("grey70", "grey40")) +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 3),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        plot.margin = unit(c(0,0,1,1),"cm")
  )

#-------------
df_dec_turn <- df %>%
  filter(grepl("deciduous",deciduousness)) %>%
  filter(site_years_with_leaf_turnover > 0)

df_dec_turn$phase <- ifelse(df_dec_turn$corr_turnover_tmax_JR_timing == "0", "in-phase",
                            ifelse(df_dec_turn$corr_turnover_tmax_JR_timing %in% c("-1","-2","-3","-4","-5"), "lag",
                                   # ifelse(df_dec_turn$corr_turnover_tmax_JR_timing %in% c("1","2","3"), "lead",
                                   NA))

df_dec_turn$relation <- ifelse(df_dec_turn$corr_turnover_tmax_JR > 0, "pos",
                               ifelse(df_dec_turn$corr_turnover_tmax_JR < 0, "neg", NA))
# tapply(df_dec_turn$corr_dormancy_tmax_JR, list(df_dec_turn$phase,df_dec_turn$relation), length)
counts_dec_turn <- length(df_dec_turn$species_full)
df_dec_turn <- df_dec_turn[!(is.na(df_dec_turn$phase)),]

p_dec_turn_tmax <- ggplot(data = df_dec_turn,
                     aes(x = phase,
                         group = relation,
                         fill = relation,
                         na.rm = TRUE)) +
  geom_histogram(aes(y=..count../counts_dec_turn),#sum(..count..)),
                 stat = "count",
                 position = position_dodge2(preserve = "single")) +
  scale_y_continuous(limits = c(0,0.46),
                     breaks = seq(0,0.46,0.1)) +
  labs(y = "",
       x = "",
       title = "") + #"Turnover") +
  scale_fill_manual(values = c("grey70", "grey40")) +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        plot.margin = unit(c(0,0,1,0),"cm")
  )

#-------------
df_ever_dorm <- df %>%
  filter(grepl("evergreen",deciduousness)) %>%
  filter(site_years_with_leaf_dormancy > 0)

df_ever_dorm$phase <- ifelse(df_ever_dorm$corr_dormancy_tmax_JR_timing == "0", "in-phase",
                             ifelse(df_ever_dorm$corr_dormancy_tmax_JR_timing %in% c("-1","-2","-3","-4","-5"), "lag",
                                    # ifelse(df_ever_dorm$corr_dormancy_tmax_JR_timing %in% c("1","2","3"), "lead",
                                    NA))

df_ever_dorm$relation <- ifelse(df_ever_dorm$corr_dormancy_tmax_JR > 0, "pos",
                                ifelse(df_ever_dorm$corr_dormancy_tmax_JR < 0, "neg", NA))
# tapply(df_ever_dorm$corr_dormancy_tmax_JR, list(df_ever_dorm$phase,df_ever_dorm$relation), length)
counts_ever_dorm <- length(df_ever_dorm$species_full)
df_ever_dorm <- df_ever_dorm[!(is.na(df_ever_dorm$phase)),]

p_ever_dorm_tmax <- ggplot(data = df_ever_dorm,
                      aes(x = phase,
                          group = relation,
                          fill = relation)) +
  geom_histogram(aes(y=..count../counts_ever_dorm),#sum(..count..)),
                 stat = "count",
                 position = position_dodge2(preserve = "single")) +
  scale_y_continuous(limits = c(0,0.46),
                     breaks = seq(0,0.46,0.1)) +
  labs(y = "",
       x = "",
       title = "") + #"Dormancy") +
  scale_fill_manual(values = c("grey70", "grey40")) +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        plot.margin = unit(c(0,0,1,0),"cm")
  )

#-------------
df_ever_turn <- df %>%
  filter(grepl("evergreen",deciduousness)) %>%
  filter(site_years_with_leaf_turnover > 0)

df_ever_turn$phase <- ifelse(df_ever_turn$corr_turnover_tmax_JR_timing == "0", "in-phase",
                             ifelse(df_ever_turn$corr_turnover_tmax_JR_timing %in% c("-1","-2","-3","-4","-5"), "lag",
                                    # ifelse(df_ever_turn$corr_turnover_tmax_JR_timing %in% c("1","2","3"), "lead",
                                    NA))

df_ever_turn$relation <- ifelse(df_ever_turn$corr_turnover_tmax_JR > 0, "pos",
                                ifelse(df_ever_turn$corr_turnover_tmax_JR < 0, "neg", NA))
# tapply(df_ever_turn$corr_dormancy_tmax_JR, list(df_ever_turn$phase,df_ever_turn$relation), length)
counts_ever_turn <- length(df_ever_turn$species_full)
df_ever_turn <- df_ever_turn[!(is.na(df_ever_turn$phase)),]

p_ever_turn_tmax <- ggplot(data = df_ever_turn,
                      aes(x = phase,
                          group = relation,
                          fill = relation)) +
  geom_histogram(aes(y=..count../counts_ever_turn),#sum(..count..)),
                 stat = "count",
                 position = position_dodge2(preserve = "single")) +
  scale_y_continuous(limits = c(0,0.46),
                     breaks = seq(0,0.46,0.1)) +
  labs(y = "",
       x = "",
       title = "") + #"Turnover") +
  scale_fill_manual(values = c("grey70", "grey40")) +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(colour = "grey89", size = 0.3),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(hjust = 0),
        axis.line.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        plot.margin = unit(c(0,0,1,0),"cm")
  )



grid.arrange(arrangeGrob(p_dec_dorm_tmax, p_dec_turn_tmax,top=textGrob("Deciduous",gp=gpar(fontsize=15)),ncol=2,widths = c(1,0.8)),
             arrangeGrob(p_ever_dorm_tmax, p_ever_turn_tmax,top=textGrob("Evergreen",gp=gpar(fontsize=15)), ncol=2),
             ncol = 2,widths = c(1,0.8))



#-------------
df_dec_dorm <- df %>%
  filter(grepl("deciduous",deciduousness)) %>%
  filter(site_years_with_leaf_dormancy > 0)

apply(df_dec_dorm, 2, function(x){sum(!is.na(x))})

test <- df_dec_dorm %>%
  filter(corr_dormancy_tmax_JR_timing == "0")

test <- df_dec_dorm %>%
  filter(corr_dormancy_tmax_JR_timing %in% c("-1","-2","-3","-4","-5"))
mean(as.numeric(test$corr_dormancy_tmax_JR_timing), na.rm=TRUE)

test <- df_dec_dorm %>%
  filter(corr_dormancy_tmax_JR_timing %in% c("1","2","3","4","5"))
mean(as.numeric(test$corr_dormancy_tmax_JR_timing), na.rm=TRUE)


#-------------
df_dec_turn <- df %>%
  filter(grepl("deciduous",deciduousness)) %>%
  filter(site_years_with_leaf_turnover > 0)

apply(df_dec_turn, 2, function(x){sum(!is.na(x))})

test <- df_dec_turn %>%
  filter(corr_turnover_tmax_JR_timing == "0")

test <- df_dec_turn %>%
  filter(corr_turnover_tmax_JR_timing %in% c("-1","-2","-3","-4","-5"))
mean(as.numeric(test$corr_turnover_tmax_JR_timing), na.rm=TRUE)

test <- df_dec_turn %>%
  filter(corr_turnover_tmax_JR_timing %in% c("1","2","3","4","5"))
mean(as.numeric(test$corr_turnover_tmax_JR_timing), na.rm=TRUE)

#-------------
df_ever_dorm <- df %>%
  filter(grepl("evergreen",deciduousness)) %>%
  filter(site_years_with_leaf_dormancy > 0)

apply(df_ever_dorm, 2, function(x){sum(!is.na(x))})

test <- df_ever_dorm %>%
  filter(corr_dormancy_tmax_JR_timing == "0")

test <- df_ever_dorm %>%
  filter(corr_dormancy_tmax_JR_timing %in% c("-1","-2","-3","-4","-5"))
mean(as.numeric(test$corr_dormancy_tmax_JR_timing), na.rm=TRUE)

test <- df_ever_dorm %>%
  filter(corr_dormancy_tmax_JR_timing %in% c("1","2","3","4","5"))
mean(as.numeric(test$corr_dormancy_tmax_JR_timing), na.rm=TRUE)

#-------------
df_ever_turn <- df %>%
  filter(grepl("evergreen",deciduousness)) %>%
  filter(site_years_with_leaf_turnover > 0)

apply(df_ever_turn, 2, function(x){sum(!is.na(x))})

test <- df_ever_turn %>%
  filter(corr_turnover_tmax_JR_timing == "0")

test <- df_ever_turn %>%
  filter(corr_turnover_tmax_JR_timing %in% c("-1","-2","-3","-4","-5"))
mean(as.numeric(test$corr_turnover_tmax_JR_timing), na.rm=TRUE)

test <- df_ever_turn %>%
  filter(corr_turnover_tmax_JR_timing %in% c("1","2","3","4","5"))
mean(as.numeric(test$corr_turnover_tmax_JR_timing), na.rm=TRUE)
