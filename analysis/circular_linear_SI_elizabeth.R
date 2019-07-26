#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
# graphics.off()
#----- load required packages --------------------------------------------------#
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(gridExtra)
#----- source required functions -----------------------------------------------#
# source("analysis/circular_plot_leaf_turnover_elizabeth.r")
# source("analysis/linear_plot_senescence_elizabeth.R")
source("analysis/circular_linear_function_elizabeth.R")
#-------------------------------------------------------------------------------#

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

overview <- read.csv("data/SI_table2.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)

overview_dec <- overview %>%
  filter(grepl("deciduous",deciduousness)) %>%
  filter(percentage_site_years_with_leaf_dormancy > 0 | percentage_site_years_with_leaf_turnover > 0) %>%
  arrange_at("species_full")

sp_dec <- overview_dec$species_full
sp_dec1 <- paste(sp_dec[1:5],collapse = "|")
sp_dec2 <- paste(sp_dec[6:10],collapse = "|")
sp_dec3 <- paste(sp_dec[11:15],collapse = "|")
sp_dec4 <- paste(sp_dec[16:20],collapse = "|")
sp_dec5 <- paste(sp_dec[21:25],collapse = "|")
sp_dec6 <- paste(sp_dec[26:30],collapse = "|")
sp_dec7 <- paste(sp_dec[31:33],collapse = "|")



# sequence 1
plot_dec1 <- circular_linear_plot(data,
                                  species_name = sp_dec1,
                                  leg_gradient = c(0,0.3,1),
                                  title_name = "(b) Deciduous")
pdf("~/Desktop/deciduous1.pdf",8.85,12.1)
plot(plot_dec1)
dev.off()

# sequence 2
plot_dec2 <- circular_linear_plot(data,
                                  species_name = sp_dec2,
                                  leg_gradient = c(0,0.3,1),
                                  title_name = "(b) Deciduous - continued")
pdf("~/Desktop/deciduous2.pdf",8.85,12.1)
plot(plot_dec2)
dev.off()

# sequence 3
plot_dec3 <- circular_linear_plot(data,
                                  species_name = sp_dec3,
                                  leg_gradient = c(0,0.3,1),
                                  title_name = "(b) Deciduous - continued")
pdf("~/Desktop/deciduous3.pdf",8.85,12.1)
plot(plot_dec3)
dev.off()

# sequence 4
plot_dec4 <- circular_linear_plot(data,
                                  species_name = sp_dec4,
                                  leg_gradient = c(0,0.3,1),
                                  title_name = "(b) Deciduous - continued")
pdf("~/Desktop/deciduous4.pdf",8.85,12.1)
plot(plot_dec4)
dev.off()

# sequence 5
plot_dec5 <- circular_linear_plot(data,
                                  species_name = sp_dec5,
                                  leg_gradient = c(0,0.3,1),
                                  title_name = "(b) Deciduous - continued")
pdf("~/Desktop/deciduous5.pdf",8.85,12.1)
plot(plot_dec5)
dev.off()

# sequence 6
plot_dec6 <- circular_linear_plot(data,
                                  species_name = sp_dec6,
                                  leg_gradient = c(0,0.3,1),
                                  title_name = "(b) Deciduous - continued")
pdf("~/Desktop/deciduous6.pdf",8.85,12.1) # freq. with 1 digit
plot(plot_dec6)
dev.off()

# sequence 7
# changed legend position
plot_dec7 <- circular_linear_plot(data,
                                  species_name = sp_dec7,
                                  leg_pos = c(1,0.175),
                                  leg_gradient = c(0,0.3,1),
                                  title_name = "(b) Deciduous - continued")
pdf("~/Desktop/deciduous7.pdf",8.85, 7.5)
plot(plot_dec7)
dev.off()


#------------------------------------------------------------------------
overview_ever <- overview %>%
  filter(grepl("evergreen",deciduousness)) %>%
  filter(percentage_site_years_with_leaf_dormancy > 0 | percentage_site_years_with_leaf_turnover > 0) %>%
  arrange_at("species_full")

sp_ever <- overview_ever$species_full
sp_ever1 <- paste(sp_ever[1:4],collapse = "|") # up to 4 for the first page which includes figure description
sp_ever2 <- paste(sp_ever[5:9],collapse = "|")
sp_ever3 <- paste(sp_ever[10:14],collapse = "|")
sp_ever4 <- paste(sp_ever[15:19],collapse = "|")
sp_ever5 <- paste(sp_ever[20:24],collapse = "|")
sp_ever6 <- paste(sp_ever[25:29],collapse = "|")
sp_ever7 <- paste(sp_ever[30:34],collapse = "|")
sp_ever8 <- paste(sp_ever[35:39],collapse = "|")
sp_ever9 <- paste(sp_ever[40:44],collapse = "|")
sp_ever10 <- paste(sp_ever[45:49],collapse = "|")
sp_ever11 <- paste(sp_ever[50:54],collapse = "|")
sp_ever12 <- paste(sp_ever[55:59],collapse = "|")
sp_ever13 <- paste(sp_ever[60:64],collapse = "|")


# sequence 1
# changed legend position
plot_ever1 <- circular_linear_plot(data,
                                   species_name = sp_ever1,
                                   leg_pos = c(1,0.125),
                                   title_name = "(a) Evergreen")
pdf("~/Desktop/evergreen1.pdf",8.85,9.8) # freq. with 3 digits
plot(plot_ever1)
dev.off()

# sequence 2
plot_ever2 <- circular_linear_plot(data,
                                   species_name = sp_ever2,
                                   title_name = "(a) Evergreen - continued")
pdf("~/Desktop/evergreen2.pdf",8.85,12.1)
plot(plot_ever2)
dev.off()

# sequence 3
plot_ever3 <- circular_linear_plot(data,
                                   species_name = sp_ever3,
                                   title_name = "(a) Evergreen - continued")
pdf("~/Desktop/evergreen3.pdf",8.85,12.1)
plot(plot_ever3)
dev.off()

# sequence 4
plot_ever4 <- circular_linear_plot(data,
                                   species_name = sp_ever4,
                                   title_name = "(a) Evergreen - continued")
pdf("~/Desktop/evergreen4.pdf",8.85,12.1)
plot(plot_ever4)
dev.off()

# sequence 5
plot_ever5 <- circular_linear_plot(data,
                                   species_name = sp_ever5,
                                   title_name = "(a) Evergreen - continued")
pdf("~/Desktop/evergreen5.pdf",8.85,12.1)
plot(plot_ever5)
dev.off()

# sequence 6
plot_ever6 <- circular_linear_plot(data,
                                   species_name = sp_ever6,
                                   title_name = "(a) Evergreen - continued")
pdf("~/Desktop/evergreen6.pdf",8.85,12.1)
plot(plot_ever6)
dev.off()


# sequence 7
plot_ever7 <- circular_linear_plot(data,
                                   species_name = sp_ever7,
                                   title_name = "(a) Evergreen - continued")
pdf("~/Desktop/evergreen7.pdf",8.85,12.1)
plot(plot_ever7)
dev.off()

# sequence 8
plot_ever8 <- circular_linear_plot(data,
                                   species_name = sp_ever8,
                                   title_name = "(a) Evergreen - continued")
pdf("~/Desktop/evergreen8.pdf",8.85,12.1)
plot(plot_ever8)
dev.off()

# sequence 9
plot_ever9 <- circular_linear_plot(data,
                                   species_name = sp_ever9,
                                   title_name = "(a) Evergreen - continued")
pdf("~/Desktop/evergreen9.pdf",8.85,12.1)
plot(plot_ever9)
dev.off()

# sequence 10
plot_ever10 <- circular_linear_plot(data,
                                   species_name = sp_ever10,
                                   title_name = "(a) Evergreen - continued")
pdf("~/Desktop/evergreen10.pdf",8.85,12.1)
plot(plot_ever10)
dev.off()

# sequence 11
plot_ever11 <- circular_linear_plot(data,
                                    species_name = sp_ever11,
                                    title_name = "(a) Evergreen - continued")
pdf("~/Desktop/evergreen11.pdf",8.85,12.1)
plot(plot_ever11)
dev.off()

# sequence 12
plot_ever12 <- circular_linear_plot(data,
                                    species_name = sp_ever12,
                                    title_name = "(a) Evergreen - continued")
pdf("~/Desktop/evergreen12.pdf",8.85,12.1)
plot(plot_ever12)
dev.off()

# sequence 13
plot_ever13 <- circular_linear_plot(data,
                                    species_name = sp_ever13,
                                    title_name = "(a) Evergreen - continued")
pdf("~/Desktop/evergreen13.pdf",8.85,12.1)
plot(plot_ever13)
dev.off()
