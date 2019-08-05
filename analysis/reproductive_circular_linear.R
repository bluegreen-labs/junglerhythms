#----- reset your R session. ---------------------------------------------------#
rm(list=ls())
# graphics.off()
#----- load required packages --------------------------------------------------#
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(gridExtra)
#----- source required functions -----------------------------------------------#
source("analysis/reproductive_circular_linear_function.R")
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


sp <- overview$species_full
sp1 <- paste(sp[1:5],collapse = "|") # paste(sp[c(1,2,4,5)],collapse = "|")
sp2 <- paste(sp[6:10],collapse = "|")
sp3 <- paste(sp[11:15],collapse = "|")
sp4 <- paste(sp[16:20],collapse = "|")
sp5 <- paste(sp[21:25],collapse = "|")
sp6 <- paste(sp[26:30],collapse = "|")
sp7 <- paste(sp[31:35],collapse = "|")
sp8 <- paste(sp[36:40],collapse = "|")
sp9 <- paste(sp[41:45],collapse = "|")
sp10 <- paste(sp[46:50],collapse = "|")
sp11 <- paste(sp[51:55],collapse = "|")
sp12 <- paste(sp[56:60],collapse = "|")
sp13 <- paste(sp[61:65],collapse = "|")
sp14 <- paste(sp[66:70],collapse = "|")
sp15 <- paste(sp[71:75],collapse = "|")
sp16 <- paste(sp[76:80],collapse = "|")
sp17 <- paste(sp[81:85],collapse = "|")
sp18 <- paste(sp[86:90],collapse = "|")
sp19 <- paste(sp[91:95],collapse = "|")
sp20 <- paste(sp[96:100],collapse = "|")
sp21 <- paste(sp[101:105],collapse = "|")
sp22 <- paste(sp[106:110],collapse = "|")
sp23 <- paste(sp[111:115],collapse = "|")
sp24 <- paste(sp[116:120],collapse = "|")
sp25 <- paste(sp[121:125],collapse = "|")
sp26 <- paste(sp[126:128],collapse = "|")

# sequence 1
plot1 <- circular_linear_plot(data,
                              species_name = sp1,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro1.pdf",8.85,12.1)
plot(plot1)
dev.off()

# sequence 2
plot2 <- circular_linear_plot(data,
                              species_name = sp2,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro2.pdf",8.85,12.1)
plot(plot2)
dev.off()

# sequence 3
plot3 <- circular_linear_plot(data,
                              species_name = sp3,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro3.pdf",8.85,12.1)
plot(plot3)
dev.off()

# sequence 4
plot4 <- circular_linear_plot(data,
                              species_name = sp4,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro4.pdf",8.85,12.1)
plot(plot4)
dev.off()

# sequence 5
plot5 <- circular_linear_plot(data,
                              species_name = sp5,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro5.pdf",8.85,12.1)
plot(plot5)
dev.off()

# sequence 6
plot6 <- circular_linear_plot(data,
                              species_name = sp6,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro6.pdf",8.85,12.1)
plot(plot6)
dev.off()

# sequence 7
plot7 <- circular_linear_plot(data,
                              species_name = sp7,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro7.pdf",8.85,12.1)
plot(plot7)
dev.off()

# sequence 8
plot8 <- circular_linear_plot(data,
                              species_name = sp8,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro8.pdf",8.85,12.1)
plot(plot8)
dev.off()

# sequence 9
plot9 <- circular_linear_plot(data,
                              species_name = sp9,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro9.pdf",8.85,12.1)
plot(plot9)
dev.off()

# sequence 10
plot10 <- circular_linear_plot(data,
                              species_name = sp10,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro10.pdf",8.85,12.1)
plot(plot10)
dev.off()

# sequence 11
plot11 <- circular_linear_plot(data,
                              species_name = sp11,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro11.pdf",8.85,12.1)
plot(plot11)
dev.off()

# sequence 12
plot12 <- circular_linear_plot(data,
                              species_name = sp12,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro12.pdf",8.85,12.1)
plot(plot12)
dev.off()

# sequence 13
plot13 <- circular_linear_plot(data,
                              species_name = sp13,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro13.pdf",8.85,12.1)
plot(plot13)
dev.off()

# sequence 14
plot14 <- circular_linear_plot(data,
                              species_name = sp14,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro14.pdf",8.85,12.1)
plot(plot14)
dev.off()

# sequence 15
plot15 <- circular_linear_plot(data,
                              species_name = sp15,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro15.pdf",8.85,12.1)
plot(plot15)
dev.off()

# sequence 16
plot16 <- circular_linear_plot(data,
                              species_name = sp16,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro16.pdf",8.85,12.1)
plot(plot16)
dev.off()

# sequence 17
plot17 <- circular_linear_plot(data,
                              species_name = sp17,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro17.pdf",8.85,12.1)
plot(plot17)
dev.off()

# sequence 18
plot18 <- circular_linear_plot(data,
                              species_name = sp18,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro18.pdf",8.85,12.1)
plot(plot18)
dev.off()

# sequence 19
plot19 <- circular_linear_plot(data,
                              species_name = sp19,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro19.pdf",8.85,12.1)
plot(plot19)
dev.off()

# sequence 20
plot20 <- circular_linear_plot(data,
                              species_name = sp20,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro20.pdf",8.85,12.1)
plot(plot20)
dev.off()

# sequence 21
plot21 <- circular_linear_plot(data,
                              species_name = sp21,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro21.pdf",8.85,12.1)
plot(plot21)
dev.off()

# sequence 22
plot22 <- circular_linear_plot(data,
                              species_name = sp22,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro22.pdf",8.85,12.1)
plot(plot22)
dev.off()

# sequence 23
plot23 <- circular_linear_plot(data,
                              species_name = sp23,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro23.pdf",8.85,12.1)
plot(plot23)
dev.off()

# sequence 24
plot24 <- circular_linear_plot(data,
                              species_name = sp24,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro24.pdf",8.85,12.1)
plot(plot24)
dev.off()

# sequence 25
plot25 <- circular_linear_plot(data,
                              species_name = sp25,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro25.pdf",8.85,12.1)
plot(plot25)
dev.off()

# sequence 26
plot26 <- circular_linear_plot(data,
                              species_name = sp26,
                              leg_gradient = c(0,0.5,1),
                              title_name = "")
pdf("~/Desktop/repro26.pdf",8.85,12.1)
plot(plot26)
dev.off()
