#----- reset your R session. ------------------------------------------
rm(list=ls())
graphics.off()
#----- load required packages -----------------------------------------
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(gridExtra)
#----- source required functions --------------------------------------
source("R/circular_linear_plots.R")
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#--------   Phenology data --------------------------------------------
#----------------------------------------------------------------------
data <- readRDS("data/jungle_rhythms_data_Pierlot_species-96.rds")
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#--- get selected species DECIDUOUS label
#----------------------------------------------------------------------
overview <- read.csv("data/Pierlot_summ_species_pheno_characteristics.csv",
                     header = TRUE,
                     sep = ",",
                     stringsAsFactors = FALSE)

overview_dec <- overview %>%
  filter(grepl("deciduous",deciduousness)) %>%
  filter(!species_full == 'Alstonia boonei') %>% # this species messes up the color scale for the others. Is done separately
  filter(site_years_with_phenophase > 0) %>%
  arrange_at("species_full")



sp_dec <- unique(overview_dec$species_full)

# groups of 5 species (5 species per page in the supplementary info)
sp_dec1 <- paste(sp_dec[1:5],collapse = "|")
sp_dec2 <- paste(sp_dec[6:10],collapse = "|")
sp_dec3 <- paste(sp_dec[11:15],collapse = "|")
sp_dec4 <- paste(sp_dec[16:20],collapse = "|")
sp_dec5 <- paste(sp_dec[21:25],collapse = "|")
sp_dec6 <- paste(sp_dec[26:30],collapse = "|")
sp_dec7 <- paste(sp_dec[31:35],collapse = "|")
sp_dec8 <- paste(sp_dec[36:39],collapse = "|")
sp_dec9 <- 'Alstonia boonei'

#------
# plots
#------
# sequence 1
plot_dec1 <- circular_linear_plot(data = data,
                                  species_name = sp_dec1,
                                  title_name = "(b) Deciduous")
pdf("~/Desktop/deciduous1.pdf",8.85,12.1)
plot(plot_dec1)
dev.off()

# sequence 2
plot_dec2 <- circular_linear_plot(data = data,
                                  species_name = sp_dec2,
                                  title_name = "(b) Deciduous - continued")
pdf("~/Desktop/deciduous2.pdf",8.85,12.1)
plot(plot_dec2)
dev.off()

# sequence 3
plot_dec3 <- circular_linear_plot(data = data,
                                  species_name = sp_dec3,
                                  title_name = "(b) Deciduous - continued")
pdf("~/Desktop/deciduous3.pdf",8.85,12.1)
plot(plot_dec3)
dev.off()

# sequence 4
plot_dec4 <- circular_linear_plot(data = data,
                                  species_name = sp_dec4,
                                  title_name = "(b) Deciduous - continued")
pdf("~/Desktop/deciduous4.pdf",8.85,12.1)
plot(plot_dec4)
dev.off()

# sequence 5
plot_dec5 <- circular_linear_plot(data = data,
                                  species_name = sp_dec5,
                                  title_name = "(b) Deciduous - continued")
pdf("~/Desktop/deciduous5.pdf",8.85,12.1)
plot(plot_dec5)
dev.off()

# sequence 6
plot_dec6 <- circular_linear_plot(data = data,
                                  species_name = sp_dec6,
                                  title_name = "(b) Deciduous - continued")
pdf("~/Desktop/deciduous6.pdf",8.85,12.1)
plot(plot_dec6)
dev.off()

# sequence 7
plot_dec7 <- circular_linear_plot(data = data,
                                  species_name = sp_dec7,
                                  title_name = "(b) Deciduous - continued")
pdf("~/Desktop/deciduous7.pdf",8.85, 12.1)
plot(plot_dec7)
dev.off()

# sequence 8
plot_dec8 <- circular_linear_plot(data = data,
                                  species_name = sp_dec8,
                                  leg_pos = c(1,0.25),  # changed legend position
                                  title_name = "(b) Deciduous - continued")
pdf("~/Desktop/deciduous8.pdf",8.85, 9.8)
plot(plot_dec8)
dev.off()

# sequence 9
plot_dec9 <- circular_linear_plot(data = data,
                                  species_name = sp_dec9,
                                  leg_pos = 'none', # changed legend position
                                  title_name = "(b) Deciduous - continued")
pdf("~/Desktop/deciduous9.pdf",8.85,3)
plot(plot_dec9)
dev.off()

#-------------------------------------------------------------------
#--- get selected species EVERGREEN
#-------------------------------------------------------------------
overview_ever <- overview %>%
  filter(grepl("evergreen",deciduousness)) %>%
  filter(site_years_with_phenophase > 0) %>%
  arrange_at("species_full")

sp_ever <- unique(overview_ever$species_full)

# groups of 5 species (5 species per page in the supplementary info)
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
sp_ever13 <- paste(sp_ever[60:63],collapse = "|")


#------
# plots
#------
# sequence 1
# changed legend position
plot_ever1 <- circular_linear_plot(data = data,
                                   species_name = sp_ever1,
                                   leg_pos = c(1,0.25),  # changed legend position
                                   title_name = "(a) Evergreens (with a recorded phenophase)")
pdf("~/Desktop/evergreen1.pdf",8.85,9.8)
plot(plot_ever1)
dev.off()

# sequence 2
plot_ever2 <- circular_linear_plot(data = data,
                                   species_name = sp_ever2,
                                   title_name = "(a) Evergreens - continued")
pdf("~/Desktop/evergreen2.pdf",8.85,12.1)
plot(plot_ever2)
dev.off()

# sequence 3
plot_ever3 <- circular_linear_plot(data = data,
                                   species_name = sp_ever3,
                                   title_name = "(a) Evergreens - continued")
pdf("~/Desktop/evergreen3.pdf",8.85,12.1)
plot(plot_ever3)
dev.off()

# sequence 4
plot_ever4 <- circular_linear_plot(data = data,
                                   species_name = sp_ever4,
                                   title_name = "(a) Evergreens - continued")
pdf("~/Desktop/evergreen4.pdf",8.85,12.1)
plot(plot_ever4)
dev.off()

# sequence 5
plot_ever5 <- circular_linear_plot(data = data,
                                   species_name = sp_ever5,
                                   title_name = "(a) Evergreens - continued")
pdf("~/Desktop/evergreen5.pdf",8.85,12.1)
plot(plot_ever5)
dev.off()

# sequence 6
plot_ever6 <- circular_linear_plot(data = data,
                                   species_name = sp_ever6,
                                   title_name = "(a) Evergreens - continued")
pdf("~/Desktop/evergreen6.pdf",8.85,12.1)
plot(plot_ever6)
dev.off()


# sequence 7
plot_ever7 <- circular_linear_plot(data = data,
                                   species_name = sp_ever7,
                                   title_name = "(a) Evergreens - continued")
pdf("~/Desktop/evergreen7.pdf",8.85,12.1)
plot(plot_ever7)
dev.off()

# sequence 8
plot_ever8 <- circular_linear_plot(data = data,
                                   species_name = sp_ever8,
                                   title_name = "(a) Evergreens - continued")
pdf("~/Desktop/evergreen8.pdf",8.85,12.1)
plot(plot_ever8)
dev.off()

# sequence 9
plot_ever9 <- circular_linear_plot(data = data,
                                   species_name = sp_ever9,
                                   title_name = "(a) Evergreens - continued")
pdf("~/Desktop/evergreen9.pdf",8.85,12.1)
plot(plot_ever9)
dev.off()

# sequence 10
plot_ever10 <- circular_linear_plot(data = data,
                                    species_name = sp_ever10,
                                    title_name = "(a) Evergreens - continued")
pdf("~/Desktop/evergreen10.pdf",8.85,12.1)
plot(plot_ever10)
dev.off()

# sequence 11
plot_ever11 <- circular_linear_plot(data = data,
                                    species_name = sp_ever11,
                                    title_name = "(a) Evergreens - continued")
pdf("~/Desktop/evergreen11.pdf",8.85,12.1)
plot(plot_ever11)
dev.off()

# sequence 12
plot_ever12 <- circular_linear_plot(data = data,
                                    species_name = sp_ever12,
                                    title_name = "(a) Evergreens - continued")
pdf("~/Desktop/evergreen12.pdf",8.85,12.1)
plot(plot_ever12)
dev.off()

# sequence 13
plot_ever13 <- circular_linear_plot(data = data,
                                    species_name = sp_ever13,
                                    leg_pos = c(1,0.25),  # changed legend position
                                    title_name = "(a) Evergreens - continued")
pdf("~/Desktop/evergreen13.pdf",8.85,9.8)
plot(plot_ever13)
dev.off()

#-------------------------------------------------------------------
#--- get selected species UNCLASSIFIED
#-------------------------------------------------------------------
overview_unclassed <- overview %>%
  filter(is.na(deciduousness)) %>%
  filter(site_years_with_phenophase > 0) %>%
  arrange_at("species_full")

unclassed <- unique(overview_unclassed$species_full)

# groups of 5 species (5 species per page in the supplementary info)
sp_unclassed1 <- paste(unclassed[1:4],collapse = "|") # up to 4 for the first page which includes figure description
sp_unclassed2 <- paste(unclassed[5:6],collapse = "|")


#------
# plots
#------
# sequence 1
# changed legend position
plot_unclassed1 <- circular_linear_plot(data = data,
                                   species_name = sp_unclassed1,
                                   leg_pos = c(1,0.25),
                                   title_name = "(c) Unclassified")
pdf("~/Desktop/unclassed1.pdf",8.85,9.8)
plot(plot_unclassed1)
dev.off()

# sequence 2
plot_unclassed2 <- circular_linear_plot(data = data,
                                   species_name = sp_unclassed2,
                                   leg_pos = c(1,0.5),
                                   title_name = "(c) Unclassified - continued")
pdf("~/Desktop/unclassed2.pdf",8.5,4.95)
plot(plot_unclassed2)
dev.off()




