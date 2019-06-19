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

metadata <- read.csv("data/phenology_archives_species_long_format_20190619.csv",
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

# specieslist based on selection made in advance
sp_evergreen <- c("Scorodophloeus zenkeri",
                  "Panda oleosa",
                  # "Anonidium mannii",
                  "Staudtia kamerunensis",
                  # "Pancovia harmsiana",
                  # "Garcinia punctata",
                  # "Carapa procera",
                  "Synsepalum subcordatum",
                  # "Pycnanthus angolensis",
                  "Prioria oxyphylla")
                  # "Blighia welwitschii",
                  # "Dacryodes edulis",
                 # "Phyllocosmus africanus")
sp_deciduous <- c("Petersianthus macrocarpus",
                  # "Greenwayodendron suaveolens",
                  "Irvingia grandifolia",
                  "Erythrophleum suaveolens",
                  "Autranella congolensis",
                  "Pericopsis elata")
                  # "Vitex ferruginea",
                  # "Ricinodendron heudelotii",
                  # "Milicia excelsa",
                  # "Lovoa trichilioides")

query_evergreen <- paste(sp_evergreen,collapse = "|")
query_deciduous <- paste(sp_deciduous,collapse = "|")

p_deciduous <- circular_linear_plot(data,
                                    species_name = query_deciduous,
                                    viridis_rescaling = 0.15,
                                    title_name = "(b) Deciduous")
pdf("~/Desktop/deciduous.pdf",8.4,12)
plot(p_deciduous)
dev.off()

p_evergreen <- circular_linear_plot(data,
                                    species_name = query_evergreen,
                                    viridis_rescaling = 0.05,
                                    title_name = "(a) Evergreen")
pdf("~/Desktop/evergreen.pdf",8.4,12)
plot(p_evergreen)
dev.off()

# p_circular_evergreen <- circle_plot(data, species_name = query_evergreen)
# pdf("~/Desktop/evergreen_selection.pdf",10,20)
# plot(p_circular_evergreen)
# dev.off()
#
# p_circular_deciduous <- circle_plot(data, species_name = query_deciduous)
# pdf("~/Desktop/deciduous_selection.pdf",10,20)
# plot(p_circular_deciduous)
# dev.off()
#
# # query_evergreen1 <- paste(sp_evergreen[c(1:8)],collapse = "|")
# # query_evergreen2 <- paste(sp_evergreen[c(8:14)],collapse = "|")
# # query_deciduous1 <- paste(sp_deciduous[c(1:6)],collapse = "|")
# # query_deciduous2 <- paste(sp_deciduous[c(6:10)],collapse = "|")
#
# p_linear_evergreen <- linear_plot_senescence(data = data, species_name = query_evergreen)
# pdf("~/Desktop/evergreen_selection_linear.pdf",5,8)
# plot(p_linear_evergreen)
# dev.off()
#
# p_linear_deciduous <- linear_plot_senescence(data = data, species_name = query_deciduous)
# pdf("~/Desktop/deciduous_selection_linear.pdf",5,8)
# plot(p_linear_deciduous)
# dev.off()


