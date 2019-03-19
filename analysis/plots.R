# query
# sp <- na.omit(read.table("~/Downloads/Specieslist.csv",header = TRUE, sep = ",",
#                          stringsAsFactors = FALSE))
# query <- paste(sp$Species,collapse = "|")

# # # read data
# data <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
# data$species_full <- paste(data$genus, data$species)
#
# # # test
# data <- data %>%
#   group_by(species_full, year) %>%
#   mutate(selection = if(all(value == 0)){
#    NA
#   }else{
#     value
#   }) %>%
#   na.omit() %>%
#   ungroup()

# # create plot
# p <- circle_plot_senescence(data = data,
#                             species_name = query,
#                             weeks = 48)
#
# pdf("~/Desktop/test.pdf",25,25)
#  plot(p)
# dev.off()

#data <- readRDS("data/luki_weekly_annotations.rds")

# create plot
p <- circle_plot_senescence(data = data,
                            species_name = "Carapa procera",
                            weeks = 36)

#pdf("~/Desktop/test.pdf",25,25)
plot(p)
#dev.off()

