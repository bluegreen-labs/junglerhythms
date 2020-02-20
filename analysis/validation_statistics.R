library(tidyverse)
library(cowplot)

# install inter rater reliability for validation stats
library(irr)

# read the data
yg <- as.tibble(readRDS("data/jungle_rhythms_weekly_annotations.rds"))
val <- readRDS("data/yangambi_validation_annotations.rds")

# convert data type for correct join
yg$image <- as.character(yg$image)

# join the data in one data frame
df <- left_join(val, yg, by = c("id", "phenophase", "year", "week", "image"))

# combine the ratings into one matrix
ratings <- na.omit(cbind(df$value.x, df$value.y))

# On the matrix, calculate the
# agreement and kappa indices for two raters
# 1. the zooniverse team
# 2. the validation team (gold standard)
kappa2(ratings)
agree(ratings)
