#' Calculates event length and timing
#'
#' This sub-routine will take vectors from a data frame
#' for a single individual and phenophase and returns the timing of
#' phenophases and the length of the .
#'
#' @param year year vector
#' @param week week vector
#' @param value vector with observations
#' @export
#' @return timing of when an event (value) switches between no event (0)
#' and a observation of an event (1)

event_length <- function(df){

  # full frame
  full_df <- data.frame(year = sort(rep(1930:1965, 48)),
                        week = rep(1:48, 36),
                        value = rep(0, 36))

  full_df$date <- paste(full_df$year,
                        full_df$week,
                        sep = "-")

  # date vector from provided data frame
  date_df <- paste(df$year, df$week, sep = "-")

  # correction
  loc <- which(full_df$date %in% date_df)

  if(length(loc)!= length(df$value)){
    print(unique(df$image))
    print(unique(df$id))
  }

  full_df$value[loc] <- df$value

  # get first differences
  diff_values <- diff(full_df$value)

  # get matching info
  start <- full_df[which(diff_values == 1) + 1,]
  end <- full_df[which(diff_values == -1),]
  start$index <- as.numeric(rownames(start))
  end$index <- as.numeric(rownames(end))

  # the index is sequential, so take the difference
  # of this index, which counts weeks regardless
  # of date line transitions
  phenophase_length <- end$index - start$index + 1

  # return data as data frame
  return(data.frame(year_start = start$year,
                    week_start = start$week,
                    year_end = end$year,
                    week_end = end$week,
                    index_start = start$index,
                    index_end = end$index,
                    phenophase_length = phenophase_length))
}

library(tidyverse)

# read in the weekly data
df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
df <- df[which(df$value != 0),]

test <- df %>%
  group_by(genus, species, id, phenophase) %>%
  do(event_length(.))
print(test)
