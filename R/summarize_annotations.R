#' Summarize annotations into weekly values
#'
#' weekly values assume a 48 week year
#'
#' @param df junglerhythms raw annotation file
#' @param image_index image index data meta-data
#' @param plot plot output TRUE or FALSE
#' @param internal output results to R workspace, TRUE or FALSE
#' @param output_path output path where to store data
#' @export
#' @return ggplot object

summarize_annotations <- function(
  df = "data/jungle_rhythms_daily_annotations.rds",
  image_index = "data/phenology_archives_species_long_format_20190319.csv",
  plot = TRUE,
  internal = FALSE,
  output_path = "./data/"
  ){

  # read in data directly from file path if
  # not a data frame
  if(!is.data.frame(df)){
    if (file.exists(df)){
    df <- readRDS(df)
    }
  }

  # read in phenology archive image index
  index <- read.csv2(image_index,
                      sep = ",",
                      header = TRUE,
                      stringsAsFactors = FALSE)

  # grab image name
  image_names <- paste(df$image_nr,
                       df$image_row,
                       df$image_col)

  # progress bar settings
  pb <- txtProgressBar(0, length(unique(image_names)), style = 3)
  env <- environment()
  i = 1

  #output <- do.call("rbind",lapply("1066198 6 2",
  output <- do.call("rbind",lapply(unique(image_names),
                                   function(image_name){

    # progress bar
    tmp <- get("i", envir = env)
    tmp <- tmp + 1
    assign("i", tmp, envir = env)
    setTxtProgressBar(pb, tmp)

    # subset the data based upon image name
    # includes all phenology observations
    subset <- df[which(image_names == image_name),]

    # get labels
    labels <- subset$labels

    # loop over all phenology phenophases
    values <- do.call(
      "cbind",
        lapply(labels, function(label){

      # subset based upon label and doy
      x <- subset[which(subset$labels == label),
                  grep("doy",names(df))]

      # convert to numeric
      x <- as.vector(apply(x, 1:2, as.numeric))

      if(all(x == 0)){
        return(rep(0,48))
      }

      # calculate relative majority vote
      # assumes "unity", the assumption that if there
      # is a visible marking all people will mark it with
      # a varying degree of accuracy for a particular row (phenophase)
      v <- floor(max(x, na.rm = TRUE) * 0.5)

        if(v <= 2 | is.na(v)){
          x[] <- NA
        } else {
          x[x<v] <- NA
          x[!is.na(x)] <- 1
        }

      # rescale to 48 week year from a 336 DOY notation
      weeks <- sort(rep(1:48,7))

      x_week <- as.vector(aggregate(unlist(x),
                          by = list(weeks),
                          FUN = function(x_subset){
        c <- length(which(!is.na(x_subset)))
        ifelse(c >= 3, 1, 0)
      }))$x

      # return data
      return(x_week)
    }))

  values <- as.data.frame(values)
  colnames(values) <- labels

  # find meta-data
  loc <- as.numeric(unlist(strsplit(image_name," ")))
  image_nr <- loc[1]
  image_row = loc[2]
  image_col = loc[3]

  # define final location of meta-data
  loc <- which(grepl(image_nr, index$image) &
                 index$row == image_row)

  # skip if no meta-data if found (empty slots)
  if(nrow(index[loc,]) == 0){
    return(NULL)
  }

  # grab starting year
  starting_year <- as.numeric(index$starting_year[loc])

  # correct starting year for double image tags
  img_tag <- unlist(strsplit(index$image[loc],"/"))
  starting_year <- 10 * (grep(image_nr, img_tag) - 1) + starting_year

  # collate meta-data
  values$year <- starting_year + (image_col - 1)
  values$week <- 1:48
  values$family <- index$family_plantlist[loc]
  values$genus <- index$genus_plantlist[loc]
  values$species <- index$species_plantlist[loc]
  values$plantlist_status <- index$plantlist_status[loc]
  values$image <- image_nr
  values$image_col <- image_col
  values$image_row <- image_row
  values$id <- index$id[loc]

  # check if there are annotated values, even if
  # thresholds are not reached
  x <- as.numeric(unlist(subset[,grep("doy",names(subset))]))
  values$annotations <- ifelse(any(x > 0), TRUE, FALSE)

  return(values)
  }))

  # close progress bar
  close(pb)

  # convert to a tidy format
  output <- tidyr::gather(output,
                key = phenophase,
                value = value,
                -year,
                -week,
                -genus,
                -species,
                -family,
                -image,
                -image_col,
                -image_row,
                -id,
                -plantlist_status,
                -annotations)

  # return data
  if(!internal){
    saveRDS(output, paste0(output_path, "/jungle_rhythms_weekly_annotations.rds"))
  } else {
    return(output)
  }
}
