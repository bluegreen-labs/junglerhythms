#' Reads in Jungle Rhythms classification
#' output file (CSV) and cleans up the data
#' by dropping early values (trials) etc
#'
#' Uses data.table for speed given the large
#' size of the data (700MB).
#'
#' @param file: path to a zooniverse generated CSV
#' @param date: date before which to discard all data (demo data)
#' @keywords data, io, transformation
#' @export

read_zoo_data <- function(file = NULL,
                          date = "2015-12-13"){

  # check if file exists
  # if so read in the data
  if (file.exists(file)){
    dt <- fread(file)
  } else {
    stop("file does not exist, check path!")
  }

  # subset the data based upon a date
  # 13 Dec. was the launch date of the project
  dt <- dt[which(as.Date(dt$created_at)>as.Date(date)),]

  # return the data table
  return(dt)
}
