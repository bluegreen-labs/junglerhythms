#' Reads in Jungle Rhythms classification
#' output file (CSV) and cleans up the data
#' by dropping early values (trials) etc
#'
#' Uses data.table for speed given the large
#' size of the data (700MB).
#'
#' @param data: a vector with doy values 1 - 365(6)
#' @keywords data, io, transformation
#' @export
#' @examples
#'
#' \dontrun{
#' }


read_data <- function(file = NULL){

  # check if file exists
  # if so read in the data
  if (file.exists(file)){
    dt <- fread(file)
  } else {
    stop("file does not exist, check path!")
  }

  # subset the data based upon a date


}
