#' point to point distance (pythagoras theorem) (2D)
#'
#' @param coord1 first point coordinate (x,y)
#' @param coord2 second point coordinate (x,y)
#' @export
#' @return distance between points

dist_point_to_point = function(coord1 = NULL,
                               coord2 = NULL){

  if (any(is.null(c(coord1,coord2)))){
    stop("Four coordinates must be provided !")
  }

  # grab the coordinates, and put them in
  # more readable variables
  x1 = coord1[1] ; y1 = coord1[2]
  x2 = coord2[1] ; y2 = coord2[2]

  d = sqrt((x2 - x1)^2 + (y2 - y1)^2)
  return(d)
}
