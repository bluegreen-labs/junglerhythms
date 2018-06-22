#' Calculate line intersection when
#' provided four coordinaes defining two lines
#'
#' @param coord1: first coordinate line 1
#' @param coord2: second coordinate line 1
#' @param coord3: first coordinate line 2
#' @param coord4: second coordinate line 2
#'
#' @keywords data, transformation, geometry
#' @export

line_intersection = function(coord1 = NULL,
                             coord2 = NULL,
                             coord3 = NULL,
                             coord4 = NULL){

  if (any(is.null(c(coord1,coord2,coord3,coord4)))){
    stop("Four coordinates must be provided !")
  }

  # put coordinates in easier variable names
  x1 = coord1[1]; y1 = coord1[2]
  x2 = coord2[1]; y2 = coord2[2]
  x3 = coord3[1]; y3 = coord3[2]
  x4 = coord4[1]; y4 = coord4[2]

  x_int = ((x1 * y2 - y1 * x2)*(x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4))/((x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4))
  y_int = ((x1 * y2 - y1 * x2)*(y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4))/((x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4))

  # return coordinates of the intersection
  return(cbind(x_int,y_int))
}

