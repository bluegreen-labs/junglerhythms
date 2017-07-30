#' orthogonal projection of a point onto a line (2D)
#'
#'
#'

proj_point_on_line = function(coord1, coord2, coord3){

  # grab the coordinates, and put them in
  # more readable variables
  x1 = coord1[1] ; y1 = coord1[2]
  x2 = coord2[1] ; y2 = coord2[2]
  x0 = coord3[1] ; y0 = coord3[2]

  # calculate slope
  m = (y2 - y1) / (x2 - x1)

  # calculate intercept
  b = y1 - (m * x1)

  # calculate projected coordinates
  x = (m * y0 + x0 - m * b)/(m^2+1)
  y = (m^2 * y0 + m*x0 + b)/(m^2+1)

  # return the value
  return(cbind(x,y))
}
