#' distance between a point and a line (orthogonal) (2D)


dist_point_to_line = function(coord1, coord2, coord3){

  # grab the coordinates, and put them in
  # more readable variables
  x1 = coord1[1] ; y1 = coord1[2]
  x2 = coord2[1] ; y2 = coord2[2]
  x0 = coord3[1] ; y0 = coord3[2]

  # calculate the distance
  d = abs((y2 - y1)*x1 - (x2-x1)*y0 + x2*y1 - y2*x1)/sqrt((y2-y1)^2+(x2-x1)^2)

  # return the value
  return(d)
}
