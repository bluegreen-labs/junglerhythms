#' point to point distance (pythagoras theorem) (2D)


dist_point_to_point = function(coord1,coord2){

  # grab the coordinates, and put them in
  # more readable variables
  x1 = coord1[1] ; y1 = coord1[2]
  x2 = coord2[1] ; y2 = coord2[2]

  d = sqrt((x2 - x1)^2 + (y2 - y1)^2)
  return(d)
}
