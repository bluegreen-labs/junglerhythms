#' distance between a point and a line (orthogonal) (2D)
#'
#' @param coord1 start line coordinate
#' @param coord2 end line coordinate
#' @param coord3 point coordinate
#' @export
#' @return distance
#'

dist_point_to_line = function(coord1, coord2, coord3){

  # grab the coordinates, and put them in
  # more readable variables
  x1 = coord1[1]
  y1 = coord1[2]
  x2 = coord2[1]
  y2 = coord2[2]
  x0 = coord3[1]
  y0 = coord3[2]

  # calculate the distance
  d = abs((y2 - y1)*x1 - (x2-x1)*y0 + x2*y1 - y2*x1)/sqrt((y2-y1)^2+(x2-x1)^2)

  # return the value
  return(d)
}

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
  x1 = coord1[1]
  y1 = coord1[2]
  x2 = coord2[1]
  y2 = coord2[2]

  d = sqrt((x2 - x1)^2 + (y2 - y1)^2)
  return(d)
}

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

#' line circle intersection (2D)
#' calculates line intersections with a
#' circle with origin at start, only
#' intersections with y > start_y are
#' reported
#'
#' @param start location of the center of a circle and the start of a line
#' @param end location of the end of a line intersecting a circle
#' @param radius radius of the circle to consider
#' @export
#' @return coordinates
#'

intersect_line_circle = function(start, end, radius){

  # circle variables
  c = start[1]
  d = start[2]
  r = radius

  # line equation / variables
  m = (start[2] - end[2])  / (start[1] - end[1])

  # if the slope is infinite
  # x = #
  if(is.infinite(m)){
    x1 = start[1]
    x2 = start[1]
    y1 = start[2] + r
    y2 = start[2] - r
  }else{
    b = line_intersection(coord1 = start,
                          coord2 = end,
                          coord3 = c(0,0),
                          coord4 = c(0,1))[2] # grab intercept

    x1 = (-sqrt(-(b^2)-(2*b*c*m)+(2*b*d)-(c^2)*(m^2)+(2*c*d*m)-(d^2)+(m^2)*(r^2)+r^2)-b*m+c+d*m)/(m^2+1)
    x2 = (sqrt(-(b^2)-(2*b*c*m)+(2*b*d)-(c^2)*(m^2)+(2*c*d*m)-(d^2)+(m^2)*(r^2)+r^2)-b*m+c+d*m)/(m^2+1)
    y1 = m*x1 + b
    y2 = m*x2 + b
  }

  # put in one data.frame
  x = rbind(x1,x2)
  y = rbind(y1,y2)
  xy = data.frame(x,y)

  # subset high value
  xy = xy[which(xy$y >= d),]
  return(xy)
}

#' orthogonal projection of a point onto a line (2D)
#'
#' @param coord1 start coordinate of a line
#' @param coord2 end coordinate fo a line
#' @param coord3 coordinates of a point
#' @export
#' @return coordinates x,y

proj_point_on_line = function(coord1, coord2, coord3){

  # grab the coordinates, and put them in
  # more readable variables
  x1 = coord1[1]
  y1 = coord1[2]
  x2 = coord2[1]
  y2 = coord2[2]
  x0 = coord3[1]
  y0 = coord3[2]

  # calculate slope
  m = (y2 - y1) / (x2 - x1)

  # calculate intercept
  b = y1 - (m * x1)

  # calculate projected coordinates
  x = (m * y0 + x0 - m * b)/(m^2+1)
  y = (m^2 * y0 + m*x0 + b)/(m^2+1)

  # return the value
  return(data.frame("x" = x,
                    "y" = y))
}
