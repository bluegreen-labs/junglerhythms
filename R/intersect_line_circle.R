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
