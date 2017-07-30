#' calculate the approximate location of the rows, based upon
#' the corners of the yearly section (provided as a vector)
#' in order: bottom_left, top_left, bottom_middle,top_middle,
#' bottom_right, top_right

row_locations = function(coord){

  # start of the year
  d_left_start = sqrt((coord[1]-coord[3])^2+(coord[2]-coord[4])^2)
  d_left_start = d_left_start/5
  d_right_start = sqrt((coord[5]-coord[7])^2+(coord[6]-coord[8])^2)
  d_right_start = d_right_start/5
  left_start = matrix(unlist(lapply(1:4 * d_left_start,function(x)linecirc(coord[1:2],coord[3:4],x))),4,2,byrow=T)
  right_start = matrix(unlist(lapply(1:4 * d_right_start,function(x)linecirc(coord[5:6],coord[7:8],x))),4,2,byrow=T)
  d_start = mean(c(d_left_start,d_right_start))

  # end of the year
  d_left_end = sqrt((coord[5]-coord[7])^2+(coord[6]-coord[8])^2)
  d_left_end = d_left_end/5
  d_right_end = sqrt((coord[9]-coord[11])^2+(coord[10]-coord[12])^2)
  d_right_end = d_right_end/5
  left_end = right_start
  right_end = matrix(unlist(lapply(1:4 * d_right_end,function(x)linecirc(coord[9:10],coord[11:12],x))),4,2,byrow=T)
  d_end = mean(c(d_left_end,d_right_end))

  # bind
  t_start = cbind(left_start,right_start)
  t_end = cbind(left_end,right_end)

  # output as nested list
  return(list(t_start,t_end,d_start,d_end))

}
