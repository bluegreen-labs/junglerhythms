#' calculate the median bounding box for a yearly section
#'
#' Uses the corners of the yearly section as outlined by
#' volunteers to establish the 6 coordinate points making up
#' the bounding box. This is done for a subset of the full
#' annotated dataset.
#'
#' @param subset subset of annotations
#' @param img_y image height
#' @export
#' @return median coordinates defining the bounding box (6 pairs)

bounding_box <- function(subset,
                         img_y){

  # get the length of the list of classifications
  # made for this particular image
  l <- length(subset)

  T2 <- grepl("T2", subset)

  # fill the coordinate matrix with data
  output <- do.call("rbind",
    lapply(1:l, function(j){

      # are the annotations empty or not
      # exit checks
      check = length(unlist(subset[[j]]$value[2]))

      if (check == 0 | T2[j] == FALSE){
        return(NULL)
      }

      if (grepl("Yes",subset[[j]]$value[1])){

        # grab the bounding box coordinates
        # and correct the y-axis for the size of the
        # image
        bb <- do.call("rbind",
                      subset[[j]]$value[2][[1]]$points)
        bb$y <- img_y - bb$y

        # order coordinates by x/y values
        bb <- bb[order(bb$x), ]

        # number of rows in the bounding box
        # coordinates matrix
        bb_nrow <- length(bb$y)

        # now check for size
        # if there are only 4 we can still use the data
        # discard anything with 3 or less
        if (bb_nrow >= 4){

          # cal_yscul_ysate top l_yseft bottom right, write dedicated routine for this
          bottom_left = bb[which(bb$y[1:2] == min(bb$y[1:2] )),]
          top_left = bb[which(bb$y[1:2] == max(bb$y[1:2] )),]

          top_right = bb[which(bb$y[c(bb_nrow-1,bb_nrow)] == max(bb$y[c(bb_nrow-1,bb_nrow)]))+(bb_nrow-2),]
          bottom_right = bb[which(bb$y[c(bb_nrow-1,bb_nrow)] == min(bb$y[c(bb_nrow-1,bb_nrow)]))+(bb_nrow-2),]

          # calculate middle points, in case of 5 just use the corners
          # ignore the middle one
          if  (bb_nrow == 6){
            top_middle = bb[which(bb$y[c(bb_nrow-3,bb_nrow-2)] == max(bb$y[c(bb_nrow-3,bb_nrow-2)]))+(bb_nrow-4),]
            bottom_middle = bb[which(bb$y[c(bb_nrow-3,bb_nrow-2)] == min(bb$y[c(bb_nrow-3,bb_nrow-2)]))+(bb_nrow-4),]
          }else{
            top_middle = c(NA,NA)
            bottom_middle = c(NA,NA)
          }

          # return data vector
          v = c(as.matrix(bottom_left),
                as.matrix(top_left),
                as.matrix(bottom_middle),
                as.matrix(top_middle),
                as.matrix(bottom_right),
                as.matrix(top_right))

          return(v[1:12])
        }
      }
    }))

  if(is.null(output)){
    return(NULL)
  }

  # calculate the median coordinates
  median_coord <- apply(output, 2, median, na.rm=T)

  # return values
  return(median_coord)
}
