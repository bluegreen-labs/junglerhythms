#' Extract line sections from annotations
#'
#' @param subset subset of annotations
#' @param img_y image height
#' @param bb bounding box coordinates from bounding_box()
#' @export
#' @return returns a matrix of line section locations (empty 0 if nothing returned)

extract_line_sections <- function(subset,
                                  bb,
                                  img_y){

  # create output matrix containing the 0/1 values
  # of observations during the year
  observations = matrix(0,4,336)

  if (!is.null(bb)){

    # calculate the approximate location of the rows, based upon
    # the corners of the yearly section
    row_coordinates = row_locations(bb)
    t_start = row_coordinates[[1]]
    t_end = row_coordinates[[2]]
    d_start = row_coordinates[[3]]
    d_end = row_coordinates[[4]]

    # get the length of the list of classifications
    # made for this particular image
    l <- length(subset)

    T2 <- grepl("T2", subset)

    # fill the coordinate matrix with data
    for (j in 1:l){

      # skip if no annotations are recorded
      check = length(unlist(subset[[j]]$value[3]))

      if (check == 0 | T2[j] == FALSE){
        next
      }

      # split out raw annotation data
      raw_annot <- subset[[j]]$value[3][[1]]

      # if there are no lines skip
      if(!any(which(raw_annot$tool %in% c(0,1)))){
        next
      }

      # get the labels of the tool used
      tool_labels <- subset[[j]]$value[3][[1]]$tool_label

      # split out line sections
      line_sections <- raw_annot[which(raw_annot$tool %in% c(0,1)),
                                c("x1","x2","y1","y2")]

      line_sections[,grep("y",names(line_sections))] <- img_y -
        line_sections[,grep("y",names(line_sections))]

      # check if there are coordinates
      if (nrow(line_sections)!=0){

        # loop over rows
        for (z in 1:4){

          # split in first and second half year
          for (c in c('start','end')){

            if (c == 'start'){
              coord1 = t_start[z,1:2] # check order
              coord2 = t_start[z,3:4]
            }else{
              coord1 = t_end[z,1:2]
              coord2 = t_end[z,3:4]
            }

            for (r in 1:nrow(line_sections)){

              # grab the coordinates of the projection of
              # the points defining a line onto the perfect
              # line

              # first coordinate distance
              d1 = dist_point_to_line(coord1,
                                      coord2,
                                      line_sections[r,c(1,3)])

              # second coordinate
              d2 = dist_point_to_line(coord1,
                                      coord2,
                                      line_sections[r,c(2,4)])

              # if the return is NA skip
              if(any(is.na(c(d1,d2)))){
                #print("no distance calculated")
                next
              }

              # if the coordinate is within 1/2 the inter line distance
              # I retain the point (if not, skip)
              if(d1 > (d_start/2) & d2 > (d_start/2)){
                #print("value out of range")
                next # kick out the values that are out of range
              }else{

                p1 = proj_point_on_line(coord1,
                                        coord2,
                                        line_sections[r,c(1,3)])
                p2 = proj_point_on_line(coord1,
                                        coord2,
                                        line_sections[r,c(2,4)])

                # if the return is NA skip
                if(any(is.na(c(p1,p2)))){
                  next
                }

                # swap the order of the coordinates if not in x sequence
                # p1 should be smaller than p2
                if (p1[1] > p2[1]){
                  tmp1 = p1
                  tmp2 = p2
                  p1 = tmp2
                  p2 = tmp1
                }

                # skip lines which are outside the first segment
                if (p1[1] > coord2[1] | p2[1] < coord1[1]){
                  next
                }

                # limit the range of the markings
                if (p1[1] < coord1[1]){
                  p1 = coord1
                }

                if (p2[1] > coord2[1]){
                  p2 = coord2
                }

                # calculate the index along the matrix in doy
                td = dist_point_to_point(coord1,coord2)
                pp1 = dist_point_to_point(coord1,p1)/td * 167
                pp2 = dist_point_to_point(coord1,p2)/td * 167

                # correction for second half of the year
                if (c == 'end'){
                  pp1 = pp1 + 168
                  pp2 = pp2 + 168
                }

                # round data to the nearest (DOY)
                pp1 = round(pp1)
                pp2 = round(pp2)

                tmp_marks = matrix(0,1,336)
                tmp_marks[pp1:pp2] = 1
                observations[z,] = observations[z,] + tmp_marks

              }
            }
          }
        }
      }
    }
  }

  # asign column names
  colnames(observations) <- paste("doy",1:ncol(observations),sep="_")

  # return observations (as data frame)
  return(as.data.frame(observations))
}
