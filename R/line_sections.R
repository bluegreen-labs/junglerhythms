#' Extract line sections from annotations
#'
#' This is the main routine used in the extraction of
#' data from the Zooniverse data frame (converted JSON or
#' nested list). The input consist of all classifications
#' for a particular subset, the bounding box associated
#' with these classications (based upon the original image)
#' and the size of the image along the y-axis for referencing
#' and scaling of the data.
#'
#' NOTE: This might be a bit less unwieldy in tidyverse
#' syntax but it isn't worth the effort to rewrite it.
#' Plenty of exceptions might also get in the way of
#' holding true.
#'
#' @param subset subset of annotations
#' @param img_y image height
#' @param bb bounding box coordinates from bounding_box()
#' @export
#' @return returns a matrix of line section locations
#' (empty 0 if nothing returned)

line_sections <- function(subset,
                                  bb,
                                  img_y){

  # Create an output matrix containing the values
  # of observations during the year, in this case
  # I specify a year with 48 weeks, corresponding
  # to the format used in the historical data
  # # doy = 7 x 48. This value can be changed,
  # but resampling is easier in this format.
  observations <- matrix(0,5,336)

  # Check if there is a valid bounding box + inter row
  # width values. This data acts as a reference to
  # calculate the relative length of the annotated
  # sections
  if (is.null(bb)){
    return(NULL)
  }

  # If there is a bounding box, calculate the approximate
  # location of the rows, based upon the corners of the
  # yearly section and "ideal" inter row distance
  row_coordinates = row_locations(bb)

  # spit out the coordinates and row distance (d_)
  # for further processing
  t_start = row_coordinates[[1]]
  t_end = row_coordinates[[2]]
  d_start = row_coordinates[[3]]
  d_end = row_coordinates[[4]]

  # Get the length of the list of classifications
  # made for this particular image. I'll itterate
  # over these values for processing, evaluating
  # every classification
  l <- length(subset)

  # Grab the answers to question two for all
  # classifications. This asks if there are
  # visible marks in the yearly section, if FALSE
  # the processing will be skipped (see below)
  T2 <- grepl("T2", subset)

  # MAIN LOOP: working through all classifications
  # enumerated from 1:l (skip if not values are present)
  for(j in 1:l){

    # Skip if no annotations are recorded, yet
    # the answer to question 2 (T2) was TRUE
    if (length(unlist(subset[[j]]$value[3])) == 0){
      next
    }

    # Skip if the answer to question two was FALSE
    if (T2[j] == FALSE){
      next
    }

    # Skip if there are not annotations for lines
    # (e.g. a marking of the end of the series but
    # not phenology anymore)
    raw_annot <- subset[[j]]$value[3][[1]]
    if(!any(which(raw_annot$tool %in% c(0,1)))){
      next
    }

    # Split out all tool labels. Will be used below
    # for processing of the crosshatched lines in
    # particular.
    tools <- raw_annot$tool_label

    # Split out the line sections, without meta-data
    # including only the coordinates for further processing.
    line_sections <- raw_annot[, c("x1","x2","y1","y2")]
    line_sections[,grep("y",names(line_sections))] <- img_y -
      line_sections[,grep("y",names(line_sections))]

    # Check if there are coordinates available and
    # trap any empty arrays.
    if (nrow(line_sections) == 0){
      next
    }

    # Loop over the 4 rows with values
    # i.e. phenology observations
    for (z in 1:4){

      # To account for skewness and warping of the
      # paper data is processed on a half yearly basis
      # evaluating the start and end separately.
      for (c in c('start','end')){

        # Set bounding box coordinates of the start or
        # end section for post processing.
        if (c == 'start'){
          coord1 = t_start[z,1:2]
          coord2 = t_start[z,3:4]
        }else{
          coord1 = t_end[z,1:2]
          coord2 = t_end[z,3:4]
        }

        # Loop over all available line sections and
        # project the line section on the ideal line
        # taking into account the distance to this ideal
        # line (too far will result in a rejection of the value).
        for (r in 1:nrow(line_sections)){

          # Grab the coordinates of the projection of
          # the points defining a line onto the ideal line.
          # Do this for both the beginning and end points.
          d1 = dist_point_to_line(coord1,
                                  coord2,
                                  line_sections[r,c(1,3)])

          d2 = dist_point_to_line(coord1,
                                  coord2,
                                  line_sections[r,c(2,4)])

          # If the routine is not valide (NA)
          # skip to the next line section.
          if(any(is.na(c(d1,d2)))){
            next
          }

          # If the coordinates are within 1/2 the inter-line distance
          # retain the point (if not, skip).
          if(d1 > (d_start/2) & d2 > (d_start/2)){
            next
          }else{

            # Project the retained points onto the ideal row
            # location.
            p1 = proj_point_on_line(coord1,
                                    coord2,
                                    line_sections[r,c(1,3)])
            p2 = proj_point_on_line(coord1,
                                    coord2,
                                    line_sections[r,c(2,4)])

            # Skip if this fails, should not occur often
            # if at all.
            if(any(is.na(c(p1,p2)))){
              next
            }

            # Swap the order of the coordinates if not in x sequence
            # p1 should be smaller than p2. The order of the points
            # defining the line can be arbitrary.
            if (p1[1] > p2[1]){
              tmp1 = p1
              tmp2 = p2
              p1 = tmp2
              p2 = tmp1
            }

            # Skip lines which are outside the other segment
            if (p1[1] > coord2[1] | p2[1] < coord1[1]){
              next
            }

            # Limit the range of the markings to the extent of
            # the bounding box (or half-yearly section in this case).
            if (p1[1] < coord1[1]){
              p1 = coord1
            }

            if (p2[1] > coord2[1]){
              p2 = coord2
            }

            # After trimming the points we can calculate the
            # distance between the two points and map this to
            # a DOY values.
            td = dist_point_to_point(coord1,coord2)
            pp1 = dist_point_to_point(coord1,p1)/td * 167
            pp2 = dist_point_to_point(coord1,p2)/td * 167

            # Correct values for second half of the year
            if (c == 'end'){
              pp1 = pp1 + 168
              pp2 = pp2 + 168
            }

            # Round data to the nearest (DOY)
            pp1 = as.vector(as.matrix(round(pp1)))
            pp2 = as.vector(as.matrix(round(pp2)))

            # Create a temporary matrix of values
            # and fill those values at the respective
            # DOY with a value of 1.
            tmp_marks = matrix(0,1,336)
            tmp_marks[pp1:pp2] = 1

            # Depending on the tool used add these values
            # to the lines corresponding to their respective
            # phenological event. This creates a cummulative
            # sum of all annotations.
            if (tools[r] == "crosshatched lines"){
              observations[5,] = observations[5,] + tmp_marks
            } else {
              observations[5-z,] = observations[5-z,] + tmp_marks
            }
          }
        }
      }
    }
  }

  # Reverse the order of rows to correspond to this of the images
  # to make visualizations easier on a row by row basis.
  observations <- observations[nrow(observations):1,]

  print(str(observations))

  # Assign numbered column names
  colnames(observations) <- paste("doy",1:ncol(observations),sep="_")

  # add labels to the file
  labels <- rev(c("flowers",
                  "fruit",
                  "fruit_drop",
                  "leaf_dormancy",
                  "leaf_turnover"))

  # return the data
  return(cbind(labels,observations))
}
