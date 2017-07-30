#' returns classification results for Jungle Rhythms
#' data

phenology_dates = function(data,
                           plot = FALSE,
                           image_path = NULL){

  # grab the annotation data from the annotation column
  subset = lapply(as.character(data$annotations),fromJSON)

  # check if have more than 5 annotations
  # if not skip
  l = length(subset)

  # if less than 10 annotations are made, skip
  if (l < 1){
    print("not enough data")
    return(NULL)
  }

  # are the annotations empty or not
  T2 = grepl("T2",subset)

  # If more than 20% is flagged as missing,
  # consider empty
  if(mean(T2) < 0.8){
    return(NULL)
  }

  # extract image name and calculate the y dimension
  subject_data = lapply(as.character(data$subject_data),fromJSON)
  image_name = unique(unlist(lapply(subject_data,function(x)unlist(x)[grep("#Filename",names(unlist(x)))])))
  r = brick(sprintf("%s/%s",image_path,image_name))
  size_y = dim(r)[1]

  # plot data if necessary
  # top image is the clean reference
  # bottom image is the 'annotated image'
  if (plot == TRUE){
    par(mfrow=c(3,1))
    plotRGB(r,main=image_name)
    plotRGB(r,main=image_name)
  }

  # create output matrix
  output = matrix(0,l,12)

  # fill the coordinate matrix with data
  for (j in 1:l){

    if (T2[j] == FALSE){
      next
    }

    if (grepl("Yes",subset[[j]]$value)){

      # grab the points marking the yearly section
      ys_p = subset[j][[3]]$value[[2]]$points
      l_ys = length(ys_p)

      ys_p = data.frame(matrix(unlist(ys_p),l_ys,2,byrow=T))
      colnames(ys_p) = c('x','y') # makes manipulations easier

      # grab the max y value to correct for the flipped
      # axis, otherwise use the provided value
      if(length(size_y)==0){
        size_y = max(ys_p$y)
      }

      ys_p[,2] = size_y - ys_p[,2]
      ys_p = ys_p[order(ys_p$x),]

      # plot data if necessary
      if (plot==TRUE){
        points(ys_p,col='green',pch=4,cex=2)
      }

      # now check for size
      # if there are only 4 we can still use the data
      # discard anything with 3 or less
      if (l_ys > 3){

        # cal_yscul_ysate top l_yseft bottom right, write dedicated routine for this
        bottom_left = ys_p[which(ys_p$y[1:2] == min(ys_p$y[1:2] )),]
        top_left = ys_p[which(ys_p$y[1:2] == max(ys_p$y[1:2] )),]

        top_right = ys_p[which(ys_p$y[c(l_ys-1,l_ys)] == max(ys_p$y[c(l_ys-1,l_ys)]))+(l_ys-2),]
        bottom_right = ys_p[which(ys_p$y[c(l_ys-1,l_ys)] == min(ys_p$y[c(l_ys-1,l_ys)]))+(l_ys-2),]

        # calculate middle points, in case of 5 just use the corners
        # ignore the middle one
        if  (l_ys == 6){
          top_middle = ys_p[which(ys_p$y[c(l_ys-3,l_ys-2)] == max(ys_p$y[c(l_ys-3,l_ys-2)]))+(l_ys-4),]
          bottom_middle = ys_p[which(ys_p$y[c(l_ys-3,l_ys-2)] == min(ys_p$y[c(l_ys-3,l_ys-2)]))+(l_ys-4),]
        }else{
          top_middle = c(NA,NA)
          bottom_middle = c(NA,NA)
        }

        # return data vector
        output[j,1:2] = as.matrix(bottom_left)
        output[j,3:4] = as.matrix(top_left)
        output[j,5:6] = as.matrix(bottom_middle)
        output[j,7:8] = as.matrix(top_middle)
        output[j,9:10] = as.matrix(bottom_right)
        output[j,11:12] = as.matrix(top_right)
      }
    }
  }

  # calculate the median coordinates
  median_coord = apply(output,2,median,na.rm=T)
  pp = matrix(median_coord,6,2,byrow=T)

  # calculate the approximate location of the rows, based upon
  # the corners of the yearly section

  t = row.sections(median_coord)
  t_start = t[[1]]
  t_end = t[[2]]
  d_start = t[[3]]
  d_end = t[[4]]

  if(plot==TRUE){
    # plot the approximate (ideal) row locations
    # intersections
    points(t_start[,1:2],col='red',pch=19)
    points(t_start[,3:4],col='red',pch=19)
    points(t_end[,3:4],col='red',pch=19)
    points(pp[,1],pp[,2],col='red',pch=19)
  }

  # create output matrix containing the 0/1 values
  # of observations during the year
  observations = matrix(0,4,336)

  # fill the coordinate matrix with data
  for (j in 1:l){
    #print(j)
    if (T2[j] == FALSE){
      next
    }
    if (length(subset[j][[1]][[3]]$value)!=0){

      # grab the points marking the yearly section
      ys_p_l = length(subset[j][[1]][[3]]$value)

      line_sections=matrix(NA,ys_p_l,4)
      for (k in 1:ys_p_l){
        if(grepl("preceding",subset[j][[1]][[3]]$value[[k]]$tool_label)){
          #print('coupe')
          next
        }else{
          line_sections[k,1] = subset[j][[1]][[3]]$value[[k]]$x1
          line_sections[k,3] = size_y - subset[j][[1]][[3]]$value[[k]]$y1
          line_sections[k,2] = subset[j][[1]][[3]]$value[[k]]$x2
          line_sections[k,4] = size_y - subset[j][[1]][[3]]$value[[k]]$y2
        }
      }

      for (z in 1:4){
        for (c in c('start','end')){

          if (c == 'start'){
            coord1 = t_start[z,1:2] # check order
            coord2 = t_start[z,3:4]
          }else{
            coord1 = t_end[z,1:2]
            coord2 = t_end[z,3:4]
          }

          # swap line coordinates
          # to they read left to right
          if (coord1[1] > coord2[1]){
            tmp1 = coord1
            tmp2 = coord2
            coord1 = tmp2
            coord2 = tmp1
          }

          for (r in 1:dim(line_sections)[1]){

            # grab the coordinates of the projection of
            # the points defining a line onto the perfect
            # line

            # first coordinate distance
            d1 = point.line.distance(coord1,coord2,line_sections[r,c(1,3)])

            # second coordinate
            d2 = point.line.distance(coord1,coord2,line_sections[r,c(2,4)])

            # if the return is NA skip
            if(is.na(c(d1,d2))){
              #print('bla')
              next
            }

            # if the coordinate is within 1/2 the inter line distance
            # I retain the point (if not, skip)
            if(d1 > (d_start/2) & d2 > (d_start/2)){
              next # kick out the values that are out of range
            }else{

              p1 = point.line.projection(coord1,coord2,line_sections[r,c(1,3)])
              p2 = point.line.projection(coord1,coord2,line_sections[r,c(2,4)])

              # swap the order of the coordinates if not in sequence
              # p1 should be smaller than p2
              if (p1[1] > p2[1]){
                tmp1 = p1
                tmp2 = p2
                p1 = tmp2
                p2 = tmp1
              }

              # skip lines which are outside the first segment
              if (p1[1] > coord2[1] | p2[1] < coord1[1]){
                break
              }

              # limit the range of the markings
              if (p1[1] < coord1[1]){
                p1 = coord1
              }

              if (p2[1] > coord2[1]){
                p2 = coord2
              }

              # calculate the index along the matrix in doy
              td = point.point.distance(coord1,coord2)
              pp1 = point.point.distance(coord1,p1)/td * 167
              pp2 = point.point.distance(coord1,p2)/td * 167

              # correction for second half of the year
              if (c == 'end'){
                pp1 = pp1 + 168
                pp2 = pp2 + 168
              }

              # round data to the nearest (DOY)
              pp1 = round(pp1)
              pp2 = round(pp2)

              if(plot==TRUE){
                # visualize lines
                lines(c(p1[1],p2[1]),c(p1[2],p2[2]),col=colours[z],lwd=2)
              }

              tmp_marks = matrix(0,1,336)
              tmp_marks[pp1:pp2] = 1
              observations[z,] = observations[z,] + tmp_marks

            }
          }
        }
      }
    }
  } # loop over all classifications in a subset

  # set colour scheme for plotting
  colours = c('blue','green','black','purple')


  if (plot == TRUE){
    par(mar=c(5,6,5,4))
    if (max(observations,na.rm=T)==0){
      plot(1,type='n',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
    }else{
      x = observations
      x = cbind(rep(0,4),x,rep(0,4))

      plot(rep(1,337),type='n',xlab='DOY (48 week year)',ylab='',ylim=c(0,50),
           yaxt='n',
           main=sprintf('summary statistics (max. # class. %s out of 10) - %s',max(observations,na.rm=T),image_name))
      abline(v=seq(0,337,28),lwd=2,lty=2,col='grey')
      abline(v=168,lwd=2,lty=2,col='black')
      for (i in 1:4){
        s = x[i,]+(i*10)
        polygon(c(0:337),s,c(337:0),rep(i*10,337),col=colours[i])
      }
      text(y = seq(10, 40, by=10), par("usr")[1],labels = rev(c('flowers','fruit','fr. ground', 'senescence')),
           at = seq(10,40, by=10),srt = 45, pos=2, xpd=TRUE)
    }
  }

  # return the observations
  return(observations)
}
