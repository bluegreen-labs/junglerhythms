#' orthogonal projection of a point onto a line (2D)
#'
#' @param data data to process
#' @param plot plot stuff
#' @param image_subset list of dimensions for image subsets
#' @param internal logical, save to disk or not
#' @export
#'

phenology_dates = function(
  data = "~/Dropbox/Research_Projects/code_repository/bitbucket/junglerhythms/data/jungle_rhythms_annotations.rds",
  plot = TRUE,
  image_path = "~/Dropbox/Research_Projects/working/congo_phenology/data/yearly_subsets_jpg/",
  output_path  = "/scratch/tmp/",
  internal = FALSE
  ){

  # read in data
  if(file.exists(data)){
    cat("reading in data ... \n")
    # read in the data
    data <- readRDS(data)
  }

  if(length(data)!=9){
    stop("not a jungle rhythms data file")
  }

  # list unique files
  files <- unique(data$image_file)
  cat(sprintf("processing %s files ... \n", length(files)))

  # progress bar settings
  pb <- txtProgressBar(0, length(files), style = 3)
  env <- environment()
  i = 1

  # process all 30K files
  output <- lapply(files,
                   function(file){

    tmp <- get("i", envir = env)
    tmp <- tmp + 1
    assign("i", tmp, envir = env)
    setTxtProgressBar(pb, tmp)

    # which annotations consitute a a subset
    subset_loc <- which(data$image_file == file)

    # subset based upon image (a year to process)
    subset <- data$annotations[subset_loc]

    # ancillary data
    # image nr
    image_nr <- unique(data$image_nr[subset_loc])

    # col
    image_col <- unique(data$image_col[subset_loc])

    # row
    image_row <- unique(data$image_row[subset_loc])

    # image height
    img_y <- data$image_data$height[grep(file,
                                         data$image_data$file_name)]

    # get the median coordinates of all classified
    # bounding values
    median_coord <- bounding_box(subset = subset,
                                 img_y = img_y)

    # extrac observations
    observations <- extract_line_sections(subset = subset,
                                          bb = median_coord,
                                          img_y = img_y)

    if (plot){

      # set colour scheme for plotting
      colours = c('blue','green','black','purple')

      output_file <- paste0(output_path,file)
      #cat(sprintf("saving data to: %s\n", output_file))

      # plot original image
      r <- raster::brick(paste0(image_path,
                                file))

      jpeg(output_file, 710, 950)
      par(mfrow=c(2,1))
      par(mar = c(5,5,5,5))
      raster::plotRGB(r,
                      main = file,
                      axes = TRUE)

      if(!is.null(median_coord)){

        row_coordinates = row_locations(median_coord)
        t_start = row_coordinates[[1]]
        t_end = row_coordinates[[2]]

        # plot the median bounding box location
        pp <- matrix(median_coord, 6, 2, byrow=T)
        points(pp, col = "red" , pch = 20, cex = 3)

        # plot the approximate (ideal) row locations
        # intersections
        points(t_start[, 1:2], col = 'blue', pch = 19)
        points(t_start[, 3:4], col = 'green', pch = 19)
        points(t_end[, 3:4], col = 'blue', pch = 19)
      }

      # plot classified results
      par(mar=c(5,6,5,4))
      if (max(observations, na.rm =TRUE) == 0){
        plot(1,type='n',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
      }else{
        x = observations
        x = cbind(rep(0,4),x,rep(0,4))
        plot(rep(1,337),type='n',xlab='DOY (48 week year)',ylab='',ylim=c(0,50),
             yaxt='n',
             main=sprintf('summary statistics (max. # class. %s out of 10) - %s',
                          max(observations,na.rm=T),file))
        abline(v=seq(0,337,28),lwd=2,lty=2,col='grey')
        abline(v=168,lwd=2,lty=2,col='black')
        for (i in 1:4){
          s = x[i,]+(i*10)
          polygon(c(0:337),s,c(337:0),rep(i*10,337),col=colours[i])
        }
        text(y = seq(10, 40, by=10),
             par("usr")[1],
             labels = rev(c('flowers','fruit','fr. ground', 'senescence')),
             srt = 45, pos=2, xpd=TRUE)
      }

      dev.off()
    }

    # add labels to the file
    labels <- rev(c("flowers","fruit","fruit_drop","senescence"))

    return(data.frame(image_nr,
                      image_col,
                      image_row,
                      observations))

  }) # loop over images

  # close progress bar
  close(pb)

  # concat output data
  output <- do.call("rbind", output)

  if(!internal){
    saveRDS(output,
            paste0(output_path,"jungle_rhythms_raw_converted_annotations.rds"))
  } else {
    return(output)
  }
}

