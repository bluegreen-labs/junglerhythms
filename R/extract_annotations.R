#' Extract all phenology dates
#'
#' Wrapper around the extract_line_sections() function, to loop
#' over all available annotations and extract the relevant line
#' sections.
#'
#' @param data raw annotation data to process (pre-processed with
#' format_zoo_data())
#' @param plot plot the output to a particular directory
#' @param image_path location of the images used in the project,
#' for visualization only.
#' @param internal logical, save to disk or not
#' @param output_path where to store images and generated data
#' @export
#'

extract_annotations = function(
  data = "./data-raw/jungle_rhythms_annotations.rds",
  plot = FALSE,
  image_path = "~/Dropbox/Research_Projects/working/congo_phenology/data/yearly_subsets_jpg/",
  output_path  = "/scratch/tmp/",
  internal = TRUE
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

  files <- "R1066197_c_4_9.jpg"

  # process all +30K files
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
    observations <- line_sections(subset = subset,
                                          bb = median_coord,
                                          img_y = img_y)

    print(observations[3,])

    if(!is.null(observations)){
    # Should the conversions be plotted for review?
    if (plot){
      panel_plot(observations = observations,
                       median_coord = median_coord,
                       file = file,
                       image_path = image_path,
                       output_path = output_path)
    }

    # bind meta-data
    return(data.frame(image_nr,
                      image_col,
                      image_row,
                      observations))
    }
  })

  # close progress bar
  close(pb)

  # concat output data
  output <- do.call("rbind", output)

  if(!internal){
    saveRDS(output,
            paste0(output_path,"jungle_rhythms_daily_annotations.rds"))
  } else {
    return(output)
  }
}


#' Panel plot of annotations
#'
#' Visualizes the extracted annotations (for internal use only)
#'
#' @param observations annotations for a particular yearly section
#' @param median_coord corner point coordinates
#' @param file yearly section (unique filename) to process
#' @param image_path location of the images used in the project,
#' for visualization only.
#' @param output_path where to store images and generated data
#' @export
#'

panel_plot <- function(observations,
                             median_coord,
                             file,
                             image_path,
                             output_path){

  # plot original image
  r <- raster::brick(paste0(image_path,
                            file))

  # jpeg(file.path(output_path, file),
  #      710,
  #      950)
  par(mfrow=c(2,1))
  par(mar = c(5,5,5,5))
  raster::plotRGB(r,
                  main = file,
                  axes = TRUE)

  if(!is.null(median_coord)){
    row_coordinates <- row_locations(median_coord)
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
  par(mar=c(5,6,5,4))
  plot_raw_annotations(observations)

  # dev.off()
}


