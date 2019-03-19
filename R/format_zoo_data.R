#' Reads in Jungle Rhythms classification
#' output file (CSV) and cleans up the data
#' by dropping early values (trials) etc
#'
#' @param file: path to a zooniverse generated CSV
#' @param date: date before which to discard all data (demo data)
#' @param image_file file with image subset details (filename, dimensions)
#' @param path path where to save the data
#' @keywords data, io, transformation
#' @export

format_zoo_data = function(file = NULL,
                         image_file = "data/jungle_rhythms_subsets.rds",
                         date = "2015-12-10",
                         path = NULL){

  # read subset dimensions
  cat("ingesting ancilllary image data ... \n")
  image_subsets = readRDS(image_file)

  # check if file exists
  # if so read in the data
  if (file.exists(file)){
    cat("reading in zooniverse data ... \n")
    df = read.table(file,
                     header = TRUE,
                     sep = ",",
                     colClasses = "character")
  } else {
    stop("file does not exist, check path!")
  }

  # subset the data based upon a date
  # 13 Dec. was the launch date of the project
  df = df[which(as.Date(df$created_at)>as.Date(date)),]
  df = df[df$workflow_id == 125,] # only use the final workflow

  # convert certain classes from character
  # to integers
  cat("converting data types ... \n")
  df$subject_ids = as.integer(df$subject_ids)
  df$workflow_id = as.integer(df$workflow_id)
  df$user_id = as.integer(df$user_id)
  df$classification_id  = as.integer(df$classification_id)

  # add a column for the image names, use regexpr
  # rather than converting from JSON and subsetting (faster)
  # image names are the names of the yearly segment files
  # not the large sheets, extract thosse as the image id
  # used in the species reference file
  cat("extracting image file names ... \n")

  image_file = unlist(lapply(df$subject_data,
      function(x){

        f = jsonlite::fromJSON(x)[[1]]

        # potential names
        name = c("Filename",
                 "#Filename",
                 "filename",
                 "#filename")

        # final name
        name = name[which(name %in% names(f))]

        return(f[name])
      }))

  cat("extracting image numbers ... \n")
  image_file_tmp = gsub("-","_",image_file)
  image_nr = as.numeric(substr(image_file_tmp,2,8))
  image_row = as.numeric(unlist(lapply(strsplit(image_file_tmp,"_"),"[[",3)))
  image_col = as.numeric(gsub(".jpg",
                              "",
                              unlist(lapply(strsplit(image_file_tmp,"_"),"[[",4))))

  # extract the annotation data from the json string
  # and put in a nested list (data frames)
  cat("converting annotation data ... \n")
  annotations = lapply(df$annotations, function(x)jsonlite::fromJSON(x))

  # put everything into a nested list and save as an RDS file
  data_list = list("subject_id" = df$subject_ids,
                   "user_id" = df$user_id,
                   "classification_id" = df$classification_id,
                   "image_file" = image_file,
                   "image_nr" = image_nr,
                   "image_row" = image_row,
                   "image_col" = image_col,
                   "annotations" = annotations,
                   "image_data" = image_subsets)

  # return the data table
  if (is.null(path)){
    return(data_list)
  } else {
    saveRDS(data_list, sprintf("%s/jungle_rhythms_annotations.rds", path))
  }
}
