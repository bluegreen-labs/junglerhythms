#' Reads in Jungle Rhythms classification
#' output file (CSV) and cleans up the data
#' by dropping early values (trials) etc
#'
#' @param file: path to a zooniverse generated CSV
#' @param date: date before which to discard all data (demo data)
#' @keywords data, io, transformation
#' @export

read_zoo_data = function(file = NULL,
                         date = "2015-12-13",
                         path = NULL){

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
                                  function(x)jsonlite::fromJSON(x)[[1]]$'#Filename'))

  cat("extracting image numbers ... \n")
  image_file = gsub("-","_",image_file)
  image_nr = as.numeric(substr(image_file,2,8))
  image_row = as.numeric(unlist(lapply(strsplit(image_file,"_"),"[[",3)))

  # extract the annotation data from the json string
  # and put in a nested list (data frames)
  cat("converting annotation data ... \n")
  annotations = lapply(df$annotations, jsonlite::fromJSON)

  # put everything into a nested list and save as an RDS file
  data_list = list("subject_id" = df$subject_ids,
       "user_id" = df$user_id,
       "classification_id" = df$classification_id,
       "image_file" = image_file,
       "image_nr" = image_nr,
       "image_row" = image_row,
       "annotations" = annotations)

  # return the data table
  if (is.null(path)){
    return(data_list)
  } else {
    saveRDS(df,sprintf("%s/jungle_rhythms_annotations.rds"))
  }
}
