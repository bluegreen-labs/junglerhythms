# processing
image_index = "data/phenology_archives_species.csv"

# read in phenology archive image index
index <- read.csv2(image_index,
                   sep = "|",
                   header = TRUE,
                   stringsAsFactors = FALSE)

out <- do.call("rbind", apply(index, 1, function(x){
  x <- as.data.frame(t(as.matrix(x)), stringsAsFactors = FALSE)
  x$l <- NA
  splits <- as.vector(unlist(strsplit(x$image,"\\/")))

  l  <- length(splits)
   if(l == 1){

     return(x)
   } else {

     x <- matrix(rep(x, l), nrow = l, byrow = TRUE)
     colnames(x) <- c(colnames(index),"l")
     x <- as.data.frame(x, stringsAsFactors = FALSE)
     x$image <- splits
     x$l <- 1:l
     x$starting_year[x$l > 1] <- NA
     return(x)
   }
}))

out <- as.matrix(out)

write.table(out, "data/phenology_archives_species_long_format.csv",sep = ",",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
