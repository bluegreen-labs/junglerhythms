source("R/bounding_box.R")
source("R/geometry_functions.R")
source("R/row_locations.R")
source("R/extract_annotations.R")
source("R/extract_line_sections.R")
library(tidyverse)

df <- extract_annotations(internal = TRUE)

print(df)
