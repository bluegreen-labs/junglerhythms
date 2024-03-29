% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_annotations.R
\name{extract_annotations}
\alias{extract_annotations}
\title{Extract all phenology dates}
\usage{
extract_annotations(
  data = "./data-raw/jungle_rhythms_annotations.rds",
  plot = FALSE,
  image_path =
    "~/Dropbox/Research_Projects/working/congo_phenology/data/yearly_subsets_jpg/",
  output_path = "/scratch/tmp/",
  internal = TRUE
)
}
\arguments{
\item{data}{raw annotation data to process (pre-processed with
format_zoo_data())}

\item{plot}{plot the output to a particular directory}

\item{image_path}{location of the images used in the project,
for visualization only.}

\item{output_path}{where to store images and generated data}

\item{internal}{logical, save to disk or not}
}
\description{
Wrapper around the extract_line_sections() function, to loop
over all available annotations and extract the relevant line
sections.
}
