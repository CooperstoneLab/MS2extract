% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batch_import.R
\name{batch_import_mzxml}
\alias{batch_import_mzxml}
\title{Batch import mzxml files}
\usage{
batch_import_mzxml(compounds_dt)
}
\arguments{
\item{compounds_dt}{a data frame containing the following columns.
\describe{
\item{Name}{The name of the compound}
\item{Formula}{The compound's chemical formula}
\item{Ionization_mode}{The ionization mode set in data collection
(only Positive and Negative mode allowed).}
\item{File}{The filename of the \code{.mzXML}/\code{.mzML} file inluding the path}
\item{COLLISIONENERGY}{Collision energy applied in MS/MS fragmentation}
}

Additionally, you can provide the ROI by adding two columns
in the data frame.

\describe{
\item{min_rt}{a double, minumim retention time to keep}
\item{max_rt}{a double maximum retention time to keep}
}}
}
\value{
a list with \emph{n} elements, where n is the number of compounds
provided in the \code{compounds_dt}
}
\description{
This function imports multiple \code{.mzXML}/\code{.mzML} files into a named list.
It takes a data frame containing the basic information of
metabolites such as file name (including the file path), Chemical formula,
Ionization mode, and Collision Energy, as required fields.
Additionally, Region Of Interest (ROI) can be provided to narrow down
the elution window.
}
\examples{

# Select the csv file name and path
batch_file <- system.file("extdata", "batch_read.csv",
  package = "MS2extract"
)
# Read the data frame
batch_data <- read.csv(batch_file)

# File paths for Procyanidin A2 and Rutin
ProcA2_file <- system.file("extdata",
  "ProcyanidinA2_neg_20eV.mzXML",
  package = "MS2extract"
)
Rutin_file <- system.file("extdata",
  "Rutin_neg_20eV.mzXML",
  package = "MS2extract"
)

# Add file path - User should specified the file path -
batch_data$File <- c(ProcA2_file, Rutin_file)

# Checking batch_data data frame
batch_data

# Using batch import to import multiple compounds
batch_compounds <- batch_import_mzxml(batch_data)

# Checking dimension by compound
# Procyanidin A2: 24249 ions
# Rutin: 22096 ions
purrr::map(batch_compounds, dim)
}
