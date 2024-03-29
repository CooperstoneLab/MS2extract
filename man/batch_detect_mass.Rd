% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batch_detect_mass.R
\name{batch_detect_mass}
\alias{batch_detect_mass}
\title{Detect masses for multiple spectra}
\usage{
batch_detect_mass(batch_spect, normalize = TRUE, min_int = 1)
}
\arguments{
\item{batch_spect}{a list of MS/MS scans imported with \code{\link[=batch_import_mzxml]{batch_import_mzxml()}}
function.}

\item{normalize}{a boolean indicating if the MS/MS spectra is normalized by
the base peak before proceeding to filter out low intensity signals
(normalize  = TRUE), if normalize = FALSE  user has to provide the
minimum ion count.}

\item{min_int}{an integer with the minimum ion intensity}
}
\description{
This function takes multiple spectra imported with the \code{\link[=batch_import_mzxml]{batch_import_mzxml()}}
function. Then, it parses the argument to \code{\link[=detect_mass]{detect_mass()}} function.
Briefly, similarly to  the mass detection step in the MZmine workflow,
\code{\link[=detect_mass]{detect_mass()}} filters out low intensity signals.
Here, you can opt to detect masses by the raw ion intensity (ion counts),
or normalize the spectra  to the most abundant ion, and detect as an
intensity percentage of that ion.
The default is set to include all ions that are at least 1\% of the most
abundant ion.
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

# Extracting MS2 spec
batch_extracted <- batch_extract_MS2(batch_compounds,
  verbose = TRUE,
  out_list = FALSE
)

# Batch detect mass
batch_mass_detected <- batch_detect_mass(batch_extracted, # Compound list
  normalize = TRUE, # Normalize
  min_int = 1
) # Minimum intensity

# Checking dimension by compound
# Procyanidin A2: 38 ions
# Rutin: 4 ions
purrr::map(batch_mass_detected, dim)
}
