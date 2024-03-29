% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_mgf.R
\name{check_gnps_metadata}
\alias{check_gnps_metadata}
\title{Checking names in the provided GNPS metadata}
\usage{
check_gnps_metadata(met_metadata)
}
\arguments{
\item{met_metadata}{library metadata provided according to the GNPS
\href{https://ccms-ucsd.github.io/GNPSDocumentation/batchupload}{template}.

You can check a working example of this template retrieving the
\code{GNPS_template.xlsx} file by using:
\code{system.file("extdata", "GNPS_template.xlsx", package = "MS2extract")}.}
}
\description{
\itemize{
\item Internal Function *
}
}
\details{
This function evaluates if the names in the provided table  (\code{met_metadata})
matches the required column names.
}
