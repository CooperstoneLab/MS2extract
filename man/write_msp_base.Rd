% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_msp.R
\name{write_msp_base}
\alias{write_msp_base}
\title{Export MS/ms spectra to a \code{.msp} file}
\usage{
write_msp_base(spec = NULL, spec_metadata = NULL)
}
\arguments{
\item{spec}{a data frame containing the extracted MS/MS spectra, the following
columns are required:

\describe{
\item{mz_precursor}{}
\item{rt}{}
\item{mz}{}
\item{intensity}{}
}}

\item{spec_metadata}{a data frame containing the values to be including
in the resulting \code{.msp} file. The following column are required as vital
information for the \code{.msp} output.

\describe{
\item{NAME}{}
\item{PRECURSORTYPE}{}
\item{FORMULA}{}
\item{RETENTIONTIME}{}
\item{IONMODE}{}
}

The following fields are also included in the resulting msp files, but are
not required to be present in the metadata table. If the column does not
exist in the column or the value is missing, it will export a blank field.
\describe{
\item{INCHIKEY}{}
\item{SMILES}{}
\item{CCS}{}
\item{COLLISIONENERGY}{}
\item{INSTRUMENTTYPE}{}
}}
}
\description{
\emph{Internal function}
}
\details{
This functions takes an extracted MS/MS spectra and writes it to a \code{.msp}
file format. This function incorporates the extracted MS/ms spectra along
with metadata for the compound.
}
