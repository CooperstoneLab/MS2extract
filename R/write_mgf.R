
#' Create the mgf backbone file format
#'
#' *Internal function*
#'
#' This functions takes an extracted MS2 spectra and creates the backbone
#' of a typical mgf file based on
#' [MATRIXSCIENCE](http://www.matrixscience.com/help/data_file_help.html)
#' information.
#'
#' @param spec a data frame containing the extracted MS2 spectra, the following
#' colums are required:
#'
#' \describe{
#'  \item{mz_precursor}{}
#'  \item{rt}{}
#'  \item{mz}{}
#'  \item{intensity}{}
#' }
#'
#' @param spec_metadata a data frame containing the values to be including
#' in the resulting mgf file. In this case, this is the minimum mandatory
#' information to be included.
#'
#'  \describe{
#'   \item{INCHIKEY}{}
#'   \item{SMILES}{}
#'   \item{CCS}{}
#'   \item{COLLISIONENERGY}{}
#'   \item{INSTRUMENTTYPE}{}
#' }
write_mgf_gnps <- function(spec = NULL, spec_metadata = NULL, mgf_name = NULL) {

  # https://fiehnlab.ucdavis.edu/projects/lipidblast/mgf-files
  # the minimum mgf definition contains precursor mass, charge and m/z abundance

  #creating .mgf file
  mgf_file_name <- paste0(mgf_name, unique(spec_metadata$IONMODE),
                          collapse = "_")
  mgf_file_name <- paste0(mgf_file_name, ".mgf")


  # MGF top core ----
  mgf_core_top <- paste0(
    "BEGIN IONS", "\n",
    "PEPMASS=", round(as.numeric(unique(spec$mz_precursor)), 5), "\n",
    "CHARGE=1",  "\n",#Only supporting single charged
    "MSLEVEL=2",  "\n",# Only supporting MS2 data
    "FILENAME=", mgf_file_name, "\n",
    "SEQ=", "*..*", "\n",
    "IONMODE=", spec_metadata$IONMODE, "\n",
    "ORGANISM=", spec_metadata$SPECIES, "\n",
    "NAME=", spec_metadata$COMPOUND_NAME, "\n",
    "PI=", spec_metadata$PI, "\n",
    "DATACOLECTOR=", spec_metadata$DATACOLLECTOR, "\n",
    "SMILES=", spec_metadata$SMILES, "\n",
    "INCHI=", spec_metadata$INCHI, "\n",
    "INCHIAUX=", spec_metadata$INCHIAUX, "\n",
    "PUBMED=", spec_metadata$PUBMED, "\n",
    #SUBMITUSER=USER
    "LIBRARYQUALITY=", spec_metadata$LIBQUALITY, "\n"
  )

  # Add more attributes

  # Writing ions ----
  # Number of rows equals number of peaks
  n_peaks <- nrow(spec)
  for (ion in seq(1, n_peaks)) {
    mgf_core_top <- paste0(
      mgf_core_top,
      round(spec[ion, "mz"], 5), "\t", spec[ion, "intensity"], "\n"
    )
  }

  mgf_entry <- paste0(mgf_core_top, "END IONS")
  return(mgf_entry)
}

