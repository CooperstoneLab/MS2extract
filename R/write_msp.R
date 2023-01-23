#' Evaluate additional msp fields
#'
#' This internal functions evaluate non-critical msp fields to be included in the
#' msp file. This function will create an empty field if any of the
#' evaluated fields are either NA value in the metadata table or the
#' column is missing.
#'
#' The evaluated additional msp fields are:
#'
#'  \describe{
#'   \item{INCHIKEY}{}
#'   \item{SMILES}{}
#'   \item{CCS}{}
#'   \item{COLLISIONENERGY}{}
#'   \item{INSTRUMENTTYPE}{}
#' }
#' @param msp_attribute a string with the name of the msp field
#' @param spec_metadata a data frame containing the msp field values
#' @param msp_backbone the constructed msp string to export
add_attributes <- function(msp_attribute, spec_metadata, msp_backbone) {

  # Eval attribute in metadata table
  if( msp_attribute %in% names(spec_metadata) ) {

    # Eval if metadtaa column contains data
    if ( !is.na(spec_metadata[,msp_attribute]) ) {
      msp_backbone <- paste0(msp_backbone,
                             msp_attribute, ": ",
                             spec_metadata[,msp_attribute],
                             "\n"
      )
      return(msp_backbone)
      # If does not have metadata, file attribute blank
    } else {
      msp_backbone <- paste0(msp_backbone,
                             msp_attribute, ": ",
                             "\n")
      return(msp_backbone)
    }
    # Write attribute blank if not found in metadata
  } else {
    msp_backbone <- paste0(msp_backbone,
                           msp_attribute, ": ", "\n")
    return(msp_backbone)
  }
}

#' Export MS2 spectra to a msp file
#'
#' Thins functions takes an extracted MS2 spectra and writes in a msp file
#' format. This fuctions takes the extracted MS2 spectra and the metadata
#' for the compound. It will draw information such as the precursor m/z
#' as well as the retention time of the scan from the spec data frame and
#' the rest of information from the provided metadata.
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
#' in the resulting msp file. The following column are required as vital
#' information for the msp output.
#'
#' \describe{
#'  \item{NAME}{}
#'  \item{PRECURSORTYPE}{}
#'  \item{FORMULA}{}
#'  \item{RETENTIONTIME}{}
#'  \item{IONMODE}{}
#' }
#'
#' The follosing fields are also included in the resulting msp files, but are
#' not requiered to be present in the metadata table. If the column does not
#' exist in the column or the value is missing, it will export a blank field.
#'  \describe{
#'   \item{INCHIKEY}{}
#'   \item{SMILES}{}
#'   \item{CCS}{}
#'   \item{COLLISIONENERGY}{}
#'   \item{INSTRUMENTTYPE}{}
#' }
#'
#' @param msp_name a string with the name of the msp file not containing (.msp)
#' @export
#' @examples
#' \dontrun{
#' # Importing the Spectrum of Procyanidin A2 in negative ionization mode
#' # and 20 eV as the collision energy
#' ProcA2_file <- system.file("extdata",
#'                        "ProcyanidinA2_neg_20eV.mzXML",
#'                         package = "MS2extract")
#'
#' ProcA2_raw <- import_mzxml(ProcA2_file)
#'
#' # Detecting masses with the normalized spectra and ions with
#' # intensities greater than 5%
#' ProcA2_detected <- detect_mass(ProcA2_raw, normalize = TRUE, min_int = 1)
#'
#' # Reading the metadata
#' metadata_file <- system.file("extdata",
#'                              "msp_metadata.csv",
#'                               package = "MS2extract")
#' metadata <- reac.csv(metadata_file)
#'
#' # Exporting msp file
#' write_msp(spec = ProcA2_detected,
#'           spec_metadata = metadata,
#'           msp_name = "Procyanidin_A2")
#' }


write_msp <- function(spec = NULL, spec_metadata = NULL, msp_name = NULL) {

  # Number of rows equals the numbers of peaks
  n_peaks <- nrow(spec)

  #if(!dir.exists('spectra')) {
  # dir.create('spectra')
  #}

  # Writing the spectra metadata
  msp_entry <- paste0(
    "NAME: ", spec_metadata$NAME, "\n",
    "PRECURSORMZ: ", unique(spec$mz_precursor),"\n",
    "PRECURSORTYPE: ", spec_metadata$PRECURSORTYPE,"\n",
    "FORMULA: ", spec_metadata$FORMULA,"\n",
    "RETENTIONTIME: ", unique(spec$rt),"\n",
    "IONMODE: ", spec_metadata$IONMODE,"\n",
    "COMMENT: ", "Spectra extracted by MS2extract R package", "\n",
    "Num Peaks: ", n_peaks,"\n"
  )

  # Eval if INCHIKEY is available
  msp_entry <- add_attributes(msp_attribute = "INCHIKEY",
                              spec_metadata = spec_metadata,
                              msp_backbone = msp_entry)

  # Eval if SMILES is available
  msp_entry <- add_attributes(msp_attribute = "SMILES",
                              spec_metadata = spec_metadata,
                              msp_backbone = msp_entry)

  # Eval if CCS is available
  msp_entry <- add_attributes(msp_attribute = "CCS",
                              spec_metadata = spec_metadata,
                              msp_backbone = msp_entry)

  # Eval if COLLISIONENERGY is available
  msp_entry <- add_attributes(msp_attribute = "COLLISIONENERGY",
                              spec_metadata = spec_metadata,
                              msp_backbone = msp_entry)

  # Eval if INSTRUMENTTYPE is available
  msp_entry <- add_attributes(msp_attribute = "INSTRUMENTTYPE",
                              spec_metadata = spec_metadata,
                              msp_backbone = msp_entry)

  #Writing the MS2 spectra on the msp entry
  for (i  in seq(1, n_peaks)) {
    msp_entry <- paste0(msp_entry,
                        round(spec[i, "mz"], 4),
                        " ",
                        spec[i, "intensity"],
                        "\n"
    )
  }

  # Writing in msp file the built string
  sink(paste0(msp_name,".msp"))
  cat(msp_entry)
  sink()
}