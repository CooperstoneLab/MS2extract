#' Evaluate additional .msp fields
#'
#' *Internal function*
#'
#' This internal functions evaluate non-critical msp fields to be included in the
#' msp file. This function will create an empty field if any of the
#' evaluated fields are either listed as NA or missing in the metadata table.
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
#' @param msp_attribute a string with the name of an msp field
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
#' *Internal function*
#'
#' This functions takes an extracted MS2 spectra and writes it to a .msp
#' file format. This function incorporates the extracted MS2 spectra along
#' with metadata for the compound.
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
#' @export
#'
write_msp_base <- function(spec = NULL, spec_metadata = NULL) {

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
    "COMMENT: ", "Spectra extracted with MS2extract R package", "\n"
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

  # Adding the number of peaks
  msp_entry <- paste0(
    msp_entry,
    "Num Peaks: ", n_peaks,"\n"
  )

  #Writing the MS2 spectra on the msp entry
  for (i  in seq(1, n_peaks)) {
    msp_entry <- paste0(msp_entry,
                        round(spec[i, "mz"], 5),
                        " ",
                        spec[i, "intensity"],
                        "\n"
    )
  }

  return(msp_entry)
}

#' Export MS2 spectra to a .msp file
#'
#' This functions takes an extracted MS2 spectra and writes it to a .msp
#' file format. This function incorporates the extracted MS2 spectra along
#' with metadata for the compound.
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
#' in the resulting .msp file. The following columns are required as vital
#' information for the .msp output.
#'
#' \describe{
#'  \item{NAME}{}
#'  \item{PRECURSORTYPE}{}
#'  \item{FORMULA}{}
#'  \item{RETENTIONTIME}{}
#'  \item{IONMODE}{}
#' }
#'
#' The following fields can be included in the resulting .msp file, but are
#' not required to be present in the metadata table. If the column does not
#' exist in the column or the value is missing, a blank field will be exported
#'  \describe{
#'   \item{INCHIKEY}{}
#'   \item{SMILES}{}
#'   \item{CCS}{}
#'   \item{COLLISIONENERGY}{}
#'   \item{INSTRUMENTTYPE}{}
#' }
#'
#' @param msp_name a string with the name of the .msp file excluding the file extension
#' @export
#' @examples
#'
#' \dontrun{
#'
#' # Example with single MS2 spectra -----
#'
#' # Importing the Spectrum of Procyanidin A2 in negative ionization mode
#' # and 20 eV as the collision energy
#' ProcA2_file <- system.file("extdata",
#'                        "ProcyanidinA2_neg_20eV.mzXML",
#'                         package = "MS2extract")
#'
#'  # Region of interest table (rt in seconds)
#' ProcA2_data <- data.frame(Formula = "C30H24O12",Ionization_mode = "Negative",
#'                      min_rt = 163, max_rt = 180)
#'
#' ProcA2_raw <- import_mzxml(ProcA2_file, ProcA2_data)
#'
#' # Extracting the most intense MS2 scan
#' ProcA2_ext <- extract_MS2(ProcA2_raw)
#'
#' # Detecting masses with the normalized spectra and ions with
#' # intensities greater than 1%
#' ProcA2_detected <- detect_mass(ProcA2_ext$MS2_spec,
#'                                normalize = TRUE, # Allow normalization
#'                                min_int = 1) # 1% as minimum intensity
#'
#' # Reading the metadata
#' metadata_file <- system.file("extdata",
#'                              "msp_metadata.csv",
#'                               package = "MS2extract")
#' metadata <- read.csv(metadata_file)
#'
#' # Exporting msp file
#' write_msp(spec = ProcA2_detected,
#'           spec_metadata = metadata,
#'           msp_name = "Procyanidin_A2")
#'
#'
#'
#' # Example with batch spectra ----
#'
#'
#' # Select the csv file name and path
#' batch_file <- system.file("extdata", "batch_read.csv",
#'                           package = "MS2extract")
#' # Read the data frame
#' batch_data <- read.csv(batch_file)
#'
#' # File paths for Procyanidin A2 and Rutin
#' ProcA2_file <- system.file("extdata",
#'                        "ProcyanidinA2_neg_20eV.mzXML",
#'                         package = "MS2extract")
#' Rutin_file <- system.file("extdata",
#'                        "Rutin_neg_20eV.mzXML",
#'                         package = "MS2extract")
#'
#' # Add file path - User should specified the file path -
#' batch_data$File <- c(ProcA2_file, Rutin_file)
#'
#' # Checking batch_data data frame
#' batch_data
#'
#' # Reading batch metadata
#' metadata_msp_file <- system.file("extdata",
#'                                 "batch_msp_metadata.csv",
#'                                  package = "MS2extract")
#'
#' metadata_msp <- read.csv(metadata_msp_file)
#'
#' # Using batch import to import multiple compounds
#' batch_compounds <- batch_import_mzxml(batch_data)
#' # Checking dimension by compound
#' # Procyanidin A2: 24249 ions
#' # Rutin: 22096 ions
#' purrr::map(batch_compounds, dim)
#'
#' batch_extracted_compounds <- batch_extract_MS2(batch_compounds)
#'
#' #Batch detect mass
#' batch_mass_detected <- batch_detect_mass(batch_extracted_compounds, # Compound list
#'                                          normalize = TRUE, # Normalize
#'                                          min_int = 1) # Minimum intensity
#'
#' # Checking dimension by compound
#' # Procyanidin A2: 107 ions
#' # Rutin: 12 ions
#' purrr::map(batch_mass_detected, dim)
#'
#' # Writing msp file
#'  write_msp(spec = batch_mass_detected,
#'            spec_metadata = metadata_msp,
#'            msp_name = "ProcA2_Rutin_batch")
#' }


write_msp <- function(spec = NULL, spec_metadata = NULL, msp_name = NULL) {
  # Checking of all args are provided
  if( is.null(spec) ) stop("Please provide a spec data")
  if( is.null(spec_metadata) ) stop("Please provide MS2 msp metadata")
  if( is.null(msp_name) ) stop("Please provide a msp file name")

  # Writing for batch compound extraction
  if ( is.list(spec) & !is.data.frame(spec) ) {

    # Spliting metadata per compond
    spec_metadata <- split(spec_metadata, f = spec_metadata$NAME)
    msp_entry <- purrr::map2(.x = spec, # List MS2 spectra table
                              .y = spec_metadata, # Metadata for MSP entry
                              \(x, y) write_msp_base( as.data.frame(x),
                                                     as.data.frame(y) ) ) # Function

    # Pasting all characters to produce a single msp file
    msp_entry <- paste0(msp_entry, collapse = "\n")
  } else { # Write individual msp
    msp_entry <- write_msp_base(spec = spec,
                                spec_metadata = spec_metadata)
  }


  # Writing in msp file the built string
  sink(paste0(msp_name,".msp"))
  cat(msp_entry)
  sink()

}
