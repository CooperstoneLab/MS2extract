#' Extract the most MS2 intense scan from batch spectra
#'
#' This function takes multiple spectra imported with the batch_import_mzxml()
#' function. Then, it passes the spectra to a version of extract_MS2().
#'  Briefly, this function takes a series of MS2 spectra, selects the
#'  most intense scan and extracts the MS2 spectra from it.
#' Additionally, it plots the MS2 TIC chromatogram and colors
#' the most intense scan with red circle, and the precursor ion with a blue
#' diamond
#'
#' @param batch_spect a list created by batch_import_mzxml()
#' @param verbose a boolean indicating if the MS2 TIC chromatogram is displayed
#' @param out_list a boolean expressing if the output is a list containing
#' the MS2 spectra plus the TIC chromatogram (out_list = TRUE), or only
#' the data frame with the MS2 spectra (out_lists = FALSE).
#' @export
#' @examples
#'
#' # Select the csv file name and path
#' batch_file <- system.file("extdata", "batch_read.csv",
#'   package = "MS2extract"
#' )
#' # Read the data frame
#' batch_data <- read.csv(batch_file)
#'
#' # File paths for Procyanidin A2 and Rutin
#' ProcA2_file <- system.file("extdata",
#'   "ProcyanidinA2_neg_20eV.mzXML",
#'   package = "MS2extract"
#' )
#' Rutin_file <- system.file("extdata",
#'   "Rutin_neg_20eV.mzXML",
#'   package = "MS2extract"
#' )
#'
#' # Add file path - User should specified the file path -
#' batch_data$File <- c(ProcA2_file, Rutin_file)
#'
#' # Checking batch_data data frame
#' batch_data
#'
#' # Using batch import to import multiple compounds
#' batch_compounds <- batch_import_mzxml(batch_data)
#'
#' # Checking dimension by compound
#' # Procyanidin A2: 24249 ions
#' # Rutin: 22096 ions
#' purrr::map(batch_compounds, dim)
#'
#' # Use extract batch extract_MS2
#' batch_extracted <- batch_extract_MS2(batch_compounds,
#'   verbose = TRUE,
#'   out_list = FALSE
#' )
batch_extract_MS2 <- function(batch_spect, verbose = TRUE, out_list = TRUE) {
  batch_result <- purrr::map(batch_spect,
    .f = extract_MS2,
    verbose = verbose,
    out_list = out_list
  )
  return(batch_result)
}
