#' Detect masses for multiple spectra
#'
#' This function takes multiple spectra imported with batch_import_mzxml()
#' function. Then, it parses the argument to detect_mass() function.
#'
#' Briefly, similarly to MZmine, detect mass filters out low intensity signals
#' that can be attributed to background noise. Here, you can opt to detect masses
#' by the raw ion intensity or normalize the spectra by the base peak ion
#' and detect masses by intensity percentage from 1% (by default) to 100%.
#'
#' @param batch_spect a list of MS2 scans imported with batch_import_mzxml()
#'  function.
#' @param normalize a boolean indicating if the MS2 spectra is normalized by
#' the base peak before proceeding to filter out low intensity signals
#' (normalize  = TRUE), if normalize = FALSE the user has to provide the
#' minimum ion count.
#' @param min_int a integer with the minimum ion intensity
#' @export
#' @examples
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
#' # Using batch import to import multiple compounds
#' batch_compounds <- batch_import_mzxml(batch_data)
#' # Checking dimension by compound
#' # Procyanidin A2: 24249 ions
#' # Rutin: 22096 ions
#' purrr::map(batch_compounds, dim)
#'
#' # Extracting MS2 spec
#' batch_extracted <- batch_extract_MS2(batch_compounds,
#'                                      verbose = TRUE,
#'                                      out_list = FALSE)
#'
#' #Batch detect mass
#' batch_mass_detected <- batch_detect_mass(batch_extracted, # Compound list
#'                                          normalize = TRUE, # Normalize
#'                                          min_int = 1) # Minimum intensity
#'
#' # Checking dimension by compound
#' # Procyanidin A2: 38 ions
#' # Rutin: 4 ions
#' purrr::map(batch_mass_detected, dim)

batch_detect_mass <- function(batch_spect, normalize = TRUE, min_int = 1) {

  # Checking if batch spec is an embedded list with MS2 TIC plot inside
  is_dataframe_valid <- batch_spect |>
    purrr::map_lgl(is.data.frame) |>
    all()

  if( !is_dataframe_valid ) {

    # Extracting only the MS2_spec data frame
    batch_spect <- purrr::map(batch_spect, \(x) x[["MS2_spec"]])

  }

  # Parsing arguments to detect_mass()
  batch_result <- purrr::map(batch_spect,
                             .f = detect_mass,
                             normalize = normalize,
                             min_int = min_int )
  return(batch_result)
}
