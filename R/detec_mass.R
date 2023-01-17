#' Detect masses in a spectra
#'
#' Similarlly to MZmine, detect mass filters out low intensity signal that can
#' be attributed to background noise. Here, you can opt to detect masses
#' by the raw ion intensity or normalize and detect masses by intensity
#' percentabe from 1% (by default) to 100%.
#'
#' @param spec a data frame containing the MS2 spectra
#' @param normalize a boolean indicating if the MS2 spectra is normalized by
#' the base peak before proceeding to filter out low intensity signal
#' (normalize  = TRUE), if normalize = FALSE the user has to provide the
#' minimum ion count.
#' @param min_int a integer with the minimum ion intensity
#' @export
#' @examples
#'
#' # Importing the Spectrum of Procyanidin A2 in negative ionization mode
#' # and 20 eV as the collision energy
#' ProcA2_file <- system.file("extdata",
#'                        "ProcyanidinA2_neg_20eV.mzXML",
#'                         package = "MS2extract")
#'
#' ProcA2_raw <- import_mzxml(ProcA2_file)
#' # 26731 ion signals
#' dim(ProcA2_raw)
#'
#' # Detecting masses with the normalized spectra and ions with
#' # intensities greater than 5%
#' ProcA2_detected <- detect_mass(ProcA2_raw, normalize = TRUE, min_int = 5)
#' dim(ProcA2_detected)
detect_mass <- function(spec, normalize = TRUE, min_int = 1){

  if(normalize){
    spec_result <- spec |> # Normalized by the base peak
      dplyr::mutate(intensity = round(.data$intensity /
                                        max(.data$intensity), 2 ) * 100 ) |>
      dplyr::filter(intensity > min_int) # Filter intensities > than 1% by default
  } else {
    spec_result <- dplyr::filter(spec, "intensity" > min_int)
  }

  return(spec_result)
}
