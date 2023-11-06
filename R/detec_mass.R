#' Detect masses in a MS/MS spectra
#'
#' Similarly to  the mass detection step in the MZmine workflow,
#' this function filters out low-intensity signals. Here, you can opt 
#' to detect masses by the raw ion intensity (ion counts), or normalize
#' the spectra to the most abundant ion, and detect as an intensity 
#' percentage of that ion. The default is set to include all ions that
#' are at least 1% of the most abundant ion.
#'
#' @param spec a data frame containing the MS/MS spectra.
#' @param normalize a boolean indicating if the MS/MS spectra are normalized by
#' the base peak before proceeding to filter out low-intensity signal
#' (normalize  = TRUE), if normalize = FALSE the user has to provide the
#' minimum ion count.
#' @param min_int an integer with the minimum ion intensity
#' @export
#' @examples
#'
#' # Importing the Spectrum of Procyanidin A2 in negative ionization mode
#' # and 20 eV as the collision energy
#' ProcA2_file <- system.file("extdata",
#'   "ProcyanidinA2_neg_20eV.mzXML",
#'   package = "MS2extract"
#' )
#' # Region of interest table (rt in seconds)
#' ProcA2_data <- data.frame(
#'   Formula = "C30H24O12", Ionization_mode = "Negative",
#'   min_rt = 163, max_rt = 180
#' )
#' ProcA2_raw <- import_mzxml(ProcA2_file, ProcA2_data)
#' # 26731 ion signals
#' dim(ProcA2_raw)
#'
#' # Detecting masses with the normalized spectra and ions with
#' # intensities greater than 5%
#' ProcA2_detected <- detect_mass(ProcA2_raw, normalize = TRUE, min_int = 5)
#' dim(ProcA2_detected)
detect_mass <- function(spec, normalize = TRUE, min_int = 1) {
  if (normalize) {
    spec_result <- spec |> # Normalized by the base peak
      dplyr::mutate(intensity = round(.data$intensity /
        max(.data$intensity), 2) * 100) |>
      dplyr::filter(.data$intensity > min_int) # Filter intensities > than 1% by default
  } else {
    spec_result <- dplyr::filter(spec, .data$intensity > min_int)
  }

  return(spec_result)
}
