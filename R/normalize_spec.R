#' Normalization and low intensity signal filtering
#'
#' This function normalizes raw spectra by converting from raw intensities,
#' to 0-100% relative intensity based on the most intense ion.
#' This function also remove low intensity signals based on minimum cutoff.
#' The default cutoff is 1% and the spectra is normalized.
#'
#' @param spec a data frame produced by `import_mzxml()` function.
#' @param min_int an integer, minimum  intensity. If spectra is normalized,
#' this values should be in a 0-100 range. If the spectra is not
#' not normalized, this value must be ion count.
#' @export
#' @examples
#' # Importing the Spectrum of Procyanidin A2 in negative ionization mode
#' # and 20 eV as the collision energy
#' ProcA2_file <- system.file("extdata",
#'   "ProcyanidinA2_neg_20eV.mzXML",
#'   package = "MS2extract"
#' )
#'
#' # Region of interest table (rt in seconds)
#' ProcA2_data <- data.frame(
#'   Formula = "C30H24O12", Ionization_mode = "Negative",
#'   min_rt = 163, max_rt = 180
#' )
#'
#' ProcA2_raw <- import_mzxml(ProcA2_file, ProcA2_data)
#'
#' # Check raw intensities
#' # Ranges:  1.125 851439.500
#' range(ProcA2_raw$intensity)
#'
#' ProcA2_norm <- normalize_spec(ProcA2_raw)
#' # Ranges: 2 100
#' range(ProcA2_norm$intensity)
normalize_spec <- function(spec, min_int = 1) {
  spec_normalized <- spec |>
    dplyr::group_by(.data$rt) |> # Normalized by the base peak
    dplyr::mutate(intensity = round(.data$intensity /
      max(.data$intensity), 2) * 100) |>
    # Filter intensities > than 1% by default
    dplyr::filter(.data$intensity > min_int) |>
    dplyr::ungroup()

  return(spec_normalized)
}
