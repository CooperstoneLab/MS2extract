#' Normalization and low intensity signal filtering
#'
#' This function takes a raw spectra, ion abundance in total counts, and
#' normalized by the base peak. Additionally, this function filters the
#' signal by a minimum intensity to remove background noise. By default,
#' the minimum intensity is 1%.
#'
#' @param spec a data frame produced by import_mzxml function or with the same columns
#' @param min_int a integer the minimum normalized ion intensity between 1 and 100%.
#' @export
#' @examples
#' # Importing the Spectrum of Procyanidin A2 in negative ionzation mode
#' # and 20 eV as the collision energy
#' ProcA2_file <- system.file("extdata",
#'                            "ProcyanidinA2_neg_20eV.mzXML",
#'                             package = "MS2extract")
#'
#' ProcA2_raw <- import_mzxml(ProcA2_file)
#'
#' # Check raw intensities
#' range(ProcA2_raw$intensity) # Ranges:  1.125 851439.500
#'
#' ProcA2_norm <- normalize_spec(ProcA2_raw)
#' range(ProcA2_norm$intensity) # Ranges: 2 100

normalize_spec <- function(spec, min_int = 1) {

  spec_normalized <- spec |> dplyr::group_by(rt) |> # Normalized by the base peak
    dplyr::mutate(intensity = round(intensity / max(intensity), 2 ) * 100 ) |>
    # Filter intensities > than 1% by default
    dplyr::filter(intensity > min_int) |>
    dplyr::ungroup()

  return(spec_normalized)
}
