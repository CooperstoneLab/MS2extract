#' Compare two spectra based on cosine score
#'
#' Wrap function to calculate the cosine score  between two spectra.
#' This function select the *mz* and *intensity* columns before parsing
#' the data frames to OrgMassSpecR::SpectrumSimilarity().
#'
#' @param spec1 a data frame containing the spectra info.
#' @param spec2 a data fram contaning the spectra info.
#' @param output.list a boolean if the the output is returned as a list
#' @param ... arguments parsed to OrgMassSpecR::SpectrumSimilarity().
#' @export
#' @examples
#' # Importing the Spectrum of Procyanidin A2 in negative ionzation mode
#' # and 20 eV as the collision energy
#' ProcA2_file <- system.file("extdata",
#'                        "ProcyanidinA2_neg_20eV.mzXML",
#'                         package = "MS2extract")
#'
#' # Importing the MS2 of Procyanidin A2 deconvoluted by PCDL (Agilent)
#' ProcA2_pcdl_fl <- system.file("extdata",
#'                        "ProcA2_neg_20eV_PCDL.csv",
#'                         package = "MS2extract")
#'
#' # Reading the Procyanidin A2 spectra
#' ProcA2_raw <- import_mzxml(ProcA2_file)
#' ProcA2_extracted <- extract_MS2(ProcA2_raw, out_list = FALSE)
#' ProcA2_norm <- detect_mass(ProcA2_extracted, normalize = TRUE, min_int = 10)
#'
#' # Reading the MS2 spectra of Procynidin A2 by PCDL
#' ProcA2_PCDL <- read.csv(ProcA2_pcdl_fl)
#'
#' compare_spectra(ProcA2_norm, ProcA2_PCDL)

compare_spectra <- function(spec1, spec2, output.list = T, ...) {
  # select mz and intensity columns from spec1
  spec1 <- dplyr::select(.data = spec1, "mz", "intensity")

  # select mz and intensity columns from spec2
  spec2 <- dplyr::select(.data = spec2, "mz", "intensity")


  OrgMassSpecR::SpectrumSimilarity(spec1, spec2, output.list = output.list, ...)

}
