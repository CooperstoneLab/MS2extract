#' Compare two spectra based on cosine score
#'
#' A wrapper function to calculate the cosine similarity score between two spectra.
#' This function selects the *m*/*z* and *intensity* columns before parsing
#' the data frames to \code{\link[OrgMassSpecR]{SpectrumSimilarity}}.
#'
#' @param spec1 a data frame containing spectra info.
#' @param spec2 a data frame containing spectra info.
#' @param output.list a boolean, if  TRUE the output is returned as a list
#' @param ... arguments parsed to \code{\link[OrgMassSpecR]{SpectrumSimilarity}}.
#' @export
#' @examples
#' # Importing the Spectrum of Procyanidin A2 in negative ionization mode
#' # and 20 eV as the collision energy
#' ProcA2_file <- system.file("extdata",
#'   "ProcyanidinA2_neg_20eV.mzXML",
#'   package = "MS2extract"
#' )
#'
#' # Importing the MS2 of Procyanidin A2 deconvoluted by PCDL (Agilent)
#' ProcA2_pcdl_fl <- system.file("extdata",
#'   "ProcA2_neg_20eV_PCDL.csv",
#'   package = "MS2extract"
#' )
#' # Region of interest table (rt in seconds)
#' ProcA2_data <- data.frame(
#'   Formula = "C30H24O12", Ionization_mode = "Negative",
#'   min_rt = 163, max_rt = 180
#' )
#' # Reading the Procyanidin A2 spectra
#' ProcA2_raw <- import_mzxml(ProcA2_file, ProcA2_data)
#' # Extracting the most instense MS2 spectra
#' ProcA2_extracted <- extract_MS2(ProcA2_raw, out_list = FALSE)
#' # Detecting masses
#' ProcA2_norm <- detect_mass(ProcA2_extracted, normalize = TRUE, min_int = 1)
#'
#' # Plot of the resulting reference MS2 spectra using MS2extract
#' plot_MS2spectra(ProcA2_norm) +
#'   ggplot2::labs(title = "MS2extract spectra")
#'
#' # Reading the MS2 spectra of Procynidin A2 by PCDL
#' ProcA2_PCDL <- read.csv(ProcA2_pcdl_fl)
#'
#' # Plot of the reference MS2 spectra using PCDL (Agilent software)
#' ggplot2::ggplot(ProcA2_PCDL, ggplot2::aes(mz, intensity)) +
#'   ggplot2::geom_col(width = 1) +
#'   ggplot2::theme_bw()
#'
#' # Cosine comparison between MS2extract and PCDL MS2 spectra
#' compare_spectra(ProcA2_norm, ProcA2_PCDL)
compare_spectra <- function(spec1, spec2, output.list = TRUE, ...) {
  # select mz and intensity columns from spec1
  spec1 <- dplyr::select(spec1, .data$mz, .data$intensity)

  # select mz and intensity columns from spec2
  spec2 <- dplyr::select(spec2, .data$mz, .data$intensity)


  OrgMassSpecR::SpectrumSimilarity(spec1, spec2, output.list = output.list, ...)
}
