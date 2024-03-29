#' Filtering the region of interest (ROI)
#'
#' This function takes a data frame with the minimum and maximum retention time
#' (in seconds), and keeps the *scans* inside the provided boundaries. This
#' filter function aims to keep the scans between the provided ROI and remove
#' the scans outside the ROI.
#'
#' @param spec a data frame containing the MS/MS data.
#' @param roi_table a data frame with two columns min_rt and max_rt specifying
#' the minimum and maximum retention time for a specific metabolite.
#' @export
#' @examples
#' \donttest{
#' # Importing the Spectrum of Procyanidin A2 in negative ionization mode
#' # and 20 eV as the collision energy
#' ProcA2_file <- system.file("extdata",
#'   "ProcyanidinA2_neg_20eV.mzXML",
#'   package = "MS2extract"
#' )
#'
#' # Region of interest table (rt ins seconds)
#' ProcA2_data <- data.frame(
#'   Formula = "C30H24O12", Ionization_mode = "Negative",
#'   min_rt = 163, max_rt = 180
#' )
#'
#' ProcA2_raw <- import_mzxml(ProcA2_file, ProcA2_data)
#'
#' # 26731 ions detected in total
#' dim(ProcA2_raw)
#' }
roi_filter <- function(spec, roi_table) {
  # Check for the right column names
  if (!identical(names(roi_table), c("min_rt", "max_rt"))) {
    stop_message <- paste0(
      "Column names for roi_table must be min_rt and max_rt", "\n",
      names(roi_table), " \n found instead"
    )
    stop(stop_message)
  }

  roi_min <- roi_table$min_rt
  roi_max <- roi_table$max_rt
  spec_roi <- dplyr::filter(spec, .data$rt >= roi_min & .data$rt <= roi_max)
  spec_roi
}
