#' Filtering the region of interest (ROI)
#'
#' This functions takes a table with the minimum and maximum retention time,
#' in seconds, to keep the scans encapsulated in the specified time range
#'
#' @param spec a data frame containig
#' @param roi_table a data frame with two columns min_rt and max_rt specifying
#' the minimum and maximum retention time range
#' @export
#' @examples
#' \donttest{
#' # Importing the Spectrum of Procyanidin A2 in negative ionization mode
#' # and 20 eV as the collision energy
#' ProcA2_file <- system.file("extdata",
#'                        "ProcyanidinA2_neg_20eV.mzXML",
#'                         package = "MS2extract")
#'
#' # Region of interest table (rt in seconds)
#' ProcA2_data <- data.frame(Formula = "C30H24O12",Ionization_mode = "Negative",
#'                      min_rt = 163, max_rt = 180)
#'
#' ProcA2_raw <- import_mzxml(ProcA2_file, ProcA2_data)
#'
#' # 26731 ions detected in total
#' dim(ProcA2_raw)
#'
#' # Region of interest table (rt in seconds)
#' ROI_dt <- data.frame(min_rt = 163, max_rt = 180)
#' ProcA2_roi <- import_mzxml(ProcA2_file, roi_table = ROI_dt)
#'
#' # 24249 ions detected in ROI
#' dim(ProcA2_roi)
#'
#' }


roi_filter <- function(spec, roi_table) {

  # Check for the right column names
  if( !identical(names(roi_table), c("min_rt", "max_rt"))  ) {
    stop_message <- paste0(
      "Column names for roi_table must be min_rt and max_rt", "\n",
      names(roi_table), " \n found instead")
    stop(stop_message)
  }

  roi_min <- roi_table$min_rt
  roi_max <- roi_table$max_rt
  spec_roi <- dplyr::filter(spec, .data$rt >= roi_min & .data$rt <= roi_max)
  spec_roi
}
