#' Base function for MS2 spectra
#'
#' This is the base function for plotting the MS2 spectra.
#'
#' @param spec a data frame containing the MS2 spectra of the most intense
#' scan
#' @param ppm ppm error tolerance to check if the mm/z precursor is being
#' detected or not in the MS2 spectra
plot_MS2base <- function(spec, ppm) {

  # Check if the precursor ion is being detected in MS2 spectra -----
  # Round the precursor ion to 4 digits
  precursor_ion <- unique(spec$mz_precursor)
  # Get the range of ppm given a precursor ion
  precursor_range <- ppm_range(precursor_ion, ppm = ppm)

  # Filter the precursor ion in ions detected in MS2
  precursor_table <- dplyr::filter(spec,
                                   .data$mz > precursor_range[1] &
                                     .data$mz < precursor_range[2])

  precursor_ion <- round(unique(spec$mz_precursor), 5)



  # Stop if there is more than one precursor ion
  if( nrow(precursor_table) > 1) stop("More than one precursor ion detected")

  ms2_spec <- ggplot2::ggplot(spec,
                              aes(.data$mz, .data$intensity)) +
    ggplot2::geom_col(width = 1)  +
    ggplot2::theme_bw()

  if (nrow(precursor_table) > 0){
    # Adding 5% of intensity to display diamond
    precursor_table <- dplyr::mutate(precursor_table,
                                     intensity = .data$intensity * 1.05)

    ms2_spec <- ms2_spec +
      ggplot2::geom_point(data = precursor_table, shape = 23, size = 2,
                          fill = "blue") +
      ggrepel::geom_label_repel(data = precursor_table,
                                aes(label = precursor_ion))
  } else {
    repel_data <- data.frame(mz = precursor_ion, intensity = 0)
    ms2_spec <- ms2_spec +
      ggplot2::geom_point(aes(x = precursor_ion, y = 0), shape = 23,
                          size = 2, fill = "white", color = "blue") +
      ggrepel::geom_label_repel(data = repel_data,
                                aes(label = .data$mz))
  }

  return(ms2_spec)

}



#' Plot MS2 spectra
#'
#' This functions plot the resulting MS2 spectra of the most intense scan
#'
#' @param spec a data frame containing the MS2 spectra of the most intense
#' scan
#' @param compound a character, if user is using batch_* functions, user need
#'  to provide a character with the given name of the compound stated in the
#'  batch_import table.
#' @param ppm ppm error tolerance to check if the mm/z precursor is being
#'  detected or not in the MS2 spectra.
#'
#' @return a ggplot plot. Filled blue diamond means that the precursor ion
#' was detected in the MS2 spectra. If the precursor ion was not detected in
#' the MS2 spectra, the diamond is not filled.
#'
#' @export

plot_MS2spectra <- function(spec, compound = NULL, ppm = 10) {

  # Checking if spec is single object or list

  is_list <- is.list(spec)
  is_dtframe <- is.data.frame(spec)

  # If spec is single compound proceed to plot
  if( any(c(!is_list, is_dtframe)) ) {

    plot_MS2 <- plot_MS2base(spec = spec, ppm = ppm)

  } else{

    # Checking compound is providing
    if(is.null(compound)) stop("Please, provide a compound name")

    plot_MS2 <- plot_MS2base(spec = spec[[compound]], ppm = ppm)
  }

  return(plot_MS2)
}
