#' Base function for MS2 spectra
#'
#' *Internal function*
#'
#' This is the base function for plotting the MS2 spectra.
#'
#' @param spec a data frame containing the MS2 spectra of the most intense
#' scan
#' @param ppm ppm error tolerance to check if the mm/z precursor is being
#' detected or not in the MS2 spectra
plot_MS2base <- function(spec, ppm) {

  # Adding CE: for plot labeling purposes
  spec <- dplyr::mutate(spec, CE = paste0("CE: ", CE))

  # Check if the precursor ion is being detected in MS2 spectra -----
  # Round the precursor ion to 4 digits
  precursor_ion <- unique(spec$mz_precursor)
  # Get the range of ppm given a precursor ion
  precursor_range <- ppm_range(precursor_ion, ppm = ppm)

  # Filter the precursor ion in ions detected in MS2
  precursor_table <- dplyr::filter(
    spec,
    .data$mz > precursor_range[1] &
      .data$mz < precursor_range[2]
  )

  precursor_ion <- round(unique(spec$mz_precursor), 5)

  # Stop if there is more than one precursor ion
  if (nrow(precursor_table) > 1) stop("More than one precursor ion detected")

  ms2_spec <- ggplot2::ggplot(
    spec,
    aes(.data$mz, .data$intensity)
  ) +
    ggplot2::geom_col(width = 1)

  if (nrow(precursor_table) > 0) {
    # Adding 5% of intensity to display diamond
    precursor_table <- dplyr::mutate(precursor_table,
      intensity = .data$intensity * 1.05
    )

    ms2_spec <- ms2_spec +
      ggplot2::geom_point(
        data = precursor_table, shape = 23, size = 2,
        fill = "blue"
      ) +
      ggrepel::geom_label_repel(
        data = precursor_table,
        aes(label = precursor_ion)
      )
  } else {
    repel_data <- data.frame(mz = precursor_ion, intensity = 0)
    ms2_spec <- ms2_spec +
      ggplot2::geom_point(aes(x = precursor_ion, y = 0),
        shape = 23,
        size = 2, fill = "white", color = "blue"
      ) +
      ggrepel::geom_label_repel(
        data = repel_data,
        aes(label = .data$mz)
      )
  }

  if( length(unique(spec$CE)) > 1 ) {ms2_spec <- ms2_spec + facet_wrap("CE")}

  # Changing theme
  ms2_spec <- ms2_spec + theme_light() + # Using a white background theme
    theme(legend.position = "none", # Removing legend
          panel.grid.major = ggplot2::element_blank(),  # Removing grinds
          panel.grid.minor = ggplot2::element_blank())

  return(ms2_spec)
}



#' Plot MS2 spectra
#'
#' This functions plots the resulting MS2 spectra of the most intense scan
#'
#' @param spec a data frame containing the MS2 spectra of the most intense
#' scan
#' @param compound a character, if user is using batch_* functions, they need
#'  to provide a character with the identical “Name” of the compound provided
#'  in the `batch_import_mzxml()` table.
#' @param ppm mass error in ppm tolerance to check if the mm/z precursor is being
#'  detected or not in the MS2 spectra.
#'
#' @return a `ggplot` plot of the MS2 spectra. A filled blue diamond is
#'  placed above the precursor ion. If the precursor ion was not
#'  detected in the MS2 spectra, the blue diamond is not filled.
#'
#' @export
#' @examples
#'
#' Rutin_file <- system.file("extdata",
#'   "Rutin_neg_20eV.mzXML",
#'    package = "MS2extract"
#'  )
#'
#'  # Region of interest table (rt in seconds)
#' Rutin_data <- data.frame(Formula = "C27H30O16",
#'     Ionization_mode = "Negative",
#'     min_rt = 160, max_rt = 175
#'  )
#'  # Importing MS2 data
#' rutin_raw <- import_mzxml(Rutin_file, Rutin_data)
#' Rutin_extracted <- extract_MS2(rutin_raw)
#'
#' Rutin_detected <- detect_mass(Rutin_extracted,
#'    normalize = TRUE, # Allow normalization
#'    min_int = 1) # 1% as minimum intensity
#' MS2_spectra <- plot_MS2spectra(Rutin_detected)
#' print(MS2_spectra)

plot_MS2spectra <- function(spec, compound = NULL, ppm = 10) {
  # Checking if spec is single object or list

  is_list <- is.list(spec)
  is_dtframe <- is.data.frame(spec)

  # If spec is single compound proceed to plot
  if (any(c(!is_list, is_dtframe))) {
    plot_MS2 <- plot_MS2base(spec = spec, ppm = ppm)
  } else {
    # Checking compound is providing
    if (is.null(compound)) stop("Please, provide a compound name")

    plot_MS2 <- plot_MS2base(spec = spec[[compound]], ppm = ppm)
  }

  plot_MS2 <- plot_MS2 + ggplot2::labs(x = "m/z", y = "Intensity")

  return(plot_MS2)
}
