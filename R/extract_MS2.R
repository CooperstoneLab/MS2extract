#' Calculate the TIC
#'
#' This function compute the MS2 TIC of a spectra grouping by retention time.
#'
#' @param spec a data fram with two columns, mz and intensity
get_TIC <- function(spec) {
  spec_tic <- dplyr::group_by(spec, .data$rt) |>
    dplyr::summarise(TIC = sum(.data$intensity) )
  spec_tic
}

#' Plot MS2 TIC
#'
#' This function plots the MS2 TIC and will mark the most intense scan
#' that will be used to extract the MS2 spectra
#' @param TIC a data frame with the TIC
#' @return a list. The first element is a data frame with the retention time
#' with the highest intensity and the TIC plot as the second list element.

plot_tic <- function(TIC) {
  most_intense <- TIC |>
    dplyr::filter(TIC == max(.data$TIC))

  TIC_plot <- ggplot2::ggplot(data = TIC,
                              aes(x = .data$rt, y = .data$TIC)) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 2) +
    ggplot2::geom_point(data = most_intense, color = "red", size = 3) +
    ggplot2::labs(title = "MS2 TIC plot",
                  subtitle = paste0("MS2 spectra at ","rt: ",
                                    most_intense$rt,
                                    " will be exported"),
                  x = "rt (s)", y = "Intensity") +
    ggplot2::theme_bw()

  return_list <- list(most_intense = most_intense, TIC_plot = TIC_plot)
  return(return_list)
}

#' Plot MS2 spectra
#'
#' This functions plot the resulting MS2 spectra of the most intense scan
#'
#' @param spec a data frame containing the MS2 spectra of the most intense
#' scan
#' @param max_labels an integer indicating the maximum number of labels of
#' the most intense ions
#' @return a ggplot plot
#' @export

plot_MS2spectra <- function(spec, max_labels = 5) {

  if(nrow(spec) > 4){
    spec_nrow <- nrow(spec)
    most_inte_label <- dplyr::arrange(spec, - .data$intensity) |>
      dplyr::mutate(intense = seq(dplyr::n()) ) |>
      dplyr::filter(.data$intense < max_labels )
  } else {
    most_inte_label <- spec
  }

  most_inte_label <- most_inte_label |>
    dplyr::mutate(mz = round(.data$mz, 4))

  ms2_spec <- ggplot2::ggplot(spec,
                              aes(.data$mz, .data$intensity)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::geom_point(data = most_inte_label, color = "red") +
    ggplot2::theme_classic() +
    ggrepel::geom_label_repel(data = most_inte_label,
                              aes(label = .data$mz),
                              box.padding = 0.5,
                              max.overlaps = Inf)


  return(ms2_spec)
}

#' Extrac the most MS2 intense scan
#'
#' This function takes raw MS2 signal and looks for the most intense scan
#' and extracts the MS2 spectra. Additionally, it plots the MS2 TIC
#' chromatogram and color with red the most intense scan.
#'
#' @param spec a data frame with the MS2 spectra
#' @param verbose a boolean indicating if the MS2 TIC chromatogram is displayed
#' @param out_list a boolean expresing if the output is a list containing
#' the MS2 spectra plus the TIC chromatogram (verbose = TRUE), or only
#' the data frame with the MS2 spectra (verbose = FALSE).
#' @param  ... any other arguments passed to detect_mass() function
#' @export
#' @examples
#' # Importing the Spectrum of Procyanidin A2 in negative ionzation mode
#' # and 20 eV as the collision energy
#' ProcA2_file <- system.file("extdata",
#'                        "ProcyanidinA2_neg_20eV.mzXML",
#'                         package = "MS2extract")
#' # Importing MS2 data
#' ProcA2_raw <- import_mzxml(ProcA2_file)
#'
#' # Extracting most intense scan ----
#' # Returning plot + MS2 spectra
#' extract_MS2(ProcA2_raw)
#'
#' # Returning MS2 spectra only
#' extract_MS2(ProcA2_raw, out_list = FALSE)

extract_MS2 <- function(spec, verbose = TRUE, out_list = TRUE, ...) {

  # Get MS2 TIC
  TIC <- get_TIC(spec)
  # Plot MS2 TIC
  TIC_results <- plot_tic(TIC)

  # Calculate the most intense scan
  most_intense <- TIC_results$most_intense

  # Get the rt of the most intense scan
  MS2_spec <- dplyr::filter(spec, .data$rt == most_intense$rt)
  spec_plot <- plot_MS2spectra(MS2_spec)


  plot_results <- ggpubr::ggarrange(TIC_results$TIC_plot, spec_plot, ncol = 1)
  # IF verbose = T print the plot
  if(verbose) print(plot_results)

  # Modify what is returned
  if(out_list){
    return( list(MS2_spec = MS2_spec, TIC_plot = plot_results) )
  } else {
      return(MS2_spec)
    }
}
