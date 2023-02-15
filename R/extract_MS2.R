#' Calculate the Total Ion Chromatogram (TIC)
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
#'
#' @param TIC a data frame with the TIC
#'
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


#' Extract the most MS2 intense scan
#'
#' This function takes raw MS2 signal and looks for the most intense scan
#' and proceed to extracts the MS2 spectra. Additionally, it plots the MS2 TIC
#' chromatogram and color in red the most intense scan.
#'
#' @param spec a data frame with the MS2 spectra
#' @param verbose a boolean indicating if the MS2 TIC chromatogram is displayed
#' @param out_list a boolean expressing if the output is a list containing
#' the MS2 spectra plus the TIC chromatogram (verbose = TRUE), or only
#' the data frame with the MS2 spectra (verbose = FALSE).

#' @export
#' @examples
#' # Importing the Spectrum of Procyanidin A2 in negative ionzation mode
#' # and 20 eV as the collision energy
#' ProcA2_file <- system.file("extdata",
#'                        "ProcyanidinA2_neg_20eV.mzXML",
#'                         package = "MS2extract")
#'
#' # Region of interest table (rt in seconds)
#' ProcA2_data <- data.frame(Formula = "C30H24O12",Ionization_mode = "Negative",
#'                      min_rt = 163, max_rt = 180)
#' # Importing MS2 data
#' ProcA2_raw <- import_mzxml(ProcA2_file, ProcA2_data)
#'
#' # Extracting most intense scan ----
#' # Returning plot + MS2 spectra
#' extract_MS2(ProcA2_raw)
#'
#' # Returning MS2 spectra only
#' extract_MS2(ProcA2_raw, out_list = FALSE)

extract_MS2 <- function(spec, verbose = TRUE, out_list = TRUE) {

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
