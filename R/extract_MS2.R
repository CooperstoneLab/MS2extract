#' Calculate the total ion chromatogram (TIC)
#'
#' *Internal function*
#'
#' This function computes the MS/MS TIC of a given precursor ion.
#'
#' @param spec a data frame with two columns: mz and intensity
get_TIC <- function(spec) {
  spec_tic <- dplyr::group_by(spec, .data$rt, .data$CE) |>
    dplyr::summarise(TIC = sum(.data$intensity)) |>
    dplyr::group_by(CE)
  spec_tic
}

#' Plot MS2 TIC
#'
#' *Internal function*
#'
#' This function plots the MS2 TIC and will mark the most intense scan
#' that will be used to extract the MS/MS spectra
#'
#' @param spec a data frame with the MS/MS spectra
#'
#' @return a list. The first element is a data frame with the retention time
#' with the highest intensity and the EIC plot as the second list element.

plot_tic <- function(spec) {
  # Calculating the TIC (MS/MS EIC)
  TIC <- get_TIC(spec)

  #Getting the precursor m/z to display in plot
  precursorIon <- as.character(round(unique(spec$mz_precursor), 5))

  # Getting max TIC by CE
  most_intense <- TIC |> dplyr::reframe(TIC = max(TIC))

  # Boolean, more than 1 CE found
  n_mostIntense <- nrow(most_intense) > 1L

  most_intense <- TIC %>% dplyr::filter(TIC %in% most_intense$TIC)

  TIC_plot <- ggplot2::ggplot( data = TIC,
                               aes(x = .data$rt, y = .data$TIC) ) +
    ggplot2::geom_line()

  # Coloring most intense scans
  if(n_mostIntense) {
    # Creating all subtitle - based on CE availables
    subtitles <- dplyr::rowwise (most_intense)  |>
      dplyr::mutate(Subtitle = paste("CE:", CE, "at", round(rt, 2),
                                     " will be exported", collapse = " "))
    subtitles <- paste(subtitles$Subtitle, collapse  = "\n")

    TIC_plot <- TIC_plot +
      ggplot2::geom_point(aes(color = factor(CE) )) +
      ggplot2::labs(color = "CE") +
      ggrepel::geom_label_repel(data = most_intense,
                               aes(label = paste0("CE: ", CE) )) +
      ggplot2::labs(
        title = paste0(precursorIon, " MS/MS EIC plot"),
        subtitle = subtitles
      ) +
      ggsci::scale_color_aaas()

  } else {
    TIC_plot <- TIC_plot + ggplot2::geom_point(size = 2) +
      ggplot2::geom_point(data = most_intense,
                          aes(x = .data$rt, y = .data$TIC),
                          size = 3, inherit.aes = F, color = "red") +
      ggplot2::labs(
        title = paste0(precursorIon, " MS/MS ", "@",most_intense$CE," EIC plot"),
        subtitle = paste0(
          "MS/MS spectra at ", "rt: ",
          round(most_intense$rt, 0), " (s)",
          " will be exported"
        ),
        x = "rt (s)", y = "Intensity"
      )
    }

    TIC_plot <- TIC_plot +
    ggplot2::theme_bw()

  return_list <- list(most_intense = most_intense, TIC_plot = TIC_plot)
  return(return_list)
}


#' Extract the most intense MS/MS  scan
#'
#' This function takes a series of MS/MS spectra, selects the
#'  most intense scan and extracts the MS/MS  spectra from it.
#' Additionally, it plots the MS/MS EIC chromatogram and colors
#' the most intense scan with red circle, and the precursor ion with a blue
#' diamond
#'
#' @param spec a data frame with the MS/MS  spectra
#' @param verbose a boolean indicating if the MS/MS EIC chromatogram
#'  is displayed
#' @param out_list a boolean expressing if the output is a list containing
#' the MS/MS spectra plus the EIC chromatogram (verbose = TRUE), or only
#' the data frame with the MS/MS  spectra (verbose = FALSE).

#' @export
#' @examples
#' # Importing the Spectrum of Procyanidin A2 in negative ionzation mode
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
#' # Importing MS/MS data
#' ProcA2_raw <- import_mzxml(ProcA2_file, ProcA2_data)
#'
#' # Extracting most intense scan ----
#' # Returning plot + MS2 spectra
#' extract_MS2(ProcA2_raw)
#'
#' # Returning MS/MS spectra only
#' extract_MS2(ProcA2_raw, out_list = FALSE)
extract_MS2 <- function(spec, verbose = TRUE, out_list = FALSE) {
  # Get MS2 TIC
  TIC <- get_TIC(spec)
  # Plot MS2 TIC
  TIC_results <- plot_tic(spec)

  # Calculate the most intense scan
  most_intense <- TIC_results$most_intense


  # Get the rt of the most intense scan
  MS2_spec <- dplyr::filter(spec, .data$rt %in% most_intense$rt)


  # IF verbose = T print the plot
  if (verbose){
    spec_plot <- plot_MS2spectra(spec = MS2_spec)
    plot_results <- ggpubr::ggarrange(TIC_results$TIC_plot, spec_plot, ncol = 1)

    print(plot_results)
  }

  # Modify what is returned
  if (out_list) {
    return(list(MS2_spec = MS2_spec, TIC_plot = plot_results))
  } else {
    return(MS2_spec)
  }
}
