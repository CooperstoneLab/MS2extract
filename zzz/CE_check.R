#' Handle multiple CE values in spec_data in spec_metadata.
#'
#' * Internal Function *
#'
#' This function will evaluate if samples that have multiple CE values in
#' spec data and these values matches de values provided in spec_metadata
#'
#' @param spec MS/MS spectra
#' @param spec_metadata  MS/MS spectra metadata containing `COLLISIONENERGY` as
#' one of the provided columns.

multiple_CE_filter <- function(spec, spec_metadata){
  # Extract number of present CE in spectra
  n_CE_all <- sapply(spec, function(x) length(unique(x$CE)) )

  # Samples with more than 1 CE
  mutiple_CE <- data.frame(nCE = n_CE_all[n_CE_all > 1])
  # Data frame of samples with multiple CEs
  mutiple_CE <- mutiple_CE |>
    dplyr::mutate(Sample = rownames(mutiple_CE)) |>
    dplyr::mutate(Sample = stringr::str_replace(Sample, "\\.", "_") ) |>
    dplyr::mutate(Sample = stringr::str_remove(Sample, " eV") ) |>
    tidyr::separate(Sample, sep = "_", remove = FALSE,
                    into = c("NAME", "CE")) |>
    dplyr::mutate(CE = as.numeric(.data$CE))

  # Case 1 - Multiple CE files with diff mzML files and diff spec_metada rows
  # We will assume if multiple CE model
  mutiple_CE_multiple_mzML <- mutiple_CE  |>
    dplyr::group_by(NAME) |> # Group by compound
    dplyr::mutate(n = n())  |> # If compound n > 1, multiple files were pointed out
    filter(n > 1)


  # Check if all CE in multiple mzML and multiple CE are present in spec_metadata
  spec_met_check <- split(mutiple_CE_multiple_mzML, f = ~NAME)
  spec_dt_check_bl <-
    purrr::map_lgl(spec_met_check,
                   \(check, spec_dt) {
                     per_compound <- dplyr::filter(spec_dt,
                                                   .data$NAME %in% check$NAME)

                     allCE <- (readr::parse_number(per_compound$COLLISIONENERGY) %in%
                                 check$CE)
                     all(allCE)
                   }, spec_dt = spec_metadata
    )

  if(!all(spec_dt_check_bl)) {
    cli::cli_abort("")
  }


}
