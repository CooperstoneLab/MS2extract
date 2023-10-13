#' Filter Collision Energy Values
#'
#' * Internal Function *
#'
#' This function evaluate the collision energy value provided in the
#' `spe_metadata` data frame, and matches this value with the collision
#' energy value found in the MS/MS data. This function uses
#' [readr::parse_number()]  to extract the numeric value provided in
#' the `COLLISIONENERGY` column.
#'
#' @param spec MS/MS spectra
#' @param spec_metadata  MS/MS spectra metadata containing `COLLISIONENERGY` as
#' one of the provided columns.

CE_filter <- function(spec, spec_metadata) {
  # Eval if CE column is present
  names_spec_metdt <- names(spec_metadata)
  is_CEpresent <- "COLLISIONENERGY" %in% names_spec_metdt

  # CE not present
  if(!any(is_CEpresent)) {
    cli::cli_abort(c("Collision Energy column is not present",
                     "i" = "{.COLLISIONENERGY} column is required in spec_metadata",
                     "x" = "You provided {names_spec_metdt}"))
  }

  # CE value in metadata
  CE_dt_text <- unique(spec_metadata$COLLISIONENERGY)
  CE_dt_number <- readr::parse_number(CE_dt_text)

  # CE value in spectra
  CE_spec <- unique(spec$CE)

  is_CE_val_present <- any(CE_spec == CE_dt_number)

  if(!is_CE_val_present) {
    cli::cli_abort(c("CE value provided in metadata does not match CE value in data",
                     "i" = "CE value provided: {CE_dt_number}",
                     "x" = "CE value in MS/MS data: {CE_spec}"))
  }

  cli::cli_li("Filtering MS/MS scans for {CE_dt_number} CE")
  # Filter for provided CE
  spec_filtered <- dplyr::filter(spec, CE == CE_dt_number)
  spec_filtered <- dplyr::ungroup(spec_filtered)
  return(spec_filtered)
}

