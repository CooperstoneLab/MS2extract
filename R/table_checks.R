#' Sort compound table before importing data
#'
#' *Internal function*
#'
#' This function check for the compound data to have unique values (metabolites)
#' per row before importing the MS/MS data. This function performs the
#' following checks:
#'
#' * Create a unique row key: Compound name + Polarity + CE
#' * Check if there are any duplicated values
#' * Sort the compound table by Name and CE value
#'
#' @param compounds_dt  a data frame containing the following columns.
#' \describe{
#'  \item{Name}{The name of the compound}
#'  \item{Formula}{The compound's chemical formula}
#'  \item{Ionization_mode}{The ionization mode set in data collection
#'  (only Positive and Negative mode allowed).}
#'  \item{File}{The filename of the mzXML file inluding the path}
#'  \item{COLLISIONENERGY}{Collision energy applied in MS/MS fragmentation}
#' }
#'
sort_compound_table <- function(compounds_dt = NULL) {

  # Creating unique keey
   compound_check <- compounds_dt |>
    dplyr::mutate(KEY = paste(.data$Name, .data$Ionization_mode,
                               readr::parse_number(.data$COLLISIONENERGY),
                              sep = "_"))
  # Check if unique keys are repeated
  is_KEY_unique <- length(unique(compound_check$KEY)) == nrow(compounds_dt)

  if(!is_KEY_unique) {

    duplicated_KEY <- compound_check |> dplyr::group_by(.data$KEY) |>
      dplyr::add_count() |> dplyr::filter(.data$n > 1)

    n_dup_key <- nrow(duplicated_KEY)

    cli::cli_abort(
      c("Compounds are not unique in terms of Name, Polarity, or CE",
        i = "Repeated compounds that share the same Name, Polarity or CE",
        x = "Found {n_dup_key} key{?s}: {duplicated_KEY$KEY}")
    )
  }

  # Sorting by Name and CE
  compound_check <- compound_check |>
    dplyr::arrange(.data$Name, .data$COLLISIONENERGY)

  return(compound_check)

}

#' Check spec_metadata for NIST format before exporting the MS/MS library
#'
#' *Internal function*
#'
#' This functions aims to create a unique KEY to sort and iterate the compound
#' table in order to match the keys created when MS/MS was imported in
#' `sort_compound_table()`.
#'
#' The three main steps are:
#'
#' * Create a unique row key: Compound name + Polarity + CE
#' * Check if there are any duplicated values
#' * Sort the compound table by Name and CE value

check_specdt_msp <- function(spec_metadata = NULL) {

  # Checking if names are present to create the keys
  spec_names <- names(spec_metadata)
  expected_names <- c("NAME", "COLLISIONENERGY", "IONMODE")

  # All three names present in hte spec_metadata
  is_names_present <- all(expected_names %in% spec_names)

  if(!is_names_present){
    cli::cli_abort(
      c("Names to create unique keys are not present in {.field spec_metadata}",
        i = "{.field NAME COLISSIONENERGY IONMODE} are expected")
      )
  }

  spec_metadata <-  spec_metadata |>
    dplyr::mutate(KEY = paste(.data$NAME, .data$IONMODE,
                              readr::parse_number(.data$COLLISIONENERGY),
                              sep = "_")) |>
    dplyr::arrange(.data$NAME, .data$COLLISIONENERGY)

  return(spec_metadata)

}



#' Check spec_metadata for GNPS `.mgf` format before exporting the MS/MS library
#'
#' *Internal function*
#'
#' This functions aims to create a unique KEY to sort and iterate the compound
#' table in order to match the keys created when MS/MS was imported in
#' `sort_compound_table()`.
#'
#' The three main steps are:
#'
#' * Create a unique row key: Compound name + Polarity + CE
#' * Check if there are any duplicated values
#' * Sort the compound table by Name and CE value

check_specdt_mgf_gnps <- function(spec_metadata = NULL) {

  # Checking if names are present to create the keys
  spec_names <- names(spec_metadata)
  expected_names <- c("COMPOUND_NAME", "COLLISIONENERGY", "IONMODE")

  # All three names present in hte spec_metadata
  is_names_present <- all(expected_names %in% spec_names)

  if(!is_names_present){
    cli::cli_abort(
      c("Names to create unique keys are not present in {.field spec_metadata}",
        i = "{.field NAME COLISSIONENERGY IONMODE} are expected")
    )
  }

  spec_metadata <-  spec_metadata |>
    dplyr::mutate(KEY = paste(.data$COMPOUND_NAME, .data$IONMODE,
                              readr::parse_number(.data$COLLISIONENERGY),
                              sep = "_")) |>
    dplyr::arrange(.data$COMPOUND_NAME, .data$COLLISIONENERGY)

  return(spec_metadata)
}

#' Check if spec MS/MS data and spec_metadata is proper aligned
#'
#' * Internal Function *
#'
#' This function is intended to check the order of the key of the spec
#' data and the spec metadata. If both keys are aligned, the code continues.
#' If the keys are not aligned, it aborts. This function will try to align the
#' keys, as the the keys are ordered back-end and key alignment is expected.
#'
#' Checks in this function include:
#' * Check equal number of MS/MS data and number of observations in spec_metadata
#' * All key names are aligned
#'
#' @param spec a data frame containing the extracted MS2 spectra, the following
#' colums are required:
#'
#' \describe{
#'  \item{mz_precursor}{}
#'  \item{rt}{}
#'  \item{mz}{}
#'  \item{intensity}{}
#' }
#'
#' @param spec_metadata a data frame containing the values to be including
#' in the resulting .msp file. The following columns are required as vital
#' information for the .msp output.
#'
#' \describe{
#'  \item{NAME}{}
#'  \item{PRECURSORTYPE}{}
#'  \item{FORMULA}{}
#'  \item{RETENTIONTIME}{}
#'  \item{IONMODE}{}
#'  \item{COLLISIONENERGY}{}
#' }
#'
check_MS_data_order <- function(spec = NULL, spec_metadata = NULL) {

  # Number of MS/MS spectra
  n_spec <- length(spec)
  # number of observations
  n_obs <- nrow(spec_metadata)

  if( n_spec != n_obs ){
    cli::cli_abort(
      c("Number of MS/MS spectra and spec_metadata does not match",
        i = "{n_spec} MS/MS spectra found but {n_obs} metadata was provided")
    )
  } else {

    spec_names <- names(spec)
    metadata_keys <- spec_metadata$KEY
    is_key_present <- spec_names %in% metadata_keys

    if(all(is_key_present)) {
      return(TRUE)
    } else {
      keys_absent <- spec_names[!is_key_present]
      n_keys_absent <- sum(!is_key_present)
      cli::cli_abort(
        c("{.field spec_metadta} does not match MS/MS spec keys",
          i = "{n_keys_absent} not found in {.field spec_metadata}",
          x = "Missing keys: {keys_absent}")
      )
    }

  }

}








