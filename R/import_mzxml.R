#' Extract scan info out of the list
#'
#' *Internal function*
#'
#' Given a list created by `MS2extract:::read_mzxml()`, this function
#' extracts the scan info for all scans in a tidy format.
#'
#' @param mzxml a list created by masstools::read_mzxml()
#'
#' @return a data.frame containing
#'  \describe{
#'   \item{name}{the scan name (combined *m*/*z* and rt)}
#'   \item{mz}{the precursor m/z ion}
#'   \item{rt}{retention time}
#'   \item{CE}{collision energy}
#'   \item{file name}{file name}
#' }

extract_scan_info <- function(mzxml) {
  scan_id <- lapply(mzxml, function(x) x[[1]]) |> # Extract info data
    dplyr::bind_rows() |> # Create a dataset of scan info
    dplyr::mutate(Index = seq(dplyr::n())) # creating a scan index number
  scan_id
}

#' Extract MS/MS spectrum info out of the list
#'
#' *Internal function*
#'
#' This function extracts spectra per scan in a tidy format given a list
#' created by `MS2extract:::read_mzxml()`
#'
#' @param scan_list a list created by `MS2extract::read_mzxml()`
#'
#' @return a data.frame containing
#'  \describe{
#'   \item{mz}{ion *m*/*z* value}
#'   \item{intensity}{ion intensity count}
#'   \item{precursor_mz}{precursor ion}
#'   \item{rt}{retention time}
#' }

assign_scan_id <- function(scan_list) {
  scan_id <- scan_list[[1]] # Extracs scan info
  # rename mz for precursor mz
  scan_id <- dplyr::rename(.data = scan_id, mz_precursor = "mz")
  scan_data <- scan_list[[2]] # Extract scan data
  scan_data <- scan_data |> # Add rt data to the scan
    dplyr::mutate(
      mz_precursor = scan_id$mz_precursor,
      rt = scan_id$rt,
      CE = scan_id$CE
    )
  scan_data
}


#' Check compound metadata
#'
#' *Internal function*
#'
#' This function evaluates if the metadata has the expected column names.
#' Then, it calculates
#' the theoretical exact mass of the compound using the given formula and
#' ionization mode to obtain a charged mass.
#'
#' Required values
#'
#' \describe{
#' \item{Formula}{character, compound chemical formula}
#' \item{Ionization_mode}{character,
#'  only *Positive* and *Negative* values are accepted }
#' }
#'
#' @param met_metadata a data frame with at least the Formula and the
#' Ionization_mode
#' @return a double with the ionized mass of the compound.

check_metadata <- function(met_metadata) {
  compund_info <- c("Formula", "Ionization_mode")
  ionization_modes <- c("Positive", "Negative")

  # Check for Fomrula and Ionization mode in data frame
  if (all(compund_info %in% names(met_metadata))) {
    # Calculating the exact mass and the ionized mass
    exact_mass <- Rdisop::getMolecule(met_metadata$Formula)$exactmass

    # Checking for ionization mode
    if (any(met_metadata$Ionization_mode %in% ionization_modes)) {
      # Check for positive ionization mode
      if (met_metadata$Ionization_mode == "Positive") {
        exact_mass <- exact_mass + 1.00727
        # If possitive mode is false, then assume it is negative
      } else {
        exact_mass <- exact_mass - 1.00727
      }
    } else {
      cli::cli_abort(
        c("Value provided in {.Ionization_mode} is prohibited",
          "i" = "{.Ionization_mode} accepts only Positive or Negative values" )
      )
    }
  } else {
    cli::cli_abort(
      c("{.Fomrula} and {.Ionization_mode} colums are required")
    )

  }

  return(exact_mass)
}


#' Imports mzXML/mzML files with MS/MS scans
#'
#' This function reads `.mzXML` and `.mzML` files containing MS/MS.
#' This function is inspired on `masstools::read_mzxml()` which imports
#' the data as a list. Then, each element in a list represents one scan.
#' Also, each scan contains
#' two sub-lists that contain (1) the scan information and
#' (2) the spectra per scan.
#'
#' @param file file name and path of the .mzXML or .mzML MS/MS data
#' @param met_metadata a data frame with the following columns.
#'
#' **Required**:
#' \describe{
#'  \item{Formula}{A character string specifying the metabolite formula}
#'  \item{Ionization_mode}{The ionization mode employed in data collection. It
#'  can be only Positive or Negative}
#'  \item{Ionization_mode}{The ionization mode set in data collection
#'  (only Positive and Negative modes allowed).}
#'  \item{File}{The filename of the mzXML file including the path}
#'  \item{COLLISIONENERGY}{Collision energy applied in MS/MS fragmentation}
#' }
#'
#' These two columns are mandatory since the formula is used by the Rdisop
#' package to calculate the exact mass and the ionization mode
#' determine whether the mass of a proton is added or subtracted to
#' calculate the charged mass.
#'
#' Additionally, you can provide the minimum and maximum retention times
#' to look for a peak only within a given area by including the following
#' columns:
#'
#' \describe{
#'  \item{min_rt}{numeric, with the minimum retention time to keep}
#'  \item{max_rt}{numeric, with the minimum retention time to keep}
#'
#' }
#'
#' These two columns are highly recommended to be included to narrow down
#' the search window and ensure the peak you want is selected. This is
#' especially important if you have multiple peaks in the same data file.
#'
#' @param ppm the mass error in ppm. 10 ppm is the default value.
#'
#'
#' @return data.frame in a tidy format for MS/MS spectra in a tidy format.
#'  \describe{
#'   \item{mz}{ion *m*/*z* value}
#'   \item{intensity}{ion intensity count}
#'   \item{mz_precursor}{precursor ion}
#'   \item{rt}{retention time (in seconds)}
#'   \item{Formula}{Chemical formula provided in `met_metada`}
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Importing the Spectrum of Procyanidin A2 in negative ionization mode
#' # and 20 eV as the collision energy
#' ProcA2_file <- system.file("extdata",
#'   "ProcyanidinA2_neg_20eV.mzXML",
#'   package = "MS2extract"
#' )
#' # Compound metadata without ROI information
#' ProcA2_data <- data.frame(Formula = "C30H24O12", Ionization_mode = "Negative")
#' ProcA2_raw <- import_mzxml(ProcA2_file, met_metadata = ProcA2_data)
#'
#' # 26731 ions detected in total
#' dim(ProcA2_raw)
#'
#' # Region of interest table (rt in seconds)
#' ProcA2_data <- data.frame(
#'   Formula = "C30H24O12", Ionization_mode = "Negative",
#'   min_rt = 163, max_rt = 180
#' )
#' ProcA2_roi <- import_mzxml(ProcA2_file, met_metadata = ProcA2_data)
#'
#' # 24249 ions detected in ROI
#' dim(ProcA2_roi)
#' }
import_mzxml <- function(file = NULL, met_metadata = NULL, ppm = 10) {
  # Check info in met_metadta ---
  ionized_mass <- check_metadata(met_metadata)

  ppm_error <- ppm_range(mz = ionized_mass, ppm = ppm)

  mzxml_raw <- read_mzxml(file)
  scan_info <- extract_scan_info(mzxml_raw) # Scan info

  # Getting scan index number by ROI
  # mzxml_roi <- mzxml_raw[scan_info$Index]

  # Creating a tidy data out of MS2 spectra
  mzxml_tidy <- lapply(mzxml_raw, assign_scan_id) |>
    dplyr::bind_rows()

  # Keep ions with abundance greater than 0
  mzxml_tidy <- mzxml_tidy |> dplyr::filter(.data$intensity > 0)

  # Check for precursor m/z to be in the ppm range ---
  mzxml_tidy <- dplyr::filter(
    mzxml_tidy,
    .data$mz_precursor < ppm_error[2] &
      .data$mz_precursor > ppm_error[1]
  )
  #ppm_calculated <- unique(mzxml_tidy$mz_precursor)
  cli::cli_li("m/z range given {ppm} ppm: {round(ppm_error, 5)}")
  if (nrow(mzxml_tidy) == 0) {
    #metabolite <- met_metadata$Name
    formula <- met_metadata$Formula
    file <- basename(file)
    ppm_error <- round(ppm_error, 4)
    cli::cli_abort(
      c("Precursor ion not found in {file} ",
        "i" = "given {formula} and {ppm} ppm: {ppm_error} m/z range was evaluated" )
    )
  }

  # Eval if roi_table is null ----
  # Check if the columns are present
  is_roi_present <- all(c("min_rt", "max_rt") %in% names(met_metadata))

  if (is_roi_present) {
    # Creating roi table out of met_metadata
    roi_table <- dplyr::select(met_metadata, min_rt, max_rt)
    # Filtering using roi table
    mzxml_tidy <- roi_filter(mzxml_tidy, roi_table)
  }

  mzxml_tidy <- mutate(mzxml_tidy, Formula = met_metadata$Formula)
  mzxml_tidy <- dplyr::group_by(mzxml_tidy, Formula, CE)
  cli::col_green("Succesfully imported!")
  return(mzxml_tidy)
}
