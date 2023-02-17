#' Extract scan info out of the list
#'
#' Given a list created by masstools::read_mzxml(), this function
#' extract the scan info for all scans and put it in a tidy format
#'
#' @param mzxml a list created by masstools::read_mzxml()
#'
#' @return a data.frame containing
#'  \describe{
#'   \item{name}{the scan name (m/z and rt)}
#'   \item{mz}{the precursor m/z ion}
#'   \item{rt}{retention time}
#'   \item{file name}{file name}
#' }

extract_scan_info <- function(mzxml) {
  scan_id <- lapply(mzxml, function(x) x[[1]] ) |> # Extract info data
    dplyr::bind_rows()  |> # Create a dataset of scann info
    dplyr::mutate(Index = seq( dplyr::n() )) # creating a scan index number
  scan_id
}

#' Extract MS2 spectrum info out of the list
#'
#' Given a list created by masstools::read_mzxml(), this function
#' extract the spectra per scan and put it in a tidy form
#'
#' @param scan_list a list created by masstools::read_mzxml()
#'
#' @return a data.frame containing
#'  \describe{
#'   \item{mz}{ion m/z value}
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
    dplyr::mutate(mz_precursor = scan_id$mz_precursor,
                  rt = scan_id$rt)
  scan_data
}


#' Check compound metadata
#'
#' This function evaluates the information provided by the user whether the
#' right column were provided as well the right ionization modes. Then
#' it calculates the neutral exact mass of the compound by the given
#' formula and depending of the ionization mode, it adds or subtracts the mass
#' of a proton to obtain the ionized mass of the compounds.
#'
#' @param met_metadata a data frame with at least the Formula and the
#' Ionization_mode
#' @return a double with the ionized mass of the compound.

check_metadata <- function(met_metadata) {
  compund_info <- c("Formula", "Ionization_mode")
  ionization_modes <- c("Positive", "Negative")

  # Check for Fomrula and Ionization mode in data frame
  if( all(compund_info %in% names(met_metadata)) ) {

    # Calculating the exact mass and the ionized mass
    exact_mass <- Rdisop::getMolecule(met_metadata$Formula)$exactmass

    # Checking for ionization mode
    if( any(met_metadata$Ionization_mode %in% ionization_modes) ){

      # Check for positive ionization mode
      if(met_metadata$Ionization_mode == "Positive") {
        exact_mass <- exact_mass + 1.00727
        # If possitive mode is false, then assume it is negative
      } else exact_mass <- exact_mass - 1.00727

    } else stop("Only Negative and Positive ionizization modes are accepted")

  } else stop("Fomrula and Ionization_mode colums are required")

  return(exact_mass)

}


#' Imports mzXML files with MS2 scans
#'
#' This function imports MS2 data using the masstools::read_mzxml()
#' function to import the data as a list. Each element in a list *list* represents
#' one scan. Subsequently, each element in this list (scans) contains two
#' sublists that gathers the scan information and the spectra per scan.
#'
#' This function will import all the scans contained in the provided mzxml
#' file, and users can use filters on the resulting data to extract the
#' desired information.
#'
#' @param file file name of the mzXML file
#' @param met_metadata a data frame with the following columns.
#' Required:
#' \describe{
#'  \item{Formula}{A character specifying the metabolite formula}
#'  \item{Ionization_mode}{The ionization mode employed in data collection. It
#'  can be only Positive or Negative}
#' }
#'
#' This two columns are mandatory since the formula is employed with Rdisop
#' package to calculate the exact mass and the ionization mode will dictate
#' if the mass of the a proton is added or subtracted.
#'
#' Additionally, you can provide the minimum and maximum retention times
#' to look for the peak by including the following columns:
#'
#' \describe{
#'  \item{min_rt}{a double with the minimum retention time to keep}
#'  \item{max_rt}{a double with the minimum retention time to keep}
#' }
#'
#' This two columns are highly suggested to be included to narrow down the search
#' window.
#'
#' @param ppm the value of the ppm error. 10 ppm is the default value.
#'
#' @param ... extra arguments passed to  masstools::read_mzxml()
#'
#' @return data.frame in a tidy format for MS2 spectra in a tidy format.
#'  \describe{
#'   \item{mz}{ion m/z value}
#'   \item{intensity}{ion intensity count}
#'   \item{precursor_mz}{precursor ion}
#'   \item{rt}{retention time (in seconds)}
#' }
#'
#' @export
#'
#' @examples
#'
#' \donttest{
#' # Importing the Spectrum of Procyanidin A2 in negative ionization mode
#' # and 20 eV as the collision energy
#' ProcA2_file <- system.file("extdata",
#'                        "ProcyanidinA2_neg_20eV.mzXML",
#'                         package = "MS2extract")
#' # Compound metadata without ROI information
#' ProcA2_data <- data.frame(Formula = "C30H24O12",Ionization_mode = "Negative")
#' ProcA2_raw <- import_mzxml(ProcA2_file, met_metadata = ProcA2_data)
#'
#' # 26731 ions detected in total
#' dim(ProcA2_raw)
#'
#' # Region of interest table (rt in seconds)
#' ProcA2_data <- data.frame(Formula = "C30H24O12",Ionization_mode = "Negative",
#'                      min_rt = 163, max_rt = 180)
#' ProcA2_roi <- import_mzxml(ProcA2_file, met_metadata = ProcA2_data)
#'
#' # 24249 ions detected in ROI
#' dim(ProcA2_roi)
#'
#' }
import_mzxml <- function(file = NULL, met_metadata = NULL, ppm = 10, ...) {

  # Check info in met_metadta ---
  ionized_mass <- check_metadata(met_metadata)

  ppm_error <-  ppm_range(mz = ionized_mass, ppm = ppm)

  mzxml_raw <- read_mzxml(file, ...)
  scan_info <- extract_scan_info(mzxml_raw) # Scan info

  # Getting scan index number by ROI
  # mzxml_roi <- mzxml_raw[scan_info$Index]

  # Creating a tidy data out of MS2 spectra
  mzxml_tidy <- lapply(mzxml_raw, assign_scan_id ) |>
    dplyr::bind_rows()

  # Keep ions with abundance greater than 0
  mzxml_tidy <- mzxml_tidy |> dplyr::filter(.data$intensity > 0)

  # Check for precursor m/z to be in the ppm range ---
  mzxml_tidy <- dplyr::filter(mzxml_tidy,
                              .data$mz_precursor < ppm_error[2] &
                                .data$mz_precursor > ppm_error[1])

  if ( nrow(mzxml_tidy) == 0 ) stop(paste0("Precursor ion not found with the",
                                           "given formula and ppm"))

  # Eval if roi_table is null ----
  # Check if the columns are present
  is_roi_present <- all(c("min_rt", "max_rt") %in% names(met_metadata))

  if(is_roi_present){
    # Creating roi table out of met_metadata
    roi_table  <- dplyr::select(met_metadata, .data$min_rt, .data$max_rt)
    # Filtering using roi table
    mzxml_tidy <- roi_filter(mzxml_tidy, roi_table)
  }

  return(mzxml_tidy)
}
