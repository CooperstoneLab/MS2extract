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


#' Imports mzxml files with MS2 scans
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
#' @param file file name of the mzxml file
#' @param ... extra arguments passed to  masstools::read_mzxml()
#' @param roi_table a data frame with two columns min_rt and max_rt specifying
#' the minimum and maximum retention time range (in seconds).
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
#'
#' ProcA2_raw <- import_mzxml(ProcA2_file)
#'
#' # 26731 ions detected in total
#' dim(ProcA2_raw)
#'
#' # Region of interest table (rt in seconds)
#' ROI_dt <- data.frame(min_rt = 163, max_rt = 180)
#' ProcA2_roi <- import_mzxml(ProcA2_file, roi_table = ROI_dt)
#'
#' # 24249 ions detected in ROI
#' dim(ProcA2_roi)
#'
#' }
import_mzxml <- function(file, roi_table = NULL, ...) {
  mzxml_raw <- read_mzxml(file, file, ...)
  scan_info <- extract_scan_info(mzxml_raw) # Scan info

  # Getting scan index number by ROI
 # mzxml_roi <- mzxml_raw[scan_info$Index]

  # Creating a tidy data out of MS2 spectra
  mzxml_tidy <- lapply(mzxml_raw, assign_scan_id ) |>
    dplyr::bind_rows()

  # Keep ions with abundance greater than 0
  mzxml_tidy <- mzxml_tidy |> dplyr::filter(.data$intensity > 0)
  # Eval if roi_table is null
  if(!is.null(roi_table)) mzxml_tidy <- roi_filter(mzxml_tidy, roi_table)

  return(mzxml_tidy)
}
