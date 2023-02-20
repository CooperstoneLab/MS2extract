#' Batch import mzxml files
#'
#' This function imports multiple mzXML files into a named list.
#' It takes a data frame containing the basic information of
#' metabolites such as file name (including the file path), Chemical formula,
#' Ionization mode, as required fields. Additionally, Region Of Interest (ROI)
#' can be provided to narrow down the elution window.
#'
#' @param compounds_dt a data frame containing the following columns.
#' \describe{
#'  \item{Name}{The name of the compound}
#'  \item{Formula}{The compound's chemical formula}
#'  \item{Ionization_mode}{The ionization mode set in data collection
#'  (only Positive and Negative mode allowed).}
#'  \item{File}{The filename of the mzXML file inluding the path}
#' }
#'
#' Additionally, you can provide the ROI by adding two columns
#' in the data frame.
#'
#' \describe{
#'  \item{min_rt}{a double with the minumim retention time to keep}
#'  \item{max_rt}{a double with the minimum retention time to keep}
#' }
#'
#' @return a list with n elements where n is the number of compounds
#' provided in the data.frame
#' @export
#' @examples
#'
#' # Select the csv file name and path
#' batch_file <- system.file("extdata", "batch_read.csv",
#'                           package = "MS2extract")
#' # Read the data frame
#' batch_data <- read.csv(batch_file)
#'
#' # File paths for Procyanidin A2 and Rutin
#' ProcA2_file <- system.file("extdata",
#'                        "ProcyanidinA2_neg_20eV.mzXML",
#'                         package = "MS2extract")
#' Rutin_file <- system.file("extdata",
#'                        "Rutin_neg_20eV.mzXML",
#'                         package = "MS2extract")
#'
#' # Add file path - User should specified the file path -
#' batch_data$File <- c(ProcA2_file, Rutin_file)
#'
#' # Checking batch_data data frame
#' batch_data
#'
#' # Using batch import to import multiple compounds
#' batch_compounds <- batch_import_mzxml(batch_data)
#'
#' # Checking dimension by compound
#' # Procyanidin A2: 24249 ions
#' # Rutin: 22096 ions
#' purrr::map(batch_compounds, dim)
batch_import_mzxml <- function(compounds_dt) {
  # Separating compounds by name
  compounds_list <- split(compounds_dt, f = compounds_dt$Name )

  # Extracting only file names to be read
  compounds_names <- purrr::map(compounds_list, ~.x[,'File'])

  # Extracting min_rt and max_rt from the provided list
  roi_table <- purrr::map(compounds_list,
                          ~ data.frame( Formula = .x["Formula"],
                                        Ionization_mode = .x['Ionization_mode'],
                                        min_rt = .x["min_rt"],
                                        max_rt = .x["max_rt"] ) )

  # Submitting everything to import_mzxml with map
  compounds_out <- purrr::map2(compounds_names, roi_table, # Lists
                               ~ import_mzxml(file = .x,
                                              met_metadata = .y) )

  return(compounds_out)
}
