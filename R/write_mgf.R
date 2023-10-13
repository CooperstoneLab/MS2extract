#' Creates the batch upload table for GNPS batch upload
#'
#' * Internal Function *
#'
#' This function is designed to create the batch upload table required by GNPS
#' [batch upload](https://ccms-ucsd.github.io/GNPSDocumentation/batchupload)
#' job. This function takes both, `spec` data and the `spec_metadata` to create
#' the bath upload table.
#'
#' @param spec extracted spectra
#' @param spec_metadata spectra metadata
#' @param mgf_filename file name of the .mgf library. This name is created in
#' the calling function.

write_gnps_table <- function(spec, spec_metadata, mgf_filename) {

  cmp_name <- paste(spec_metadata$COMPOUND_NAME, # Compound name
                    spec_metadata$COLLISIONENERGY)

  # Ionized mass
  # Calculating the exact mass and the ionized mass
  exact_mass <- Rdisop::getMolecule(unique(spec$Formula))$exactmass

  if (spec_metadata$IONMODE == "Positive") {
    ionized_mass <- exact_mass + 1.00727
    adduct_type <- "M+H"
    # If possitive mode is false, then assume it is negative
  } else {
    ionized_mass <- exact_mass - 1.00727
    adduct_type <- "M-H"
  }



  gnps_table <- data.frame(FILENAME = mgf_filename,
                           SEQ = "*..*",
                           COMPOUND_NAME = cmp_name,
                           COLLISIONENERGY = spec_metadata$COLLISIONENERGY,
                           MOLECULEMASS = ionized_mass,
                           INSTRUMENT = spec_metadata$INSTRUMENT,
                           IONSOURCE = spec_metadata$IONSOURCE,
                           EXTRACTSCAN = spec_metadata$SCANS,
                           SMILES = spec_metadata$SMILES,
                           INCHI = spec_metadata$INCHI,
                           INCHIAUX = spec_metadata$INCHIAUX,
                           CHARGE = 1, # Only supporting single charged now
                           IONMODE = spec_metadata$IONMODE,
                           PUBMED = spec_metadata$PUBMED,
                           ACQUISITION = spec_metadata$ACQUISITION,
                           EXACTMASS = exact_mass,
                           DATACOLLECTOR = spec_metadata$DATACOLLECTOR,
                           ADDUCT = adduct_type,
                           INTEREST = spec_metadata$INTEREST,
                           LIBQUALITY = spec_metadata$LIBQUALITY,
                           GENUS = spec_metadata$GENUS,
                           SPECIES = spec_metadata$SPECIES,
                           STRAIN = spec_metadata$STRAIN,
                           CASNUMBER = spec_metadata$CASNUMBER,
                           PI = spec_metadata$PI)

  return(gnps_table)
}


#' Writes the GNPS .mgf backbone
#'
#' *Internal Function*
#'
#' This function writes the backbone of the .mgf format compatible with GNPS.
#'
#' @param spec extracted spectra
#' @param spec_metadata spectra metadata
#' @param mgf_filename file name of the .mgf library. This name is created in
#' the calling function.
#'

mgf_GNPS_core <- function(spec, spec_metadata, mgf_filename) {


  # Filter for CE
  spec <- CE_filter(spec = spec, spec_metadata = spec_metadata)

  # MGF top core ----
  mgf_core_top <- paste0(
    "BEGIN IONS", "\n",
    "PEPMASS=", round(as.numeric(unique(spec$mz_precursor)), 5), "\n",
    "CHARGE=1",  "\n",#Only supporting single charged
    "MSLEVEL=2",  "\n",# Only supporting MS2 data
    "FILENAME=", mgf_filename, "\n",
    "SCANS=", spec_metadata$SCANS,"\n",
    "SEQ=", "*..*", "\n",
    "IONMODE=", spec_metadata$IONMODE, "\n",
    "ORGANISM=", spec_metadata$SPECIES, "\n",
    "NAME=", paste(spec_metadata$COMPOUND_NAME,
                   spec_metadata$COLLISIONENERGY), "\n",
    "PI=", spec_metadata$PI, "\n",
    "DATACOLECTOR=", spec_metadata$DATACOLLECTOR, "\n",
    "SMILES=", spec_metadata$SMILES, "\n",
    "INCHI=", spec_metadata$INCHI, "\n",
    "INCHIAUX=", spec_metadata$INCHIAUX, "\n",
    "PUBMED=", spec_metadata$PUBMED, "\n",
    #SUBMITUSER=USER
    "LIBRARYQUALITY=", spec_metadata$LIBQUALITY, "\n"
  )


  # Core bottom ----
  # Writing ions ----
  # Number of rows equals number of peaks
  n_peaks <- nrow(spec)
  for (ion in seq(1, n_peaks)) {
    mgf_core_top <- paste0(
      mgf_core_top,
      round(spec[ion, "mz"], 5), "\t", spec[ion, "intensity"], "\n"
    )
  }
  mgf_entry <- paste0(mgf_core_top, "END IONS\n")

  return(mgf_entry)
}

#' Checking names in the provided GNPS metadata
#'
#' * Internal Function *
#'
#' This function evaluates if the names in the provided table  (`met_metadata`)
#' matches the required column names.
#'
#' @param met_metadata library metadata provided according to the GNPS
#' [template](https://ccms-ucsd.github.io/GNPSDocumentation/batchupload).
#'
#' You can check a working example of this template retrieving the
#' `GNPS_template.xlsx` file by using:
#' `system.file("extdata", "GNPS_template.xlsx", package = "MS2extract")`.


check_gnps_metadata <- function(met_metadata) {

  # Names of the required column names in the GNPS metadata
  MS2extract_GNPS_names <- c("COMPOUND_NAME", "INSTRUMENT", "COLLISIONENERGY",
                             "IONSOURCE", "SMILES", "INCHI", "INCHIAUX",
                             "IONMODE", "PUBMED", "ACQUISITION",
                             "DATACOLLECTOR", "INTEREST", "LIBQUALITY", "GENUS",
                             "SPECIES", "STRAIN", "CASNUMBER", "PI")

  provided_names <- names(met_metadata)

  # Checking if the provided names match the required names
  names_check <- MS2extract_GNPS_names %in% provided_names

  if(all(names_check)) {
    return(TRUE) } else {
      return(MS2extract_GNPS_names[!names_check])
    }

}

#' Create the GNPS mgf backbone file format
#'
#' *Internal function*
#'
#' This function facilitates to create the structure of the GNPS .mgf format.
#' For more information about submitting your spectra to GNPS, please
#' visit this [link](https://ccms-ucsd.github.io/GNPSDocumentation/batchupload).
#'
#' You can find the GNPS template spreadsheet in this
#'  [link](https://ccms-ucsd.github.io/GNPSDocumentation/static/Template.xlsx).

#' @param spec a data frame containing the extracted MS2 spectra, the following
#' colums are required:
#'
#' \describe{
#'  \item{mz_precursor}{precursor ion}
#'  \item{rt}{retention time}
#'  \item{mz}{*m/z*  values}
#'  \item{intensity}{intensity values}
#' }
#'
#' @param spec_metadata a data frame containing the values to be including
#' in the resulting mgf file. In this case, this is the minimum mandatory
#' information to be included.
#'
#' The full explanation about fields and field content can be found in
#' GNPS batch library upload
#' [link](https://ccms-ucsd.github.io/GNPSDocumentation/batchupload).
#'
#' For the rest of fields that are included in the final library file,
#' MS2extract will get this information for the extracted spectra.
#'
#' We highly suggest to check the `gnps_template.xlsx` file described in
#' the function example to check what fields/columns are required in
#' MS2extract.
#'
#'  \describe{
#'   \item{COMPOUND_NAME}{character, Metabolite name, it has to be the same
#'   name used in `met_metadata`, the data frame used to import your MS/MS data}
#'   \item{INSTRUMENT}{character, instrument used for data collection}
#'   \item{COLLISIONENERGY}{character, collision energy used in MS/MS
#'    fragmentation}
#'   \item{IONSOURCE}{Ionization source}
#'   \item{SMILES}{character, SMILES chemical structure of your metabolite}
#'   \item{INCHI}{character, Inchi value for the metabolite}
#'   \item{INCHIAUX}{character, INCHIAUX for the metabolite}
#'   \item{IONMODE}{character, ionization polarity}
#'   \item{PUBMED}{character, PUBMED id where you submitted the MS/MS spectra}
#'   \item{ACQUISITION}{character, Crude, Lysate, Commercial, Isolated, Other}
#'   \item{DATACOLLECTOR}{character, Person who collected the MS/MS data}
#'   \item{INTEREST}{character, interest of the MS/MS data}
#'   \item{LIBQUALITY}{character, library quality.
#'    1 for Gold, 2 for Silver, 3 for Bronze}
#'   \item{GENUS}{character, genus of the organism}
#'   \item{SPECIES}{character, species of the organism}
#'   \item{STRAIN}{character, strain of the organism}
#'   \item{CASNUMBER}{character, CAS number of the metabolite}
#'   \item{PI}{character, principal investigator}
#' }
#'
#' @param mgf_name  file name for the exported `mgf`library. It does not have
#' to contain the file extension  `.mgf`.
#'
#' @return if batch spectra are found, this function writes two files,
#' the `.mgf` library and the required `.tsv` table required by GNPS. If single
#' spectrum is detected, it will only write the `.mgf` library.
#' @export
#' @examples
#' # Example with batch spectra ----
#'
#'
#' # Select the csv file name and path
#' batch_file <- system.file("extdata", "batch_read.csv",
#'   package = "MS2extract"
#' )
#' # Read the data frame
#' batch_data <- read.csv(batch_file)
#'
#' # File paths for Procyanidin A2 and Rutin
#' ProcA2_file <- system.file("extdata",
#'   "ProcyanidinA2_neg_20eV.mzXML",
#'   package = "MS2extract"
#' )
#' Rutin_file <- system.file("extdata",
#'   "Rutin_neg_20eV.mzXML",
#'   package = "MS2extract"
#' )
#'
#' # Add file path - User should specified the file path -
#' batch_data$File <- c(ProcA2_file, Rutin_file)
#'
#' # Checking batch_data data frame
#' batch_data
#'
#' # Using batch import to import multiple compounds
#' batch_compounds <- batch_import_mzxml(batch_data)
#' # Checking dimension by compound
#' # Procyanidin A2: 24249 ions
#' # Rutin: 22096 ions
#' purrr::map(batch_compounds, dim)
#'
#' batch_extracted_compounds <- batch_extract_MS2(batch_compounds)
#'
#' # Batch detect mass
#' batch_mass_detected <- batch_detect_mass(batch_extracted_compounds,
#'   normalize = TRUE, # Normalize
#'   min_int = 1 # Minimum intensity
#' )
#'
#' # Reading metadata from GNPS template
#' template_file <- system.file("extdata", "GNPS_template.xlsx",
#'                  package = "MS2extract")
#'
#' gnps_template <- readxl::read_excel(path = template_file,
#'                 sheet = "batch_example")
#'
#' write_mgf_gnps(spec = batch_mass_detected,
#'                spec_metadata = gnps_template,
#'                mgf_name = "PhenolicsDB")
#'
write_mgf_gnps <- function(spec = NULL, spec_metadata = NULL, mgf_name = NULL) {
  # Checking for arguments ----
  if(is.null(spec))
    cli::cli_abort(c("{.field spec} is empty",
                     "i" = "Please provide a spectra extracted with MS2extract"))
  if(is.null(spec_metadata))
    cli::cli_abort(c("{.field spec_metadata} is empty",
                     "i" = "Please provide spectra metadata"))

  if(is.null(mgf_name))
    cli::cli_abort(c("{.field mgf_name} is empty",
                     "i" = "Please provide a name for the .mgf library"))
  # Checking spec_metada ----
  metadata_check <- check_gnps_metadata(met_metadata = spec_metadata)

  if(!isTRUE(metadata_check))
    cli::cli_abort(
      c("Column names provided do not match with the required column names",
        "i" = "Column{?s} {metadata_check} missing.")
    )

  # Writing MFG ----
  # https://fiehnlab.ucdavis.edu/projects/lipidblast/mgf-files
  # the minimum mgf definition contains precursor mass, charge and m/z abundance

  #creating .mgf file
  # Checking only one Ionizatio mode is provided
  uniq_ion <- unique(spec_metadata$IONMODE)  # Extract unique value
  is_uniq_ion <- uniq_ion %in% "Positive" | uniq_ion %in% "Negative"

  # If true, only one mode and matching character was found
  if(!isTRUE(is_uniq_ion))
    cli::cli_abort(
      c("Multiple or not matching ionization mode was provided",
        "i" = "{.filed IONMODE} must be either {.field Positive} or
        {.field Negative}, {uniq_ion} was provided")
    )


  # Create unique keys
  spec_metadata <- check_specdt_mgf_gnps(spec_metadata = spec_metadata)

  # Checking all MS/MS data present in spec_metadta
  check_all_compounds <- check_MS_data_order(spec = spec,
                                             spec_metadata = spec_metadata)


  mgf_filename_base <- paste(mgf_name, unique(spec_metadata$IONMODE),
                          sep = "_")
  mgf_filename <- paste0(mgf_filename_base, ".mgf")


  if (is.list(spec) & !is.data.frame(spec)) {
    # Splitting metadata per compound
    spec_metadata <- dplyr::mutate(spec_metadata, SCANS = seq(1, dplyr::n() ))

    spec_metadata <- split(spec_metadata,  f = ~spec_metadata$KEY, drop = TRUE)

    mgf_entry <- purrr::map2(
      .x = spec,
      .y = spec_metadata,
      .f = function(x, y){
        mgf_core_top <- mgf_GNPS_core(spec = x,
                                      spec_metadata = y,
                                      mgf_filename = mgf_filename)
      }
    )

    mgf_entry <- paste0(mgf_entry, collapse = "\n")

    gnps_table_out <- purrr::map2_df(
      .x = spec,
      .y = spec_metadata,
      .f = function(x, y){
        gnps_table_out <- write_gnps_table(spec = x,
                                           spec_metadata = y,
                                           mgf_filename = mgf_filename)
      }
    )

    readr::write_tsv(x = gnps_table_out, eol = "\n",
                     file = paste0(mgf_filename_base, ".tsv"), append = F)

  } else { # write individual mgf entry
    mgf_core_top <- mgf_GNPS_core(spec = spec, spec_metadata = spec_metadata,
                                  mgf_filename = mgf_filename)
  }


  # Writing in mgf file
  sink(mgf_filename)
  cat(mgf_entry)
  sink()


}

