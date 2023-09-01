

#' @title Read mzxml/mzML files containing MS2 data
#' @description This function takes a mzxml file and imports the MS2 data only
#'  to further create a list contaning to sublists with the scan information
#'  (fisrt sublits) and the MS2 spectra (second sublist)
#' @author Xiaotao Shen inspired from the masstools package https://www.tidymass.org
#' @param file The vector of names of ms2 files. MS2 file must be mzXML or mzML.
#' @param threads Thread number
#' @param mode inMemory or onDisk
#' @return Return ms2 data. This is a list.
#' @export
#' @examples
#' ProcA2_file <- system.file("extdata",
#'   "ProcyanidinA2_neg_20eV.mzXML",
#'   package = "MS2extract"
#' )
#'
#' ProcA2_raw <- read_mzxml(ProcA2_file)
read_mzxml <- function(file,
                       threads = 3,
                       mode = c("inMemory", "onDisk")) {
  #message(crayon::green( paste0("Reading MS2 data from ") ))

  ms2 <- MSnbase::readMSData(
    files = file,
    msLevel. = 2,
    mode = mode
  )
  ms2_CE <- MSnbase::collisionEnergy(ms2)

  #message(crayon::green(paste("Processing:",  base::basename(file)) ))

  new.ms2 <- ProtGenerics::spectra(object = ms2)
  unique_CE <- unique(ms2_CE)
  n_CE <- length(unique_CE)
  cli::cli_li(c(
    paste("Processing: ",  base::basename(file)),
    "Found {n_CE} CE value{?s}: {unique_CE}",
    "Remember to match CE velues  in {.emph spec_metadata}
    when exporting your library"))

  #cli::cli_text("...")


  rm(list = c("ms2"))

  new.ms2 <- seq_along(new.ms2) %>%
    purrr::map(function(idx) {
      temp.ms2 <- new.ms2[[idx]]
      info <-
        data.frame(
          name = paste("mz", temp.ms2@precursorMz,
            "rt", temp.ms2@rt,
            sep = ""
          ),
          "mz" = temp.ms2@precursorMz,
          "rt" = temp.ms2@rt,
          "CE" = ms2_CE[idx],
          "file" = basename(file[temp.ms2@fromFile]),
          stringsAsFactors = FALSE
        )
      duplicated.name <-
        unique(info$name[duplicated(info$name)])
      if (length(duplicated.name) > 0) {
        lapply(duplicated.name, function(x) {
          info$name[which(info$name == x)] <-
            paste(x, seq_len(sum(info$name == x)), sep = "_")
        })
      }

      rownames(info) <- NULL
      spec <- data.frame(
        "mz" = temp.ms2@mz,
        "intensity" = temp.ms2@intensity,
        stringsAsFactors = FALSE
      )
      list(info = info, spec = spec)
    })


  new.ms2 <- new.ms2
  new.ms2
}
