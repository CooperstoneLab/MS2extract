test_that("Importing mzxml works", {
  # Assign file name
  ProcA2_file <- system.file("extdata",
    "ProcyanidinA2_neg_20eV.mzXML",
    package = "MS2extract"
  )

  # Testing masstools read_mzxml file ----
  mzxml_raw <- MS2extract::read_mzxml(ProcA2_file)
  expect_true(is.list(mzxml_raw))
  expect_true(length(mzxml_raw) == 9)

  # Testing extract scan info ----
  scan_info <- extract_scan_info(mzxml_raw) # Scan info
  expect_true(is.data.frame(scan_info))
  expect_named(scan_info, c("name", "mz", "rt", "file", "Index"))

  # Testing scan_assign_id ----
  mzxml_tidy <- lapply(mzxml_raw, assign_scan_id) |>
    dplyr::bind_rows()
  expect_true(is.data.frame(mzxml_tidy))
  expect_named(mzxml_tidy, c("mz", "intensity", "mz_precursor", "rt"))

  # Testing import_mzxmlx ----
  ProcA2_data <- data.frame(
    Formula = "C30H24O12", Ionization_mode = "Negative",
    min_rt = 163, max_rt = 180
  )
  ProcA2_raw <- import_mzxml(ProcA2_file, ProcA2_data)
  expect_true(is.data.frame(ProcA2_raw))
  expect_named(ProcA2_raw, c("mz", "intensity", "mz_precursor", "rt"))
})
