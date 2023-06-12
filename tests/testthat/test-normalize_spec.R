test_that("Espectra normalization works", {
  # Importing the Spectrum of Procyanidin A2 in negative ionzation mode
  # and 20 eV as the collision energy
  ProcA2_file <- system.file("extdata",
    "ProcyanidinA2_neg_20eV.mzXML",
    package = "MS2extract"
  )
  ProcA2_data <- data.frame(
    Formula = "C30H24O12", Ionization_mode = "Negative",
    min_rt = 163, max_rt = 180
  )
  ProcA2_raw <- import_mzxml(ProcA2_file, ProcA2_data)

  # Check raw intensities
  no_norm_spec <- range(ProcA2_raw$intensity) # Ranges:  1.125 851439.500
  expect_true(length(no_norm_spec) == 2)

  ProcA2_norm <- detect_mass(ProcA2_raw, normalize = TRUE, min_int = 1)
  norm_spec <- range(ProcA2_norm$intensity) # Ranges: 2 100
  expect_true(length(norm_spec) == 2)
  expect_identical(norm_spec[1], 2)
  expect_identical(norm_spec[2], 100)
})
