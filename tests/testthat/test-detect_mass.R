test_that("detect_mass works", {
  # Importing the Spectrum of Procyanidin A2 in negative ionization mode
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

  ProcA2_ext <- extract_MS2(ProcA2_raw, verbose = F, out_list = F)


  # Detecting masses with the normalized spectra and ions with
  # intensities greater than 5%
  min_intensity <- 1
  ProcA2_detected <- detect_mass(ProcA2_ext,
    normalize = TRUE,
    min_int = min_intensity
  )


  expect_true(38 == dim(ProcA2_detected)[1])
  expect_equal(max(ProcA2_detected$intensity), 100)

  expect_true(min(ProcA2_detected$intensity) > min_intensity)
})
