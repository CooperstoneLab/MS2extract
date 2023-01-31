test_that("Compare spectra works", {
  # Importing the Spectrum of Procyanidin A2 in negative ionzation mode
  # and 20 eV as the collision energy
  ProcA2_file <- system.file("extdata",
                             "ProcyanidinA2_neg_20eV.mzXML",
                             package = "MS2extract")

  # Importing the MS2 of Procyanidin A2 deconvoluted by PCDL (Agilent)
  ProcA2_pcdl_fl <- system.file("extdata",
                                "ProcA2_neg_20eV_PCDL.csv",
                                package = "MS2extract")
  ProcA2_data <- data.frame(Formula = "C30H24O12",Ionization_mode = "Negative",
                            min_rt = 163, max_rt = 180)
  # Reading the Procyanidin A2 spectra
  ProcA2_raw <- import_mzxml(ProcA2_file, ProcA2_data)
  # Normalizing total ion counts (Already normalized)
  ProcA2_norm <- normalize_spec(ProcA2_raw)

  # Reading the MS2 spectra of Procynidin A2 by PCDL
  ProcA2_PCDL <- read.csv(ProcA2_pcdl_fl)

  expect_true(is.data.frame(ProcA2_PCDL))
  expect_true(nrow(ProcA2_PCDL) == 38)
  expect_true(ncol(ProcA2_PCDL) == 2)

})
