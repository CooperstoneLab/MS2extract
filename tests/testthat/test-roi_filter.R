test_that("roi_filter works", {
   # Importing the Spectrum of Procyanidin A2 in negative ionization mode
   # and 20 eV as the collision energy
   ProcA2_file <- system.file("extdata",
                          "ProcyanidinA2_neg_20eV.mzXML",
                           package = "MS2extract")

   # Region of interest table (rt in seconds)
   min_rt_test <- 163
   max_rt_test <- 180
   ProcA2_data <- data.frame(Formula = "C30H24O12",Ionization_mode = "Negative",
                             min_rt = 163, max_rt = 180)
   ProcA2_roi <- import_mzxml(ProcA2_file, ProcA2_data)

   # Eval lower rt ROI
   expect_true(min(ProcA2_roi$rt) > min_rt_test)

   # Eval upper rt ROI
   expect_true(max(ProcA2_roi$rt) < max_rt_test)

})
