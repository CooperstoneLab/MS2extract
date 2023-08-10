## code to prepare `Neg20_6545` dataset goes here

usethis::use_data(Neg20_6545, overwrite = TRUE)


devtools::install_github("CooperstoneLab/MS2extractDB", force = T)
devtools::load_all(".")
data(package = "MS2extractDB", "batchRead_Neg20")
path <- system.file("extdata","QTOF_6545/Neg/20",package = "MS2extractDB")
tmp <- mutate(batchRead_Neg20, File = paste0(path, "/",File) )
batch_compounds <- batch_import_mzxml(tmp)


# Use extract batch extract_MS2
batch_extracted <- batch_extract_MS2(batch_compounds,
                                     verbose = TRUE,
                                     out_list = FALSE )

batch_mass_detected <- batch_detect_mass(batch_extracted, # Compound list
                                         normalize = TRUE, # Normalize
                                         min_int = 1)

data(package = "MS2extractDB", "batchMetadt_neg20")

write_msp(batch_mass_detected, spec_metadata = batchMetadt_neg20,
          msp_name = "phenolics_6545_Neg20")
