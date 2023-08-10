## code to prepare `Neg20_6545` dataset goes here

usethis::use_data(Neg20_6545, overwrite = TRUE)


devtools::install_github("CooperstoneLab/MS2extractDB", force = T)
devtools::load_all(".")
data(package = "MS2extractDB", "batchRead_Neg20")
path <- system.file("extdata","QTOF_6545/Neg/20",package = "MS2extractDB")
tmp <- mutate(batchRead_Neg20, File = paste0(path, "/",File) )
batch_compounds <- batch_import_mzxml(tmp)
