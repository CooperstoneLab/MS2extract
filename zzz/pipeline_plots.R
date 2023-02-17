# Importing the Spectrum of Procyanidin A2 in negative ionzation mode
# and 20 eV as the collision energy
ProcA2_file <- system.file("extdata",
                        "ProcyanidinA2_neg_20eV.mzXML",
                         package = "MS2extract")

 # Region of interest table (rt in seconds)
 ProcA2_data <- data.frame(Formula = "C30H24O12",Ionization_mode = "Negative",
                      min_rt = 163, max_rt = 180)
 # Importing MS2 data
 ProcA2_raw <- import_mzxml(ProcA2_file, ProcA2_data)

 procA2_extracted <- extract_MS2(ProcA2_raw)

 # Get MS2 TIC
TIC <- get_TIC(ProcA2_raw)
 # Plot MS2 TIC
TIC_results <- plot_tic(TIC)


TIC_plot <- TIC_results$TIC_plot + ggplot2::theme_classic()

ggplot2::ggsave(filename = "TIC_plot.jpeg",  plot = TIC_plot,
                path = "zzz",
                width = 40, height = 30,
                units = "mm", dpi = 600, scale = 2)

ProcA2_ext <- extract_MS2(ProcA2_raw)


ProcA2_detected <- detect_mass(ProcA2_ext$MS2_spec,
                               normalize = TRUE, # Allow normalization
                               min_int = 1) # 1% as minimum intensity



MS2_spectra <- plot_MS2spectra(ProcA2_detected) + ggplot2::theme_classic()

ggplot2::ggsave(filename = "MS2_spectra.jpeg",  plot = MS2_spectra,
                path = "zzz",
                width = 40, height = 30,
                units = "mm", dpi = 600, scale = 2)
