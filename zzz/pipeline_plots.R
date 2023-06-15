# Procyanidin A2 ----

# Importing the Spectrum of Procyanidin A2 in negative ionization mode
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

ggplot2::ggsave(filename = "TIC_plot.pdf",  plot = TIC_plot,
                path = "zzz",
                width = 5.89, height = 3.26,
                units = "in", dpi = 300, scale = 0.8)

ProcA2_ext <- extract_MS2(ProcA2_raw)


ProcA2_detected <- detect_mass(ProcA2_ext,
                               normalize = TRUE, # Allow normalization
                               min_int = 1) # 1% as minimum intensity



MS2_spectra <- plot_MS2spectra(ProcA2_detected) + ggplot2::theme_classic()

ggplot2::ggsave(filename = "MS2_spectra.pdf",  plot = MS2_spectra,
                path = "zzz",
                width = 5.76, height = 3.76,
                units = "in", dpi = 300, scale = 0.8)

# Rutin ----

# Importing the Spectrum of Ruinoin negative ionization mode
# and 20 eV as the collision energy
Rutin_file <- system.file("extdata",
                           "Rutin_neg_20eV.mzXML",
                           package = "MS2extract")

# Region of interest table (rt in seconds)
Rutin_data <- data.frame(Formula = "C27H30O16",Ionization_mode = "Negative",
                          min_rt = 160, max_rt = 175)

# Importing MS2 data
rutin_raw <- import_mzxml(Rutin_file, Rutin_data)

Rutin_extracted <- extract_MS2(rutin_raw)

# Get MS2 TIC
rutin_TIC <- get_TIC(rutin_raw)
# Plot MS2 TIC
rutin_TIC_results <- plot_tic(rutin_TIC)

rutin_TIC_plot <- rutin_TIC_results$TIC_plot + ggplot2::theme_classic()

ggplot2::ggsave(filename = "rutin_TIC_plot.pdf",  plot = rutin_TIC_plot,
                path = "zzz",
                width = 5.89, height = 3.26,
                units = "in", dpi = 300, scale = 0.8)


Rutin_ext <- extract_MS2(Rutin_ext)


Rutin_detected <- detect_mass(Rutin_ext,
                               normalize = TRUE, # Allow normalization
                               min_int = 1) # 1% as minimum intensity



MS2_spectra <- plot_MS2spectra(Rutin_detected) + ggplot2::theme_classic()

ggplot2::ggsave(filename = "Rutin_MS2_spectra.pdf",  plot = MS2_spectra,
                path = "zzz",
                width = 5.76, height = 3.76,
                units = "in", dpi = 300, scale = 0.8)
