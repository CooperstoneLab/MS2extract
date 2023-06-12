test_that("ppm calculation works", {
  chlorogenic_acid_pos <- 355.1023
  ppm_error <- 10
  ppm_value <- get_ppm(mz = chlorogenic_acid_pos, ppm = ppm_error)
  # Calculate the mass diference
  expect_equal(ppm_value, 0.00355, tolerance = 0.0005)

  # Calculate ranges
  ppm_baund <- ppm_range(mz = chlorogenic_acid_pos, ppm = ppm_error)

  expect_true(length(ppm_baund) == 2)
})
