#' Calculate the  *m*/*z* difference for a given mass and ppm tolerance
#'
#' Given a theoretical  *m*/*z* value and ppm tolerance, it computes the  *m*/*z*
#' deviation.
#'
#' @param mz a double, the theoretical m/z value.
#' @param ppm a double, the ppm error tolerance.
#'
#' @return a double with the m/z value deviation from the given
#' theoretical m/z value
#'
#' @export
#'
#' @examples
#' get_ppm(mz = 355.1023, ppm = 20)
#'
#' # Chlorogenic acid [M+H]+ = 355.1023 m/z
#' chlorogenic_acid_pos <- 355.1023
#' ppm_error <- 10
#' get_ppm(mz = chlorogenic_acid_pos, ppm = ppm_error)
get_ppm <- function(mz, ppm = 10) {
  mz_observed <- (ppm * mz / 1E6) + mz
  mz_difference <- abs(mz_observed - mz)
  mz_difference
}

#' Calculate the *m*/*z* upper and lower values limits given a ppm tolerance
#'
#' Given a theoretical *m*/*z* value, and a ppm error tolerance,
#' this function calculates the lower and upper *m*/*z* values.
#'
#' @param mz a double, the theoretical *m*/*z* value.
#' @param ppm a double, the mass error in ppm
#'
#' @return a vector of two elements containing the upper and lower *m*/*z*values
#'   given the provided mass error tolerance.
#'
#' @export
#'
#' @examples
#'
#' chlorogenic_acid_pos <- 355.1023
#' ppm_error <- 10
#' ppm_range(mz = chlorogenic_acid_pos, ppm = ppm_error)
ppm_range <- function(mz, ppm) {
  ppm_error <- get_ppm(mz = mz, ppm = ppm)
  ppm_range_value <- c(mz - ppm_error, mz + ppm_error)
  ppm_range_value
}
