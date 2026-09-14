#' Package imports
#'
#' @importFrom grDevices col2rgb dev.size gray rgb
#' @importFrom graphics abline axis layout mtext par plot.new rect text
#' @importFrom stats aggregate dist predict rnorm
#' @name mand-imports
#' @noRd
#' @keywords internal
NULL

# Internal heat-style color palette for mand.
#
# This implementation uses grDevices::colorRampPalette() with control
# colors selected independently for the mand package. It is not exported.
.mand_hotmetal <- function(n = 64L) {
  if (length(n) != 1L || is.na(n) || !is.finite(n) ||
      n < 1L || n != as.integer(n)) {
    stop("n must be one positive integer.")
  }
  n <- as.integer(n)
  controls <- c(
    "#080000", "#700000", "#D80000", "#FF4800",
    "#FF9C00", "#FFE200", "#FFFF70", "#FFFFF8"
  )
  if (n == 1L) return(controls[1L])
  grDevices::colorRampPalette(
    controls,
    space = "rgb",
    interpolate = "linear"
  )(n)
}

