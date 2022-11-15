#' @title Gapfill
#' @description This function fills the gaps for the time-gas data by applying linear functions. This step should be done prior to curve-fitting, to balance to weights of the measured points.
#' @param Time Numeric vector. Time points
#' @param gasV Numeric vector. Gas volume data.
#' @param Total_t A number. Total time, can be greater than the max(Time).
#' @param gap A number. The minimum time gap.
#'
#' @return Returns a list. List[[1]] is the filled time points. List[[2]] is the filled gas volume.
#' @export
#'
#' @examples Time <- c(0,2.731,3.618,4.101,4.426,4.695,
#'                     4.996,5.286,5.602,5.933,6.274,6.612,
#'                     6.996,7.451,8.136,10.552,46.771)
#' gasV <- c(0,2.13,  4.16,  6.24,  8.40, 10.55, 12.69, 14.94, 17.25, 19.41, 21.54, 23.64,
#'           25.71, 27.80, 29.85, 31.85, 33.85)
#' Filled <- Gapfill(Time, gasV, 48, 5/60)
#' plot(Filled[[1]], Filled[[2]])

Gapfill <- function(Time, gasV, Total_t, gap){
  filled <- list()
  if (max(Time) < Total_t) {
    Time <- append(Time, Total_t)
    gasV <- append(gasV, max(gasV))
  }
  filled[[1]] <- approx(Time, gasV, n = floor(Total_t/gap))$x
  filled[[2]] <- approx(Time, gasV, n = floor(Total_t/gap))$y
  names(filled) <- c("FilledTime", "FilledgasV")
  return(filled)
}


