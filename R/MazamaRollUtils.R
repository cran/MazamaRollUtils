#' Roll NowCast
#'
#' @description Apply the EPA NowCast algorithm to a numeric vector of hourly
#' particulate matter measurements.
#'
#' @details
#'
#' The EPA NowCast is a weighted average designed to emphasize more recent
#' hourly PM values while still using up to the previous 12 hours of data.
#' The weighting depends on how much concentrations vary within the window:
#' rapidly changing conditions place more weight on the most recent hours,
#' while stable conditions allow older hours to contribute more evenly.
#'
#' For every index in the incoming vector `x`, a value is returned that is the
#' NowCast associated with that hour.
#'
#' This calculation is always right-aligned:
#'
#' \itemize{
#'   \item{`[-----------*]` where `*` marks the hour receiving the NowCast value.}
#' }
#'
#' Early values use all available data back to the beginning of the series,
#' so the first 11 positions may be calculated from fewer than 12 hours.
#'
#' Missing values are allowed, but at least 2 valid values must be present in
#' the most recent 3 hours or the result for that index will be `NA`.
#'
#' A window whose valid values are all zero returns `0` (there is no rate of
#' change to weight by).
#'
#' Returned values are rounded to one decimal place.
#'
#' @param x Numeric vector of hourly PM measurements.
#'
#' @return Numeric vector of the same length as `x`.
#'
#' @examples
#' x <- c(10, 12, 11, 13, 15, 18, 20, 25, 30, 28, 26, 24, 22)
#' roll_nowcast(x)
#'
#' @export
roll_nowcast <- function(x) {

  if ( !is.atomic(x) || !is.numeric(x) || !is.null(dim(x)) ) {
    stop("'x' must be a numeric vector.")
  }

  result <- .roll_nowcast_cpp(as.numeric(x))

  return(result)
}

#' Roll Hampel
#'
#' @description Apply a moving-window Hampel function to a numeric vector.
#'
#' @details
#'
#' The Hampel filter is a robust outlier detector using Median Absolute Deviation (MAD).
#'
#' For every index in the incoming vector `x`, a value is returned that
#' is the Hampel function of all values in `x` that fall within a window of width
#' `width`. The score at each index compares the value at that index
#' (`x[i]`) with the median and MAD of its window: `abs(x[i] - median) /
#' (1.4826 * MAD)`. Larger scores indicate values less consistent with their
#' neighborhood. With `align = "center"` the tested value sits at the middle
#' of its window; with `align = "left"` or `"right"` it sits at the window
#' edge. Outlier detection (see [findOutliers()]) is normally done with the
#' default `align = "center"`.
#'
#' The `align` parameter determines the alignment of the return value
#' within the window. Thus:
#'
#' \itemize{
#'   \item{`align = "left"   [*------]` will cause the returned vector to have width - 1 `NA` values at the right end.}
#'   \item{`align = "center" [---*---]` will cause the returned vector to have `NA` values at either end as needed for centered alignment.}
#'   \item{`align = "right"  [------*]` will cause the returned vector to have width - 1 `NA` values at the left end.}
#' }
#'
#' For large vectors, the `by` parameter can be used to force the window
#' to jump ahead `by` indices for the next calculation. Indices that are
#' skipped over will be assigned `NA` values so that the return vector still has
#' the same length as the incoming vector. This can dramatically speed up
#' calculations for high resolution time series data.
#'
#' @param x Numeric vector.
#' @param width Integer width of the rolling window.
#' @param by Integer shift by which the window is moved each iteration.
#' @param align Character position of the return value within the window. One of:
#' `"left" | "center" | "right"`.
#' @param na.rm Logical specifying whether `NA` values should be removed
#' before the calculations within each window.
#' @param min_valid Integer minimum number of non-`NA` values that must fall
#' within a window for that window to return a non-`NA` result. Supplying
#' `min_valid` implies `NA`-tolerant counting within each window (as if
#' `na.rm = TRUE`), regardless of `na.rm`. The default, `NULL`, applies no
#' minimum.
#'
#' @return Numeric vector of the same length as `x`.
#'
#' @examples
#' x <- c(0, 0, 0, 1, 1, 2, 2, 4, 6, 9, 0, 0, 0)
#' roll_hampel(x, 3)
#'
#' @export
roll_hampel <- function(
    x,
    width = 1L,
    by = 1L,
    align = c("center", "left", "right"),
    na.rm = FALSE,
    min_valid = NULL
) {

  args <- .validateRollArgs(
    x = x,
    width = width,
    by = by,
    align = align,
    na.rm = na.rm,
    min_valid = min_valid
  )

  result <- .roll_hampel_cpp(
    args$x,
    args$width,
    args$by,
    args$align,
    args$na.rm
  )

  result <- .applyMinValid(
    result, args$x, args$width, args$by, args$align, args$min_valid
  )

  return(result)
}

#' Roll MAD
#'
#' @description Apply a moving-window Median Absolute Deviation function to a numeric vector.
#'
#' @details
#'
#' For every index in the incoming vector `x`, a value is returned that
#' is the Median Absolute Deviation (MAD) of all values in `x` that fall within
#' a window of width `width`.
#'
#' The value returned is the *unscaled* MAD -- the median of the absolute
#' deviations from the window median, with no consistency constant applied.
#' This is equivalent to `stats::mad(window, constant = 1)`. Note that
#' `stats::mad()` uses `constant = 1.4826` by default, so values from
#' `roll_MAD()` are smaller than the `stats::mad()` default by that factor.
#'
#' The `align` parameter determines the alignment of the return value
#' within the window. Thus:
#'
#' \itemize{
#'   \item{`align = "left"   [*------]` will cause the returned vector to have width - 1 `NA` values at the right end.}
#'   \item{`align = "center" [---*---]` will cause the returned vector to have `NA` values at either end as needed for centered alignment.}
#'   \item{`align = "right"  [------*]` will cause the returned vector to have width - 1 `NA` values at the left end.}
#' }
#'
#' For large vectors, the `by` parameter can be used to force the window
#' to jump ahead `by` indices for the next calculation. Indices that are
#' skipped over will be assigned `NA` values so that the return vector still has
#' the same length as the incoming vector. This can dramatically speed up
#' calculations for high resolution time series data.
#'
#' @param x Numeric vector.
#' @param width Integer width of the rolling window.
#' @param by Integer shift by which the window is moved each iteration.
#' @param align Character position of the return value within the window. One of:
#' `"left" | "center" | "right"`.
#' @param na.rm Logical specifying whether `NA` values should be removed
#' before the calculations within each window.
#' @param min_valid Integer minimum number of non-`NA` values that must fall
#' within a window for that window to return a non-`NA` result. Supplying
#' `min_valid` implies `NA`-tolerant counting within each window (as if
#' `na.rm = TRUE`), regardless of `na.rm`. The default, `NULL`, applies no
#' minimum.
#'
#' @return Numeric vector of the same length as `x`.
#'
#' @examples
#' # Wikipedia example
#' x <- c(0, 0, 0, 1, 1, 2, 2, 4, 6, 9, 0, 0, 0)
#' roll_MAD(x, 3)
#' roll_MAD(x, 5)
#' roll_MAD(x, 7)
#'
#' @export
roll_MAD <- function(
    x,
    width = 1L,
    by = 1L,
    align = c("center", "left", "right"),
    na.rm = FALSE,
    min_valid = NULL
) {

  args <- .validateRollArgs(
    x = x,
    width = width,
    by = by,
    align = align,
    na.rm = na.rm,
    min_valid = min_valid
  )

  result <- .roll_MAD_cpp(
    args$x,
    args$width,
    args$by,
    args$align,
    args$na.rm
  )

  result <- .applyMinValid(
    result, args$x, args$width, args$by, args$align, args$min_valid
  )

  return(result)
}

#' Roll Max
#'
#' @description Apply a moving-window maximum function to a numeric vector.
#'
#' @details
#'
#' For every index in the incoming vector `x`, a value is returned that
#' is the maximum of all values in `x` that fall within a window of width
#' `width`.
#'
#' The `align` parameter determines the alignment of the return value
#' within the window. Thus:
#'
#' \itemize{
#'   \item{`align = "left"   [*------]` will cause the returned vector to have width - 1 `NA` values at the right end.}
#'   \item{`align = "center" [---*---]` will cause the returned vector to have `NA` values at either end as needed for centered alignment.}
#'   \item{`align = "right"  [------*]` will cause the returned vector to have width - 1 `NA` values at the left end.}
#' }
#'
#' For large vectors, the `by` parameter can be used to force the window
#' to jump ahead `by` indices for the next calculation. Indices that are
#' skipped over will be assigned `NA` values so that the return vector still has
#' the same length as the incoming vector. This can dramatically speed up
#' calculations for high resolution time series data.
#'
#' @param x Numeric vector.
#' @param width Integer width of the rolling window.
#' @param by Integer shift by which the window is moved each iteration.
#' @param align Character position of the return value within the window. One of:
#' `"left" | "center" | "right"`.
#' @param na.rm Logical specifying whether `NA` values should be removed
#' before the calculations within each window.
#' @param min_valid Integer minimum number of non-`NA` values that must fall
#' within a window for that window to return a non-`NA` result. Supplying
#' `min_valid` implies `NA`-tolerant counting within each window (as if
#' `na.rm = TRUE`), regardless of `na.rm`. The default, `NULL`, applies no
#' minimum.
#'
#' @return Numeric vector of the same length as `x`.
#'
#' @examples
#' # Example air quality time series
#' t <- example_pm25$datetime
#' x <- example_pm25$pm25
#'
#' plot(t, x, pch = 16, cex = 0.5)
#' lines(t, roll_max(x, width = 12), col = "red")
#' lines(t, roll_min(x, width = 12), col = "deepskyblue")
#' title("12-hr Rolling Max and Min")
#'
#' plot(t, x, pch = 16, cex = 0.5)
#' points(t, roll_max(x, width = 12, na.rm = TRUE),
#'        pch = 16, col = "red")
#' points(t, roll_max(x, width = 12, na.rm = FALSE),
#'        pch = 16, col = adjustcolor("black", 0.4))
#' legend("topright", pch = c(1, 16),
#'        col = c("red", adjustcolor("black", 0.4)),
#'        legend = c("na.rm = TRUE", "na.rm = FALSE"))
#' title("12-hr Rolling max with/out na.rm")
#'
#' @export
roll_max <- function(
    x,
    width = 1L,
    by = 1L,
    align = c("center", "left", "right"),
    na.rm = FALSE,
    min_valid = NULL
) {

  args <- .validateRollArgs(
    x = x,
    width = width,
    by = by,
    align = align,
    na.rm = na.rm,
    min_valid = min_valid
  )

  result <- .roll_max_cpp(
    args$x,
    args$width,
    args$by,
    args$align,
    args$na.rm
  )

  result <- .applyMinValid(
    result, args$x, args$width, args$by, args$align, args$min_valid
  )

  return(result)
}

#' Roll Mean
#'
#' @description Apply a moving-window mean function to a numeric vector.
#'
#' @details
#'
#' For every index in the incoming vector `x`, a value is returned that
#' is the mean of all values in `x` that fall within a window of width
#' `width`.
#'
#' The `align` parameter determines the alignment of the return value
#' within the window. Thus:
#'
#' \itemize{
#'   \item{`align = "left"   [*------]` will cause the returned vector to have width - 1 `NA` values at the right end.}
#'   \item{`align = "center" [---*---]` will cause the returned vector to have `NA` values at either end as needed for centered alignment.}
#'   \item{`align = "right"  [------*]` will cause the returned vector to have width - 1 `NA` values at the left end.}
#' }
#'
#' For large vectors, the `by` parameter can be used to force the window
#' to jump ahead `by` indices for the next calculation. Indices that are
#' skipped over will be assigned `NA` values so that the return vector still has
#' the same length as the incoming vector. This can dramatically speed up
#' calculations for high resolution time series data.
#'
#' The `roll_mean()` function supports an additional `weights`
#' argument that can be used to calculate a weighted moving average,
#' a convolution of the incoming data with the kernel provided in `weights`.
#'
#' @param x Numeric vector.
#' @param width Integer width of the rolling window.
#' @param by Integer shift by which the window is moved each iteration.
#' @param align Character position of the return value within the window. One of:
#' `"left" | "center" | "right"`.
#' @param na.rm Logical specifying whether `NA` values should be removed
#' before the calculations within each window.
#' @param weights Numeric vector of length `width` specifying each window
#' index weight. If `NULL`, unit weights are used.
#' @param min_valid Integer minimum number of non-`NA` values that must fall
#' within a window for that window to return a non-`NA` result. Supplying
#' `min_valid` implies `NA`-tolerant counting within each window (as if
#' `na.rm = TRUE`), regardless of `na.rm`. The default, `NULL`, applies no
#' minimum.
#'
#' @return Numeric vector of the same length as `x`.
#'
#' @examples
#' # Example air quality time series
#' t <- example_pm25$datetime
#' x <- example_pm25$pm25
#'
#' plot(t, x, pch = 16, cex = 0.5)
#' lines(t, roll_mean(x, width = 3), col = "goldenrod")
#' lines(t, roll_mean(x, width = 23), col = "purple")
#' legend("topright", lty = c(1, 1),
#'        col = c("goldenrod", "purple"),
#'        legend = c("3-hr mean", "23-hr mean"))
#' title("3- and 23-hr Rolling mean")
#'
#' @export
roll_mean <- function(
    x,
    width = 1L,
    by = 1L,
    align = c("center", "left", "right"),
    na.rm = FALSE,
    weights = NULL,
    min_valid = NULL
) {

  args <- .validateRollArgs(
    x = x,
    width = width,
    by = by,
    align = align,
    na.rm = na.rm,
    weights = weights,
    min_valid = min_valid
  )

  result <- .roll_mean_cpp(
    args$x,
    args$width,
    args$by,
    args$align,
    args$na.rm,
    args$weights
  )

  result <- .applyMinValid(
    result, args$x, args$width, args$by, args$align, args$min_valid
  )

  return(result)
}

#' Roll Median
#'
#' @description Apply a moving-window median function to a numeric vector.
#'
#' @details
#'
#' For every index in the incoming vector `x`, a value is returned that
#' is the median of all values in `x` that fall within a window of width
#' `width`.
#'
#' The `align` parameter determines the alignment of the return value
#' within the window. Thus:
#'
#' \itemize{
#'   \item{`align = "left"   [*------]` will cause the returned vector to have width - 1 `NA` values at the right end.}
#'   \item{`align = "center" [---*---]` will cause the returned vector to have `NA` values at either end as needed for centered alignment.}
#'   \item{`align = "right"  [------*]` will cause the returned vector to have width - 1 `NA` values at the left end.}
#' }
#'
#' For large vectors, the `by` parameter can be used to force the window
#' to jump ahead `by` indices for the next calculation. Indices that are
#' skipped over will be assigned `NA` values so that the return vector still has
#' the same length as the incoming vector. This can dramatically speed up
#' calculations for high resolution time series data.
#'
#' @param x Numeric vector.
#' @param width Integer width of the rolling window.
#' @param by Integer shift by which the window is moved each iteration.
#' @param align Character position of the return value within the window. One of:
#' `"left" | "center" | "right"`.
#' @param na.rm Logical specifying whether `NA` values should be removed
#' before the calculations within each window.
#' @param min_valid Integer minimum number of non-`NA` values that must fall
#' within a window for that window to return a non-`NA` result. Supplying
#' `min_valid` implies `NA`-tolerant counting within each window (as if
#' `na.rm = TRUE`), regardless of `na.rm`. The default, `NULL`, applies no
#' minimum.
#'
#' @return Numeric vector of the same length as `x`.
#'
#' @examples
#' # Example air quality time series
#' t <- example_pm25$datetime
#' x <- example_pm25$pm25
#'
#' plot(t, x, pch = 16, cex = 0.5)
#' lines(t, roll_median(x, width = 3), col = "goldenrod")
#' lines(t, roll_median(x, width = 23), col = "purple")
#' legend("topright", lty = c(1, 1),
#'        col = c("goldenrod", "purple"),
#'        legend = c("3-hr median", "23-hr median"))
#' title("3- and 23-hr Rolling median")
#'
#' @export
roll_median <- function(
    x,
    width = 1L,
    by = 1L,
    align = c("center", "left", "right"),
    na.rm = FALSE,
    min_valid = NULL
) {

  args <- .validateRollArgs(
    x = x,
    width = width,
    by = by,
    align = align,
    na.rm = na.rm,
    min_valid = min_valid
  )

  result <- .roll_median_cpp(
    args$x,
    args$width,
    args$by,
    args$align,
    args$na.rm
  )

  result <- .applyMinValid(
    result, args$x, args$width, args$by, args$align, args$min_valid
  )

  return(result)
}

#' Roll Min
#'
#' @description Apply a moving-window minimum function to a numeric vector.
#'
#' @details
#'
#' For every index in the incoming vector `x`, a value is returned that
#' is the minimum of all values in `x` that fall within a window of width
#' `width`.
#'
#' The `align` parameter determines the alignment of the return value
#' within the window. Thus:
#'
#' \itemize{
#'   \item{`align = "left"   [*------]` will cause the returned vector to have width - 1 `NA` values at the right end.}
#'   \item{`align = "center" [---*---]` will cause the returned vector to have `NA` values at either end as needed for centered alignment.}
#'   \item{`align = "right"  [------*]` will cause the returned vector to have width - 1 `NA` values at the left end.}
#' }
#'
#' For large vectors, the `by` parameter can be used to force the window
#' to jump ahead `by` indices for the next calculation. Indices that are
#' skipped over will be assigned `NA` values so that the return vector still has
#' the same length as the incoming vector. This can dramatically speed up
#' calculations for high resolution time series data.
#'
#' @param x Numeric vector.
#' @param width Integer width of the rolling window.
#' @param by Integer shift by which the window is moved each iteration.
#' @param align Character position of the return value within the window. One of:
#' `"left" | "center" | "right"`.
#' @param na.rm Logical specifying whether `NA` values should be removed
#' before the calculations within each window.
#' @param min_valid Integer minimum number of non-`NA` values that must fall
#' within a window for that window to return a non-`NA` result. Supplying
#' `min_valid` implies `NA`-tolerant counting within each window (as if
#' `na.rm = TRUE`), regardless of `na.rm`. The default, `NULL`, applies no
#' minimum.
#'
#' @return Numeric vector of the same length as `x`.
#'
#' @examples
#' # Example air quality time series
#' t <- example_pm25$datetime
#' x <- example_pm25$pm25
#'
#' plot(t, x, pch = 16, cex = 0.5)
#' lines(t, roll_max(x, width = 12), col = "red")
#' lines(t, roll_min(x, width = 12), col = "deepskyblue")
#' title("12-hr Rolling Max and Min")
#'
#' plot(t, x, pch = 16, cex = 0.5)
#' points(t, roll_min(x, width = 12, na.rm = TRUE),
#'        pch = 16, col = "deepskyblue")
#' points(t, roll_min(x, width = 12, na.rm = FALSE),
#'        pch = 16, col = adjustcolor("black", 0.4))
#' legend("topright", pch = c(16, 16),
#'        col = c("deepskyblue", adjustcolor("black", 0.4)),
#'        legend = c("na.rm = TRUE", "na.rm = FALSE"))
#' title("12-hr Rolling min with/out na.rm")
#'
#' @export
roll_min <- function(
    x,
    width = 1L,
    by = 1L,
    align = c("center", "left", "right"),
    na.rm = FALSE,
    min_valid = NULL
) {

  args <- .validateRollArgs(
    x = x,
    width = width,
    by = by,
    align = align,
    na.rm = na.rm,
    min_valid = min_valid
  )

  result <- .roll_min_cpp(
    args$x,
    args$width,
    args$by,
    args$align,
    args$na.rm
  )

  result <- .applyMinValid(
    result, args$x, args$width, args$by, args$align, args$min_valid
  )

  return(result)
}

#' Roll Product
#'
#' @description Apply a moving-window product function to a numeric vector.
#'
#' @details
#'
#' For every index in the incoming vector `x`, a value is returned that
#' is the product of all values in `x` that fall within a window of width
#' `width`.
#'
#' The `align` parameter determines the alignment of the return value
#' within the window. Thus:
#'
#' \itemize{
#'   \item{`align = "left"   [*------]` will cause the returned vector to have width - 1 `NA` values at the right end.}
#'   \item{`align = "center" [---*---]` will cause the returned vector to have `NA` values at either end as needed for centered alignment.}
#'   \item{`align = "right"  [------*]` will cause the returned vector to have width - 1 `NA` values at the left end.}
#' }
#'
#' For large vectors, the `by` parameter can be used to force the window
#' to jump ahead `by` indices for the next calculation. Indices that are
#' skipped over will be assigned `NA` values so that the return vector still has
#' the same length as the incoming vector. This can dramatically speed up
#' calculations for high resolution time series data.
#'
#' @param x Numeric vector.
#' @param width Integer width of the rolling window.
#' @param by Integer shift by which the window is moved each iteration.
#' @param align Character position of the return value within the window. One of:
#' `"left" | "center" | "right"`.
#' @param na.rm Logical specifying whether `NA` values should be removed
#' before the calculations within each window.
#' @param min_valid Integer minimum number of non-`NA` values that must fall
#' within a window for that window to return a non-`NA` result. Supplying
#' `min_valid` implies `NA`-tolerant counting within each window (as if
#' `na.rm = TRUE`), regardless of `na.rm`. The default, `NULL`, applies no
#' minimum.
#'
#' @return Numeric vector of the same length as `x`.
#'
#' @examples
#' # Example air quality time series
#' t <- example_pm25$datetime
#' x <- example_pm25$pm25
#'
#' x[1:10]
#' roll_prod(x, width = 5)[1:10]
#'
#' @export
roll_prod <- function(
    x,
    width = 1L,
    by = 1L,
    align = c("center", "left", "right"),
    na.rm = FALSE,
    min_valid = NULL
) {

  args <- .validateRollArgs(
    x = x,
    width = width,
    by = by,
    align = align,
    na.rm = na.rm,
    min_valid = min_valid
  )

  result <- .roll_prod_cpp(
    args$x,
    args$width,
    args$by,
    args$align,
    args$na.rm
  )

  result <- .applyMinValid(
    result, args$x, args$width, args$by, args$align, args$min_valid
  )

  return(result)
}

#' Roll Standard Deviation
#'
#' @description Apply a moving-window standard deviation function to a
#' numeric vector.
#'
#' @details
#'
#' For every index in the incoming vector `x`, a value is returned that
#' is the standard deviation of all values in `x` that fall within a window of
#' width `width`.
#'
#' The `align` parameter determines the alignment of the return value
#' within the window. Thus:
#'
#' \itemize{
#'   \item{`align = "left"   [*------]` will cause the returned vector to have width - 1 `NA` values at the right end.}
#'   \item{`align = "center" [---*---]` will cause the returned vector to have `NA` values at either end as needed for centered alignment.}
#'   \item{`align = "right"  [------*]` will cause the returned vector to have width - 1 `NA` values at the left end.}
#' }
#'
#' For large vectors, the `by` parameter can be used to force the window
#' to jump ahead `by` indices for the next calculation. Indices that are
#' skipped over will be assigned `NA` values so that the return vector still has
#' the same length as the incoming vector. This can dramatically speed up
#' calculations for high resolution time series data.
#'
#' @note A `na.rm` argument is not provided for `roll_sd()` because the
#' statistical meaning of standard deviation computed from partially missing
#' windows may be ambiguous.
#'
#' @param x Numeric vector.
#' @param width Integer width of the rolling window.
#' @param by Integer shift by which the window is moved each iteration.
#' @param align Character position of the return value within the window. One of:
#' `"left" | "center" | "right"`.
#'
#' @return Numeric vector of the same length as `x`.
#'
#' @examples
#' # Example air quality time series
#' t <- example_pm25$datetime
#' x <- example_pm25$pm25
#'
#' x[1:10]
#' roll_sd(x, width = 5)[1:10]
#'
#' @export
roll_sd <- function(
    x,
    width = 1L,
    by = 1L,
    align = c("center", "left", "right")
) {

  args <- .validateRollArgs(
    x = x,
    width = width,
    by = by,
    align = align
  )

  result <- .roll_sd_cpp(
    args$x,
    args$width,
    args$by,
    args$align,
    FALSE
  )

  return(result)
}

#' Roll Sum
#'
#' @description Apply a moving-window sum to a numeric vector.
#'
#' @details
#'
#' For every index in the incoming vector `x`, a value is returned that
#' is the sum of all values in `x` that fall within a window of width
#' `width`.
#'
#' The `align` parameter determines the alignment of the return value
#' within the window. Thus:
#'
#' \itemize{
#'   \item{`align = "left"   [*------]` will cause the returned vector to have width - 1 `NA` values at the right end.}
#'   \item{`align = "center" [---*---]` will cause the returned vector to have `NA` values at either end as needed for centered alignment.}
#'   \item{`align = "right"  [------*]` will cause the returned vector to have width - 1 `NA` values at the left end.}
#' }
#'
#' For large vectors, the `by` parameter can be used to force the window
#' to jump ahead `by` indices for the next calculation. Indices that are
#' skipped over will be assigned `NA` values so that the return vector still has
#' the same length as the incoming vector. This can dramatically speed up
#' calculations for high resolution time series data.
#'
#' @param x Numeric vector.
#' @param width Integer width of the rolling window.
#' @param by Integer shift by which the window is moved each iteration.
#' @param align Character position of the return value within the window. One of:
#' `"left" | "center" | "right"`.
#' @param na.rm Logical specifying whether `NA` values should be removed
#' before the calculations within each window.
#' @param min_valid Integer minimum number of non-`NA` values that must fall
#' within a window for that window to return a non-`NA` result. Supplying
#' `min_valid` implies `NA`-tolerant counting within each window (as if
#' `na.rm = TRUE`), regardless of `na.rm`. The default, `NULL`, applies no
#' minimum.
#'
#' @return Numeric vector of the same length as `x`.
#'
#' @examples
#' # Example air quality time series
#' t <- example_pm25$datetime
#' x <- example_pm25$pm25
#'
#' x[1:10]
#' roll_sum(x, width = 5)[1:10]
#'
#' @export
roll_sum <- function(
    x,
    width = 1L,
    by = 1L,
    align = c("center", "left", "right"),
    na.rm = FALSE,
    min_valid = NULL
) {

  args <- .validateRollArgs(
    x = x,
    width = width,
    by = by,
    align = align,
    na.rm = na.rm,
    min_valid = min_valid
  )

  result <- .roll_sum_cpp(
    args$x,
    args$width,
    args$by,
    args$align,
    args$na.rm
  )

  result <- .applyMinValid(
    result, args$x, args$width, args$by, args$align, args$min_valid
  )

  return(result)
}

#' Roll Variance
#'
#' @description Apply a moving-window variance function to a numeric vector.
#'
#' @details
#'
#' For every index in the incoming vector `x`, a value is returned that
#' is the variance of all values in `x` that fall within a window of width
#' `width`.
#'
#' The `align` parameter determines the alignment of the return value
#' within the window. Thus:
#'
#' \itemize{
#'   \item{`align = "left"   [*------]` will cause the returned vector to have width - 1 `NA` values at the right end.}
#'   \item{`align = "center" [---*---]` will cause the returned vector to have `NA` values at either end as needed for centered alignment.}
#'   \item{`align = "right"  [------*]` will cause the returned vector to have width - 1 `NA` values at the left end.}
#' }
#'
#' For large vectors, the `by` parameter can be used to force the window
#' to jump ahead `by` indices for the next calculation. Indices that are
#' skipped over will be assigned `NA` values so that the return vector still has
#' the same length as the incoming vector. This can dramatically speed up
#' calculations for high resolution time series data.
#'
#' @note A `na.rm` argument is not provided for `roll_var()` because the
#' statistical meaning of variance computed from partially missing windows may
#' be ambiguous.
#'
#' @param x Numeric vector.
#' @param width Integer width of the rolling window.
#' @param by Integer shift by which the window is moved each iteration.
#' @param align Character position of the return value within the window. One of:
#' `"left" | "center" | "right"`.
#'
#' @return Numeric vector of the same length as `x`.
#'
#' @examples
#' # Example air quality time series
#' t <- example_pm25$datetime
#' x <- example_pm25$pm25
#'
#' x[1:10]
#' roll_var(x, width = 5)[1:10]
#'
#' @export
roll_var <- function(
    x,
    width = 1L,
    by = 1L,
    align = c("center", "left", "right")
) {

  args <- .validateRollArgs(
    x = x,
    width = width,
    by = by,
    align = align
  )

  result <- .roll_var_cpp(
    args$x,
    args$width,
    args$by,
    args$align,
    FALSE
  )

  return(result)
}

.validateRollArgs <- function(
    x,
    width,
    by,
    align,
    na.rm = NULL,
    weights = NULL,
    min_valid = NULL
) {

  if ( !is.atomic(x) || !is.numeric(x) || !is.null(dim(x)) ) {
    stop("'x' must be a numeric vector.")
  }

  if ( length(width) != 1 || !is.numeric(width) || is.na(width) ||
       !is.finite(width) || width < 1 || width != as.integer(width) ) {
    stop("'width' must be a single positive integer.")
  }

  if ( width > length(x) ) {
    stop("'width' cannot be larger than 'x'.")
  }

  if ( length(by) != 1 || !is.numeric(by) || is.na(by) ||
       !is.finite(by) || by < 1 || by != as.integer(by) ) {
    stop("'by' must be a single positive integer.")
  }

  if ( by > length(x) ) {
    stop("'by' cannot be larger than 'x'.")
  }

  align <- match.arg(align, c("center", "left", "right"))

  if ( !is.null(na.rm) ) {
    if ( !is.logical(na.rm) || length(na.rm) != 1 || is.na(na.rm) ) {
      stop("'na.rm' must be TRUE or FALSE.")
    }
  }

  if ( !is.null(min_valid) ) {
    if ( length(min_valid) != 1 || !is.numeric(min_valid) || is.na(min_valid) ||
         !is.finite(min_valid) || min_valid < 1 ||
         min_valid != as.integer(min_valid) ) {
      stop("'min_valid' must be a single positive integer.")
    }
    if ( min_valid > as.integer(width) ) {
      stop("'min_valid' cannot be larger than 'width'.")
    }
    # A minimum-valid-count threshold implies NA-tolerant counting within each
    # window, so the calculation proceeds as if na.rm = TRUE regardless of the
    # value supplied.
    na.rm <- TRUE
  }

  if ( !is.null(weights) ) {
    if ( !is.atomic(weights) || !is.numeric(weights) || !is.null(dim(weights)) ) {
      stop("'weights' must be NULL or a numeric vector.")
    }
    if ( anyNA(weights) || any(!is.finite(weights)) ) {
      stop("'weights' must not contain NA, NaN, or infinite values.")
    }
    if ( any(weights < 0) ) {
      stop("'weights' must not contain negative values.")
    }
    if ( length(weights) != as.integer(width) ) {
      stop("'weights' must have length equal to 'width'.")
    }
  }

  return(list(
    x = x,
    width = as.integer(width),
    by = as.integer(by),
    align = align,
    na.rm = na.rm,
    weights = weights,
    min_valid = if (is.null(min_valid)) NULL else as.integer(min_valid)
  ))
}

# Set to NA any result position whose window held fewer than 'min_valid' non-NA
# values. The rolling count of valid values is computed with roll_sum() on a
# 0/1 indicator, which reuses the same window geometry (width, by, align) and
# end padding, so 'valid_count' is NA in exactly the positions where 'result'
# is already NA.
.applyMinValid <- function(result, x, width, by, align, min_valid) {

  if ( is.null(min_valid) ) {
    return(result)
  }

  valid_count <- roll_sum(
    as.numeric(!is.na(x)),
    width = width,
    by = by,
    align = align,
    na.rm = FALSE
  )

  result[!is.na(valid_count) & valid_count < min_valid] <- NA_real_

  return(result)
}
