# MazamaRollUtils

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/MazamaRollUtils)](https://cran.r-project.org/package=MazamaRollUtils)
[![Downloads](https://cranlogs.r-pkg.org/badges/MazamaRollUtils)](https://cran.r-project.org/package=MazamaRollUtils)

MazamaRollUtils provides fast rolling-window ("moving") functions for numeric
vectors, backed by compiled C++ (Rcpp). It covers the familiar rolling
statistics — mean, median, min, max, sum, product, standard deviation,
variance — along with a Median Absolute Deviation, a Hampel filter, and the
US EPA NowCast.

The package is designed for efficient processing of environmental time series
such as hourly air-quality data. It deliberately operates on plain numeric
vectors with no underlying data model, so it composes with any workflow, and
every rolling function returns a vector the same length as its input.

## Installation

Install the released version from CRAN:

```r
install.packages("MazamaRollUtils")
```

Install the development version from GitHub:

```r
remotes::install_github("MazamaScience/MazamaRollUtils")
```

## Example

Apply a rolling mean and a rolling max/min envelope to the hourly PM2.5
air-quality series included with the package:

```r
library(MazamaRollUtils)

t <- example_pm25$datetime
x <- example_pm25$pm25

plot(t, x, pch = 16, cex = 0.5, col = "gray60")
lines(t, roll_mean(x, width = 12), col = "black", lwd = 2)
lines(t, roll_max(x, width = 12), col = "salmon")
lines(t, roll_min(x, width = 12), col = "steelblue")
```

## Overview

- **Rolling statistics** — `roll_mean()`, `roll_median()`, `roll_max()`,
  `roll_min()`, `roll_sum()`, `roll_prod()`, `roll_sd()`, `roll_var()`
- **Robust measures and outlier detection** — `roll_MAD()` (Median Absolute
  Deviation), `roll_hampel()` (Hampel filter), `findOutliers()` (indices of
  outliers flagged by a rolling Hampel filter)
- **Domain-specific calculations** — `roll_nowcast()` (US EPA NowCast for
  hourly particulate matter)

The `roll_*()` functions share the arguments `width`, `by`, `align`, and, where
statistically meaningful, `na.rm` and `min_valid` (a minimum count of non-`NA`
values per window); `roll_mean()` additionally accepts `weights` for a weighted
moving average. See the
[introductory vignette](https://mazamascience.github.io/MazamaRollUtils/articles/MazamaRollUtils.html)
and the [function reference](https://mazamascience.github.io/MazamaRollUtils/reference/index.html)
for argument details and return-value conventions.

## Background

Analysis of time series data often involves "rolling" calculations such as a
moving average. These are simple to express in R but slow, so compiled versions
of the common functions are valuable. Several R packages already provide some
of this functionality:

* [zoo](https://cran.r-project.org/package=zoo) — widely used package for
  ordered observations and rolling calculations
* [seismicRoll](https://cran.r-project.org/package=seismicRoll) — rolling
  functions focused on seismology
* [RcppRoll](https://cran.r-project.org/package=RcppRoll) — rolling functions
  for basic statistics

MazamaRollUtils exists to build up a suite of rolling functions useful in
environmental time series analysis, available in a neutral environment with no
underlying data model and usable by data analysts at any level of R expertise.

## Documentation and Help

- [Package website](https://mazamascience.github.io/MazamaRollUtils/)
- [Introduction to MazamaRollUtils](https://mazamascience.github.io/MazamaRollUtils/articles/MazamaRollUtils.html)
- [Function reference](https://mazamascience.github.io/MazamaRollUtils/reference/index.html)
- [GitHub Issues](https://github.com/MazamaScience/MazamaRollUtils/issues)
- [GitHub Discussions](https://github.com/MazamaScience/MazamaRollUtils/discussions)

## Citation

For citation information, use:

```r
citation("MazamaRollUtils")
```

## Acknowledgements

This project is supported by the [USFS AirFire](https://www.airfire.org) team.

## License

MazamaRollUtils is released under the GPL-3 license.
