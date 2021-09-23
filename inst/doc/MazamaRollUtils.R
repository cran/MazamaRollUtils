## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE,
  fig.width = 7,
  fig.height = 7
)

## ----example_1----------------------------------------------------------------
library(MazamaRollUtils)

# Extract vectors from our example dataset
t <- example_pm25$datetime
x <- example_pm25$pm25

# Plot with 3- and 24-hr rolling means
layout(matrix(seq(2)))
plot(t, x, pch = 16, cex = 0.5)
lines(t, roll_mean(x, width = 3), col = 'red')
title("3-hour Rolling Mean")
plot(t, x, pch = 16, cex = 0.5)
lines(t, roll_mean(x, width = 24), col = 'red')
title("24-hour Rolling Mean")
layout(1)

## ----example_2----------------------------------------------------------------
library(MazamaRollUtils)

# Extract vectors from our example dataset
t <- example_pm25$datetime
x <- example_pm25$pm25

# Calculate the left-aligned 24-hr max every hour, ignoring NA values
max_24hr <- roll_max(x, width = 24, align = "left", by = 1, na.rm = TRUE)

# Calculate the left-aligned daily max once every 24 hours, ignoring NA values
max_daily_day <- roll_max(x, width = 24, align = "left", by = 24, na.rm = TRUE)

# Spread the max_daily_day value out to every hour with a right-aligned look "back"
max_daily_hour <- roll_max(max_daily_day, width = 24, align = "right", by = 1, na.rm = TRUE)

# Plot with 3- and 24-hr rolling means
layout(matrix(seq(3)))
plot(t, max_24hr, col = 'red')
points(t, x, pch = 16, cex = 0.5)
title("Rolling 24-hr Max")
plot(t, max_daily_day, col = 'red')
points(t, x, pch = 16, cex = 0.5)
title("Daily 24-hr Max")
plot(t, max_daily_hour, col = 'red')
points(t, x, pch = 16, cex = 0.5)
title("Hourly Daily Max")
layout(1)

## ----example_3----------------------------------------------------------------
library(MazamaRollUtils)

# Extract vectors from our example dataset
t <- example_pm25$datetime
x <- example_pm25$pm25

# Create weights for a 9-element exponentially weighted window
#   See:  https://en.wikipedia.org/wiki/Moving_average
N <- 9
alpha <- 2/(N + 1)
w <- (1-alpha)^(0:(N-1))
weights <- rev(w)          # right aligned window

EMA <- roll_mean(x, width = N, align = "right", weights = weights)

# Plot Exponential Moving Average (EMA)
plot(t, x, pch = 16, cex = 0.5)
lines(t, EMA, col = 'red')
title("9-Element Exponential Moving Average")

