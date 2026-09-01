# 'min_valid' -- minimum number of non-NA values required within a window

roll_fns <- list(
  roll_mean   = roll_mean,
  roll_median = roll_median,
  roll_max    = roll_max,
  roll_min    = roll_min,
  roll_sum    = roll_sum,
  roll_prod   = roll_prod,
  roll_MAD    = roll_MAD,
  roll_hampel = roll_hampel
)

test_that("min_valid = NULL (default) reproduces the legacy result", {
  x <- c(1, 2, NA, 4, NA, NA, 7, 8, 9, 10)

  for (nm in names(roll_fns)) {
    fn <- roll_fns[[nm]]
    expect_identical(fn(x, 3), fn(x, 3, min_valid = NULL), info = nm)
    expect_identical(
      fn(x, 3, na.rm = TRUE),
      fn(x, 3, na.rm = TRUE, min_valid = NULL),
      info = nm
    )
  }
})

test_that("min_valid = 1 is equivalent to na.rm = TRUE", {
  x <- c(1, 2, NA, 4, NA, NA, 7, 8, 9, 10)

  for (nm in names(roll_fns)) {
    fn <- roll_fns[[nm]]
    expect_equal(fn(x, 3, min_valid = 1), fn(x, 3, na.rm = TRUE), info = nm)
  }
})

test_that("min_valid = width requires a complete window", {
  x <- c(1, 2, NA, 4, NA, NA, 7, 8, 9, 10)

  for (nm in names(roll_fns)) {
    fn <- roll_fns[[nm]]
    expect_equal(fn(x, 3, min_valid = 3), fn(x, 3, na.rm = FALSE), info = nm)
  }
})

test_that("min_valid sets positions with too few valid values to NA", {
  x <- c(1, 2, NA, 4, NA, NA, 7, 8, 9, 10)

  # width-3 centered windows and their non-NA counts:
  #   idx:      1   2   3   4   5   6   7   8   9  10
  #   count:   NA   2   2   1   1   1   2   3   3  NA
  result <- roll_mean(x, 3, min_valid = 2)

  expect_equal(which(!is.na(result)), c(2, 3, 7, 8, 9))
  expect_equal(result[2], mean(c(1, 2)))
  expect_equal(result[7], mean(c(7, 8)))
  expect_equal(result[8], mean(c(7, 8, 9)))
})

test_that("min_valid implies NA-tolerant counting regardless of na.rm", {
  x <- c(1, 2, NA, 4, NA, NA, 7, 8, 9, 10)

  expect_identical(
    roll_mean(x, 3, min_valid = 2),
    roll_mean(x, 3, na.rm = FALSE, min_valid = 2)
  )
  expect_identical(
    roll_mean(x, 3, min_valid = 2),
    roll_mean(x, 3, na.rm = TRUE, min_valid = 2)
  )
})

test_that("min_valid respects align and by", {
  h <- c(1, 2, NA, NA, 5, 6, 7, 8, NA, NA, NA, 12, 13, 14, 15, 16)

  left <- roll_mean(h, 4, by = 4, align = "left", min_valid = 3)
  # window starts: [1:4] 2 valid -> NA, [5:8] 4 valid -> 6.5,
  #                [9:12] 1 valid -> NA, [13:16] 4 valid -> 14.5
  expect_equal(left[c(1, 5, 9, 13)], c(NA, 6.5, NA, 14.5))
  expect_true(all(is.na(left[-c(1, 5, 9, 13)])))

  x <- c(1, 2, NA, 4, NA, NA, 7, 8, 9, 10)
  expect_length(roll_mean(x, 3, align = "right", min_valid = 2), length(x))
})

test_that("min_valid composes with weights in roll_mean", {
  x <- c(1, 2, NA, 4, NA, NA, 7, 8, 9, 10)

  result <- roll_mean(x, 3, weights = c(1, 2, 1), min_valid = 2)

  # position 2: window (1, 2, NA), weights (1, 2) over the two valid values
  expect_equal(result[2], (1 * 1 + 2 * 2) / (1 + 2))
  expect_true(is.na(result[4]))   # only one valid value in the window
})

test_that("min_valid is validated", {
  x <- 1:10

  expect_error(roll_mean(x, 3, min_valid = 0), "single positive integer")
  expect_error(roll_mean(x, 3, min_valid = -1), "single positive integer")
  expect_error(roll_mean(x, 3, min_valid = 2.5), "single positive integer")
  expect_error(roll_mean(x, 3, min_valid = c(1, 2)), "single positive integer")
  expect_error(roll_mean(x, 3, min_valid = NA), "single positive integer")
  expect_error(roll_mean(x, 3, min_valid = "2"), "single positive integer")
  expect_error(roll_mean(x, 3, min_valid = 4), "cannot be larger than 'width'")
})

test_that("roll_sd and roll_var do not accept min_valid", {
  x <- c(1, 2, NA, 4, 5)

  expect_error(roll_sd(x, 3, min_valid = 2))
  expect_error(roll_var(x, 3, min_valid = 2))
})
