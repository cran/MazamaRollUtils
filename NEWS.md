# MazamaRollUtils 1.1.0

* Added a `min_valid` argument to `roll_mean()`, `roll_median()`, `roll_max()`, 
`roll_min()`, `roll_sum()`, `roll_prod()`, `roll_MAD()` and `roll_hampel()`. A 
window with fewer than `min_valid` non-`NA` values returns `NA`; supplying 
`min_valid` implies `NA`-tolerant counting within each window. The default, 
`NULL`, preserves the previous behavior.
* `roll_nowcast()` now returns `0` rather than `NA` for a window whose values 
are all zero.
* `width` or `by` larger than `length(x)`, and negative `weights`, are now 
rejected in R with clearer error messages.
* Clarified that `roll_MAD()` returns the unscaled median absolute deviation, 
equivalent to `stats::mad(x, constant = 1)`.
* Documentation updates: rewrote the README, refreshed the vignette, and 
expanded function-level help.

# MazamaRollUtils 1.0.0

* Review/refactor with minor bug fixes.
* Functions now match the behavior of base R equivalents when `align = "center"` 
and `width` is an even number.
* Added `roll_nowcast()`.
* More extensive test suite.

# MazamaRollUtils 0.1.4

* Address CRAN testing error.

# MazamaRollUtils 0.1.3

* Release candidate for testing.

# MazamaRollUtils 0.1.2

* Work toward release candidate.

# MazamaRollUtils 0.0.1

* Initial setup.
