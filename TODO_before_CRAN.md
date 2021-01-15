#NEED TO DO

## Real BUGS?
* aeh function returns L^3/3, shouldn't that be L^2/2 per default?

* check initialization of vectors in R, is it ok that $min_j \Delta[k]$,
and $\delta^*$ are NA until AFTER the first border? See figure
`test_internals.[R|png]`.

## Documentation: README.md, Roxygen and Vignette

* add README.md with github and cran install instructions,
and minimal text + example,
* move plots from current test/ to roxygen doc and/or Vignette.
* sanity- and spell-check of existing roxygen doc.
* reconsider which functions in R and cpp should be public
for the user (consensus_r/c, potential functions for comparison)
* Add missing autors


#EXTRA

## Potential Functions

* add more internal pre-defined potential functions,
* vectorize evaluteEquation and generate matrix over L and n.

## Test/Debug Mode

* implement the Dslow/test calculations in Rcpp, test against
R implementation and use for testing/debugging,
* use the testhat package to write proper tests, based on
the current code in `tests/`,
* `test_internals.R` and `tests/test_compilation.R`: 
   - fuse to systematically test all three implementations:
   R/slow, R/incremental, and Rcpp,
* `test_potentials.R`:
   - compare built-in vs. R vs. Rcpp user-supplied
     potential functions.
* allow to use R version from consensus wrapper for testing and
debugging, in wrapper.

## Optimize Rcpp

* benchmark for R/slow (consensus_r(...,test=TRUE) vs. R/fast vs. Rcpp,
and built-in vs. supplied potential functions,
* use int and long double instead of NumericVector to optimize
memory usage, cast properly in potential function and for debug return


## ADD FUNCTIONALITY

### Alternative Algorithms

* implement symmetrized boundary movers distance,
* implement breakpoint only code, using the conjecture
that each consensus breakpoint is also a breakpoint in
1+ of the input segmentations,

### Utilities

* convert to and from IRanges/GRanges,
* `plot_breaklist`: allow colors, eg.
    - color gradient for weights,
    - colors of input segment classes (segmenTier cluster association),
* `plot_breaklist`: dedicated segmentations vs. consensus plot,
    - indicate origin breakpoints in segmentation.

