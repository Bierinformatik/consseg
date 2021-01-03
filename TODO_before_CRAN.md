# Real BUGS?

* check back-tracking, are 1,2 breaktpoints in `consseg_potential_scan.R?
a numerical artefact of L^5 or a bug?

# Test/Debug Mode

* return internal values only with option store or remove altogether
* allow to use R version from consensus wrapper for testing and
debugging, in wrapper, allow to use both and R version with test=TRUE,
to compare all to the slow direct implementation

# Reconsider public functions

* reconsider which functions in R and cpp should be public
for the user (consensus_r/c, potential functions for comparison)

# Optimize Rcpp

* use int and long double instead of NumericVector to optimize
memory usage, cast properly in potential function and for debug return

# Potential Functions

* add more internal pre-defined potential functions,
* vectorize evaluteEquation and generate matrix over L and n

# Examples and Vignette

* move plots from current test/ to roxygen docu and/or Vignette,
* use the testhat package to write proper tests, based on
the same code as vignette?


