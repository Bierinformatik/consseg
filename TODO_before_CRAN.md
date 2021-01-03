# Real BUGS?

* check casting requirements in user-defined function,
* check back-tracking, are 1,2 breaktpoints in `consseg_potential_scan.R?
a numerical artefact of L^5 or a bug?

# Test/Debug Mode

* return internal values only with option store or remove altogether

# Optimize Rcpp

* use int and long double instead of NumericVector, and only
cast for debug return

# Potential Functions

* add more internal pre-defined potential functions,
* allow user to parse only math, and generate Rcpp function with
proper type casts and function signature.

# Examples and Vignette

* move plots from current test/ to roxygen docu and/or Vignette,
* use the testhat package to write proper tests, based on
the same code as vignette?


