# Real BUGS?

* clarify why L*L*L/3 has novel breakpoints in consensus
of data example (`tests/consseg_data.R`)
     - similar to `tests/test_potential_scan.R` one should
       meticulously compare several potential functions in
       R, Rcpp via compile, and built-in (switch in e_ptr
       function).
* check back-tracing, are 1,2 breaktpoints in `consseg_potential_scan.R?
and artefact?
* check and rework all datatypes in Rcpp.

# Test/Debug Mode

* return internal values only with option store

# Optimize Rcpp

* use long double instead of NumericVector, and only
cast for debug return

# Potential Functions

* add more internal pre-defined potential functions,

# Examples and Vignette

* move plots from current test/ to roxygen docu and/or Vignette,
* use the testhat package to write proper tests, based on
the same code as vignette?


