/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */

// All test files should include the <testthat.h>
// header file.
#include "testthat.h"
#include "consseg_types.h"
#include "Rcpp.h"


//Rcpp bindings
using namespace Rcpp;

// Some definitions for testing
long double aeh(double L, double n) {
  //long double e =  L*1.0; // implicit cast, required if using int
  return(L*L*L/3);
}

double L = 100;
double n = 5;

NumericVector res(1);
aeh(L, n);

//evaluateEquation(aeh);

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
context("Internals") {

  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("Equation") {
    expect_true(  e == res );
  }

}
