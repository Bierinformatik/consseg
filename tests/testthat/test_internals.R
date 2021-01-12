## TODO: use testthat to test equality here
## test_that throws failures due to NA in d_k and d*


# Test weight function
test_that("Weight normalization",{
    w <- (1,.01)
})

# Test weight function
# TODO: User evaluateEquation to evaluate compileEquation output
#test_that("Compile equation",{
#    expect_match(compileEquation("L*L/2"), "Weight vector does not sum up to 1, normalizing.")
#)}


# Compare breakpoints of Rcpp and base R implementation
test_that("Breakpoints of consensus in Rcpp implementation are equal to base R implementation",{
    expect_equal(cons_c$breakpoints, cons_r$breakpoints)
})


# Compare backtracing pointers between implementations
test_that("Backtracing pointer of consensus in Rcpp implementation are equal to base R implementation",{
    expect_equal(cons_c$values$ptr, cons_r$values$ptr)
})


# Compare Deltas between implementation
test_that("Delta_k of consensus in Rcpp implementation are equal to base R implementation",{
    expect_equal(cons_c$values$Dk, cons_r$values$Dk)
})


# Compare Deltas between implementation
test_that("Delta_star of consensus in Rcpp implementation are equal to base R implementation",{
    expect_equal(cons_c$values$dstar, cons_r$values$dstar)
})
