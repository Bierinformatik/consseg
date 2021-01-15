library(consseg)

# Install package from CRAN only if not installed, needed for devtools:check()
if (!require(testthat)) install.packages('testthat', repos = "http://cran.us.r-project.org")
if (!require(magrittr)) install.packages('magrittr', repos = "http://cran.us.r-project.org")

library(testthat)

#testthat::use_catch() used for rcpp testing

test_check("consseg")
