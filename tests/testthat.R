## required for live compilation, perhaps a bug also known in testthat pkg
## https://stackoverflow.com/questions/12410694/rbundler-build-error-cannot-open-file-startup-rs-no-such-file-or-directory

library(consseg)
library(testthat)
#testthat::use_catch() used for rcpp testing

test("consseg")
test_check("consseg")
