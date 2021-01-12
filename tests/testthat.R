## required for live compilation, perhaps a bug also known in testthat pkg
## https://stackoverflow.com/questions/12410694/rbundler-build-error-cannot-open-file-startup-rs-no-such-file-or-directory

library(consseg)
library(testthat)

## GENERATE RANDOM SEGMENTS
n <- 50 # SEQUENCE LENGTH
M <- 10 # number segmentations (breakpoint lists)
l <- 4 # average number of segments
set.seed(1) # for constant results
b <- random_breakpoints(m=M,n=n,lambda=l)
bl <- lapply(b, function(x) sort(unique(c(1,x,n,n+1))))

## weights
w <- rep(1/M, M)

## R implementation
aeh_r <- function(L,n) {L^5/5}
cons_r <- consensus_r(bl, n=n, w=w, e=aeh_r, store=TRUE, test=TRUE)

## Rcpp implementation
aeh_c <- compileEquation("L*L*L*L*L/5")
cons_c <- consensus_c(bl, w=w, e=aeh_c, n=n,store=TRUE)

test("consseg")
test_check("consseg")
