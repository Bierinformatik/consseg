## required for live compilation, perhaps a bug also known in testthat pkg
## https://stackoverflow.com/questions/12410694/rbundler-build-error-cannot-open-file-startup-rs-no-such-file-or-directory
##this should be done in the first line of a testthat.R script
Sys.setenv("R_TESTS" = "")

debug <- FALSE

if ( debug ) {
    setwd("~/programs/ConsSeq/tests")
    source("../R/consseg_r.R")
    library(Rcpp)
    sourceCpp("../src/consseg.cpp")
} else {
    library(ConsSeg)
}


## GENERATE RANDOM SEGMENTS
n <- 50 # SEQUENCE LENGTH
M <- 10 # number segmentations (breakpoint lists)
l <- 4 # average number of segments
set.seed(1) # for constant results
b <- random_breakpoints(m=M,n=n,lambda=l)
bl <- lapply(b, function(x) sort(unique(c(1,x,n,n+1))))

## potential
aeh <- function(L,n) {L^2/2}
## weights
w <- rep(1/M, M)

cons_r <- consensus_r(bl, n=n, w=w, e=aeh, store=TRUE, test=FALSE)

## plot segments, leave room for consensus arrows
png("consseg_test.png", units="in", width=3.5, height=7, res=200)
par(mfcol=c(6,1),mai=c(.5,.5,.1,.1), mgp=c(1.4,.3,0), tcl=-.25)
par(mai=c(.1,.5,.1,.1))
plot_breaklist(b,lwd=1)
abline(v=cons_r$breakpoints, col="#0000FF55", lwd=2)
plot_breaklist(b,add=TRUE, col=1, lwd=1)
plot(cons_r$values$F,  type="p",xlab="sequence position k", ylab="F(k)")
plot(cons_r$values$Dk, type="p",xlab="sequence position k",
     ylab=expression(min[j]~Delta(k)))
plot(cons_r$values$dstar, ylab=expression(delta*"*"))
plot(cons_r$values$dcd, ylab=expression(delta^intersect()),
     ylim=range(c(cons_r$values$dcd,cons_r$values$dcu),na.rm=TRUE))
points(cons_r$values$dcu, col=2)
legend("bottom",c(expression(delta["<"]^intersect()*(k)),
                  expression(delta[">"]^intersect()*(j+1))),
       bty="n",col=1:2,pch=1)
plot(cons_r$values$dsq, ylab=expression(delta),ylim=c(0,max(c(cons_r$values$dsq,cons_r$values$dsm))))
points(cons_r$values$dsm, col=2)
legend("bottomright",c(expression(delta[""<=""](j)),
                       expression(delta[""<""](k))),bty="n",col=1:2,pch=1)
dev.off()

## print breakpoints
cons_r$breakpoints

ec <- compileEquation("L*L/2")
cons_c <- consensus_c(bl, w=w, e=ec, n=n,store=TRUE)

## TODO: use testthat to test equality here
cons_c$breakpoints
cons_r$breakpoints

cons_c$values$Dk
cons_r$values$Dk

cons_c$values$dstar
cons_r$values$dstar ## NOTE: should be null for test=TRUE option above!
consensus_r(bl, n=n, w=w, e=aeh, store=TRUE, test=FALSE)$values$dstar
