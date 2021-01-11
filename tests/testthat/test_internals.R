Sys.setenv("R_TESTS" = "")

debug <- FALSE

if ( debug ) {
    setwd("~/programs/consseg/tests")
    source("../R/consseg_r.R")
    library(Rcpp)
    sourceCpp("../src/consseg.cpp")
} else {
    library(consseg)
}


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

## plot results
png("test_internals.png", units="in", width=3.5, height=7, res=200)
par(mfcol=c(7,1),mai=c(.5,.5,.1,.1), mgp=c(1.4,.3,0), tcl=-.25)
par(mai=c(.1,.5,.1,.1))
plot_breaklist(b,lwd=1)
abline(v=cons_c$breakpoints, col="#FF000099",lwd=2)
abline(v=cons_r$breakpoints, col="#000000", lwd=1)
plot(cons_r$values$F,  type="p",xlab="sequence position k", ylab="F(k)")
lines(cons_c$values$F, col=2)
legend("topright",c("R/slow","Rcpp"),lty=c(NA,1),pch=c(1,NA),col=1:2)
plot(cons_r$values$ptr, type="p",
     ylab=expression("pointer[k]"))
lines(cons_c$values$ptr, col=2)
plot(cons_r$values$Dk, type="p",xlab="sequence position k",
     ylab=expression(min[j]~Delta(k)))
lines(cons_c$values$Dk, col=2)
plot(cons_r$values$dstar, ylab=expression(delta*"*"))
lines(cons_c$values$dstar, col=2)
plot(cons_r$values$dcd, ylab=expression(delta^intersect()),
     ylim=range(c(cons_r$values$dcd,cons_r$values$dcu),na.rm=TRUE))
points(cons_r$values$dcu, col=4)
lines(cons_c$values$dcd, col=2)
lines(cons_c$values$dcu, col=2)
legend("bottom",c(expression(delta["<"]^intersect()*(k)),
                  expression(delta[">"]^intersect()*(j+1))),
       bty="n",col=c(1,4),pch=1)
plot(cons_r$values$dsq, ylab=expression(delta),ylim=c(0,max(c(cons_r$values$dsq,cons_r$values$dsm))))
    points(cons_r$values$dsm, col=4)
lines(cons_c$values$dsm, col=2)
lines(cons_c$values$dsq, col=2)
legend("bottomright",c(expression(delta[""<=""](j)),
                       expression(delta[""<""](k))),bty="n",col=c(1,4),pch=1)
dev.off()


## TODO: use testthat to test equality here
cons_c$breakpoints
cons_r$breakpoints

cons_c$values$ptr
cons_r$values$ptr



cons_c$values$Dk
cons_r$values$Dk

cons_c$values$dstar
cons_r$values$dstar ## NOTE: should be null for test=TRUE option above!
consensus_r(bl, n=n, w=w, e=aeh_r, store=TRUE, test=FALSE)$values$dstar
