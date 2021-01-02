
debug <- FALSE

if ( debug ) {
    setwd("~/programs/ConsSeq/tests")
    source("../R/consseg_r.R")
} else {
    library(ConsSeg)
}


## GENERATE RANDOM SEGMENTS
n <- 50 # SEQUENCE LENGTH
M <- 10 # number segmentations (breakpoint lists)
l <- 4 # average number of segments
set.seed(1) # for constant results
b <- random_breakpoints(m=M,n=n,lambda=l)

## potential
aeh <- function(L,n) {L^2/2}
## weights
w <- rep(1/M, M)

cons <- consensus_r(b, n=n, w=w, e=aeh, store=TRUE, test=FALSE)

## plot segments, leave room for consensus arrows
png("consseg_test.png", units="in", width=3.5, height=7, res=200)
par(mfcol=c(6,1),mai=c(.5,.5,.1,.1), mgp=c(1.4,.3,0), tcl=-.25)
par(mai=c(.1,.5,.1,.1))
plot_breaklist(b,lwd=1)
abline(v=cons$breakpoints, col="#0000FF55", lwd=2)
plot_breaklist(b,add=TRUE, col=1, lwd=1)
plot(cons$F,  type="p",xlab="sequence position k", ylab="F(k)")
plot(cons$Dk, type="p",xlab="sequence position k",
     ylab=expression(min[j]~Delta(k)))
plot(cons$dstar, ylab=expression(delta*"*"))
plot(cons$dcd, ylab=expression(delta^intersect()),
     ylim=range(c(cons$dcd,cons$dcu),na.rm=TRUE))
points(cons$dcu, col=2)
legend("bottom",c(expression(delta["<"]^intersect()*(k)),
                  expression(delta[">"]^intersect()*(j+1))),
       bty="n",col=1:2,pch=1)
plot(cons$dsq, ylab=expression(delta),ylim=c(0,max(c(cons$dsq,cons$dsm))))
points(cons$dsm, col=2)
legend("bottomright",c(expression(delta[""<=""](j)),
                       expression(delta[""<""](k))),bty="n",col=1:2,pch=1)
dev.off()

## print breakpoints
cons$breakpoints
