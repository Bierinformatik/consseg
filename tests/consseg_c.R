
setwd("~/programs/ConsSeq/tests")

source("../R/consseg_r.R")
library(IRanges)

library(Rcpp)
sourceCpp("../src/consseg.cpp")

## GENERATE RANDOM SEGMENTS
n <- 50 # SEQUENCE LENGTH
M <- 10 # number segmentations (breakpoint lists)
avg.sg <- n/10 # average number of segments

set.seed(1) # for constant results
segs <- simulate_ranges_r(n,M,avg.sg,FALSE, df=TRUE)


### CONVERT TO LIST OF BREAKPOINTS
## fill up -1 ends (only required for plots)
##segs$end <- segs$end+1
## each starting with 1 and ending with n+1 (for convenience)
blst <- split(segs, f=segs$type)
b <- lapply(blst, function(x) unique(x$start))

## potential
aeh <- function(L) {L^2/2}
## weight
w <- function(m){return(1/m)}
## override as in Fall's code
## ignore m and take global M
w <- function(m) return(1/M) 

cons <- consensus_r(b, n=n, w=w, aeh=aeh, store.matrix=TRUE, test.slow=FALSE)

bl <- lapply(b, function(x) sort(unique(c(1,x,n,n+1))))
conc <- consensus_c(bl, n=n)


## plot segments, leave room for consensus arrows
png("consseg_c.png", units="in", width=3.5, height=7, res=200)
par(mfcol=c(6,1),mai=c(.5,.5,.1,.1), mgp=c(1.4,.3,0), tcl=-.25)
par(mai=c(.1,.5,.1,.1))
plot.breaklist(blst,lwd=1)
abline(v=cons$breakpoints, col="#0000FF55", lwd=2)
abline(v=conc$breakpoints, col="#FF000077",lwd=2)
plot.breaklist(blst,add=TRUE, col=1, lwd=1)
plot(cons$F,  type="p",xlab="sequence position k", ylab="F(k)")
lines(conc$F, col=2)
plot(cons$Dk, type="p",xlab="sequence position k",
     ylab=expression(min[j]~Delta(k)))
plot(cons$dstar, ylab=expression(delta*"*"))
plot(cons$dcd, ylab=expression(delta^intersect()),
     ylim=range(c(cons$dcd,cons$dcu),na.rm=TRUE))
points(cons$dcu, col=4)
lines(conc$dcd, col=2)
lines(conc$dcu, col=2)
legend("bottom",c(expression(delta["<"]^intersect()*(k)),
                  expression(delta[">"]^intersect()*(j+1))),
       bty="n",col=1:2,pch=1)
plot(cons$dsq, ylab=expression(delta),ylim=c(0,max(c(cons$dsq,cons$dsm))))
points(cons$dsm, col=2)
lines(conc$dsm, col=2)
lines(conc$dsq, col=2)
legend("bottomright",c(expression(delta[""<=""](j)),
                       expression(delta[""<""](k))),bty="n",col=1:2,pch=1)
dev.off()

## print breakpoints
cons$breakpoints
