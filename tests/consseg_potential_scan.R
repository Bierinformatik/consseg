## required for live compilation, perhaps a bug also known in testthat pkg
## https://stackoverflow.com/questions/12410694/rbundler-build-error-cannot-open-file-startup-rs-no-such-file-or-directory
##this should be done in the first line of a testthat.R script
Sys.setenv("R_TESTS" = "")

debug <- FALSE

if ( debug ) {
    setwd("~/programs/ConsSeq/tests")
    source("../R/consseg_r.R")
} else {
    library(ConsSeg)
}
library(IRanges)


## GENERATE RANDOM SEGMENTS
n <- 50 # SEQUENCE LENGTH
M <- 10 # number segmentations (breakpoint lists)
avg.sg <- n/10 # average number of segments

set.seed(1) # for constant results
segs <- simulate_ranges_r(n,M,avg.sg,FALSE, df=TRUE)

## shorter names for plot
segs$type <- sub("segmentation","",segs$type)

### CONVERT SEGMENTS TO LIST OF BREAKPOINTS
## split into different types
blst <- split(segs, f=segs$type)
blst <- blst[order(as.numeric(names(blst)))] # order segmentation numbers
## each starting with 1 and ending with n+1 (for convenience)
b <- lapply(blst, function(x) unique(x$start))

## potential
aeh <- function(L) {L^2/2}
## weights
w <- rep(1/M, M)


aehs <- list(negentropy=function(L,n) (L/n)*log(L/n),
             subqu=function(L,n) L^(3/2),
             quad=function(L,n) L^2,
             cub=function(L,n) L^3/3,
             quin=function(L,n) L^5/5,
             expo=function(L,n) exp(L/2.0)-1.0) #skip /10
exprs <- list(negentropy=expression(italic(z/n)*log(italic(z/n))),
              subqu=expression(italic(z^(3/2))),
              quad=expression(italic(z^2)),
              cub=expression(italic(z^3/3)),
              quin=expression(italic(z^5/5)),
              expo=expression(exp(italic(z/2))-1))

bps <- rep(0,length(aehs))

## for clearer plots increase segment ends by +1
blst <- lapply(blst, function(x) {
    x$end[x$end<n]=x$end[x$end<n]+1;return(x)})

png("consseg_potential_scan.png", width=2*3.5, height=2.5,
    units="in",res=200)
par(mfrow=c(2,3),mai=c(.35,.05,.05,.05), mgp=c(1.4,.3,0), tcl=-.25)
for ( i in 1:length(aehs) ) {
    cons <- consensus_r(b, n=n, w=w, e=aehs[[i]])

    plot_breaklist(blst,axis1=FALSE, axis2=FALSE, col=NA)
    abline(v=cons$breakpoints, col="#0000FFCC", lwd=2)
    plot_breaklist(blst,add=TRUE, col=1, lwd=1, length=.05,
                   axis1=FALSE, axis2=FALSE)
    axis(1)
    mtext(exprs[[names(aehs)[i]]],1,1.25*par("mgp")[1], cex=.9)
    cat(paste(paste(cons$breakpoints, collapse=", "), "\n"))
}
dev.off()

## the same, using compiler

aehs <- list(negentropy="(L/50.0)*log(L/50.0)",
             subqu="pow(L,3.0/2.0)",
             quad="L*L",
             cub="L*L*L/3.0",
             quin="pow(L,5.0)/5.0",
             expo="exp(L/2.0)-1.0") 
exprs <- list(negentropy=expression(italic(z/n)*log(italic(z/n))),
              subqu=expression(italic(z^(3/2))),
              quad=expression(italic(z^2)),
              cub=expression(italic(z^3/3)),
              quin=expression(italic(z^5/5)),
              expo=expression(exp(italic(z/2))-1))

png("consseg_potential_scan_c.png", width=2*3.5, height=2.5,
    units="in",res=200)
par(mfrow=c(2,3),mai=c(.35,.05,.05,.05), mgp=c(1.4,.3,0), tcl=-.25)
for ( i in 1:length(aehs) ) {
    fnc <- paste("long double my_aeh(int L, int n) {return ",aehs[[i]],";}")
    cons <- consensus(b, n=n, w=w, e=fnc)

    plot_breaklist(blst,axis1=FALSE, axis2=FALSE, col=NA)
    abline(v=cons, col="#0000FFCC", lwd=2)
    plot_breaklist(blst,add=TRUE, col=1, lwd=1, length=.05,
                   axis1=FALSE, axis2=FALSE)
    axis(1)
    mtext(exprs[[names(aehs)[i]]],1,1.25*par("mgp")[1], cex=.9)
    cat(paste(paste(cons, collapse=", "), "\n"))
}
dev.off()
