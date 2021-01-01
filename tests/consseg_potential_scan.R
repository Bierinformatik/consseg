

library(segmenTier) #ONLY FOR plotdev !!
source("../R/consseg_r.R")

fig.type <- "png"

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


aehs <- list(negentropy=function(L,n=50) (L/n)*log(L/n),
             subqu=function(L) L^(3/2),
             quad=function(L) L^2,
             cub=function(L) L^3/3,
             quin=function(L) L^5/5,
             expo=function(L) exp(L/2)-1) #skip /10
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

plotdev("consseg_potential_scan", width=2*3.5, height=2.5,
                    type=fig.type)
par(mfrow=c(2,3),mai=c(.35,.05,.05,.05), mgp=c(1.4,.3,0), tcl=-.25)
for ( i in 1:length(aehs) ) {
    cons <- consensus_r(b, n=n, w=w, aeh=aehs[[i]])

    plot.breaklist(blst,axis1=FALSE, axis2=FALSE, col=NA)
    abline(v=cons$breakpoints, col="#0000FFCC", lwd=2)
    plot.breaklist(blst,add=TRUE, col=1, lwd=1, length=.05,
                   axis1=FALSE, axis2=FALSE)
    axis(1)
    mtext(exprs[[names(aehs)[i]]],1,1.25*par("mgp")[1], cex=.9)
    cat(paste(paste(cons$breakpoints, collapse=", "), "\n"))
}
dev.off()


