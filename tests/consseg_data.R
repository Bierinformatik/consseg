
## EXAMPLE DATASET FROM BUDDING YEAST
## generate ~ Figure 3 from Machne, Murray, Stadler 2017 (segmenTier)
## and calculate and add consensus segmentation

library(segmenTier)
data(primseg436)
source("../R/consseg_r.R")

fig.type <- "png"

## pre-process time series for clustering
tset <- processTimeseries(ts=tsd, na2zero=TRUE,
                          trafo="raw", dc.trafo="ash",
                          use.fft=TRUE, dft.range=1:7,
                          use.snr=TRUE, low.thresh=-Inf)

## cluster pre-processed time series
set.seed(15) # stable kmeans clustering
cset <- clusterTimeseries(tset, K=12, nui.thresh=0.6 ) # cluster number 12

### CALCULATE SEGMENTS
vary <- setVarySettings(
    E=c(1,3),    # scale exponent of similarity matrices csim
    S="icor",    # scoring function
    M=c(100,200),# scoring function minimal length penalty
    Mn=100,      # M for nuisance clusters
    nui=c(1,3)   #-/+ correlation of nuisance cluster with others and itself
)

### for chosen segmentation parameters
sset <- segmentCluster.batch(cset, varySettings=vary,
                             id="mysegments",
                             type.name=c("E","M","nui"), # segment type names
                             verb=1, save.matrix=FALSE) 



## convert from segment start/end table to breakpoint list
## add +1 to ends, defining segment starts as breakpoints
n <- sset$N
blst <- split(sset$segments, f=sset$segments$type)
b <- lapply(blst, function(x) c(x$start,x$end+1))
M <- length(b)

## CALCULATE CONSENSUS
w <- function(m){return(1/M)}
e <- function(L){ return(L^2/2)}
cons <- consensus_r(b, n=n, w=w, aeh=e, test.slow=FALSE) 


## convert to segment list as used in plot.breaklist
csegs <- list(consensus=bp2seg(cons$breakpoints))

## plot all
plotdev("consseg_data",res=300,width=10,height=5,type=fig.type)
layout(matrix(1:4,ncol=1), heights=c(.25,.5,.5,.15))
par(mai=c(0.1,2,0.05,0.1),mgp=c(1.3,.4,0),tcl=-.25, xaxs="i",yaxs="r")
par(cex=1) 
plot(tset,ylabh=TRUE) #, plot="timeseries")
axis(1,at=pretty(c(0,10000)))
par(cex=1.2) # increase axis labels
par(mai=c(0.01,2,0.01,0.1))
plot(sset,"segments",lwd=3)
abline(v=cons$breakpoints)
axis(1,at=pretty(c(0,10000)),labels=NA)
axis(1,at=pretty(c(0,10000)),labels=NA,tcl=-par("tcl"))
## ADD CONSENSUS BREAKPOINTS
par(mai=c(0.1,2,0.05,0.1),tcl=0)
plot.breaklist(csegs, n=n, axis1=FALSE)
abline(v=cons$breakpoints)
par(tcl=-.25)
axis(1,at=pretty(c(0,10000)),labels=NA,tcl=-par("tcl"))
dev.off()


## dump breakpoints
cons$breakpoints
