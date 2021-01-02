
## EXAMPLE DATASET FROM BUDDING YEAST
## generate ~ Figure 3 from Machne, Murray, Stadler 2017 (segmenTier)
## and calculate and add consensus segmentation

library(segmenTier)
data(primseg436)
source("../R/consseg_r.R")
library(Rcpp)
sourceCpp("../src/consseg.cpp")

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

## CALCULATE CONSENSUS
## NOTE: sset contains sequence length n, and is thus not required as argument
csegs <- consensus(sset, return="breakpoints")


## plot all
plotdev("consseg_data_c",res=300,width=10,height=5,type=fig.type)
layout(matrix(1:4,ncol=1), heights=c(.25,.5,.5,.15))
par(mai=c(0.1,2,0.05,0.1),mgp=c(1.3,.4,0),tcl=-.25, xaxs="i",yaxs="r")
par(cex=1) 
plot(tset,ylabh=TRUE) #, plot="timeseries")
axis(1,at=pretty(c(0,10000)))
par(cex=1.2) # increase axis labels
par(mai=c(0.01,2,0.01,0.1))
plot(sset,"segments",lwd=3)
abline(v=csegs)
axis(1,at=pretty(c(0,10000)),labels=NA)
axis(1,at=pretty(c(0,10000)),labels=NA,tcl=-par("tcl"))
## ADD CONSENSUS BREAKPOINTS
par(mai=c(0.1,2,0.05,0.1),tcl=0)
plot.breaklist(list(consensus=bp2seg(csegs)), n=sset$N, axis1=FALSE)
abline(v=csegs)
par(tcl=-.25)
axis(1,at=pretty(c(0,10000)),labels=NA,tcl=-par("tcl"))
dev.off()


## dump breakpoints
csegs



## convert from segment start/end table to breakpoint list
## add +1 to ends, defining segment starts as breakpoints
n <- sset$N
blst <- split(sset$segments, f=sset$segments$type)
b <- lapply(blst, function(x) c(x$start,x$end+1))
M <- length(b)

## CALCULATE CONSENSUS
bl <- lapply(b, function(x) sort(unique(c(1,x,n,n+1))))
cons <- consensus_c(bl, n=n)
