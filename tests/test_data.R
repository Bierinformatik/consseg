###TEST AREA
library(ConsSeg)
library(IRanges)


### use real data
data(primseg436_sset)

outfile <- "sset_done.rda"
if ( !file.exists(outfile) ) {
    w <- function(m){return(1/m)}
    e <- function(width){ return(width^2/2)}
    cons <- consensus(sset$segments, w, e)
    save.image(outfile)
} else {
    load(outfile)
}

## bind consensus segments to segmenTier object
nsegs <- data.frame(ID=paste0("consensus_",1:length(cons)),
                    CL=NA,
                    type=rep("consensus",length(cons)),
                    start=start(cons), end=end(cons),
                    fuse=NA,
                    color="#000000")
sset$settings <- rbind(sset$settings,
                       consensus=rep(NA,ncol(sset$segments)))
sset$segments <- rbind(sset$segments,
                       nsegs)
sset$ids <- c(sset$ids,"consensus")

png("data_consensus.png", units="in", res=200, width=2*3.5, height=3.5/2)
par(mfcol=c(1,1), mai=c(.5,2,.1,.1))
plot(sset)
axis(1,at=pretty(c(0,10000)))
dev.off()

## TODO: use segment_data.R to re-create Figure 3 from paper,
##       calculate consensus, and re-create Figure 3+consensus.
## library(segmenTier)
##plotdev("segment_data_examples",res=300,width=10,height=5,type=fig.type)
##layout(matrix(1:10,ncol=1),heights=c(.25,.5,.5,.075,.075,.075,.075,.075,.075,.075))
##par(mai=c(0.1,2,0.05,0.01),xaxs="i",yaxs="r")
##par(cex=1) 
##plot(tset,ylabh=TRUE)
##par(cex=.6) 
##plot(cset,axes=2,cex=.7,ylabh=FALSE); mtext("clustering",2,7.2,las=2,cex=1.2)
##par(cex=1.2) # increase axis labels
##par(mai=c(0.01,2,0.01,0.01))
##plot(sset,"segments",lwd=3)
##dev.off()

