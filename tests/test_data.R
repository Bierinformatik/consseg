###TEST AREA
library(ConsSeg)
library(IRanges)


### use real data
data(primseg436_sset)

cons <- consensus(sset$segments, w, e)

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
par(mfcol=c(1,1), mai=c(.5,.1,.1,.1))
plot(sset)
axis(1,at=pretty(c(0,10000)))
dev.off()

