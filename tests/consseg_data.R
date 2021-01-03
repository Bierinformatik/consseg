## required for live compilation, perhaps a bug also known in testthat pkg
## https://stackoverflow.com/questions/12410694/rbundler-build-error-cannot-open-file-startup-rs-no-such-file-or-directory
##this should be done in the first line of a testthat.R script
Sys.setenv("R_TESTS" = "")

## EXAMPLE DATASET FROM BUDDING YEAST
## from Machne, Murray, Stadler 2017 (segmenTier)
## and calculate and add consensus segmentation

debug <- FALSE

if ( debug ) {
    setwd("~/programs/ConsSeq/tests")
    source("../R/consseg_r.R")
    library(Rcpp)
    sourceCpp("../src/consseg.cpp")
    load("../data/primseg436_sset.rda")
} else {
    library(ConsSeg)
    data(primseg436_sset)
}



## CALCULATE CONSENSUS with custom potential function and weights
## to favor breakpoint before AIM9
## NOTE: sset contains sequence length n, and is thus not required as argument
e <- compileEquation("(L/n)*log(L/n)")
csegs <- consensus(sset, w=c(1,1.01,1,1), e=e, return="breakpoints")

## plot all
png("consseg_data.png",width=7,height=3.5/2,res=300,units="in")
layout(matrix(1:2,ncol=1), heights=c(.5,.2))
par(mai=c(0.2,2,0.05,0.1),mgp=c(1.3,.4,0),tcl=-.25, xaxs="i",yaxs="r")
plot_breaklist(sset, n=n, axis1=FALSE)
abline(v=csegs)
axis(1,at=pretty(c(0,10000)))
axis(1,at=pretty(c(0,10000)),labels=NA,tcl=-par("tcl"))
## ADD CONSENSUS BREAKPOINTS
par(mai=c(0.1,2,0.05,0.1),tcl=0)
plot_breaklist(csegs, n=sset$N, axis1=FALSE)
abline(v=csegs)
par(tcl=-.25)
axis(1,at=pretty(c(0,10000)),labels=NA,tcl=-par("tcl"))
dev.off()





