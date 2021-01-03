## required for live compilation, perhaps a bug also known in testthat pkg
## https://stackoverflow.com/questions/12410694/rbundler-build-error-cannot-open-file-startup-rs-no-such-file-or-directory
##this should be done in the first line of a testthat.R script
Sys.setenv("R_TESTS" = "")

debug <- FALSE

if ( debug ) {
    setwd("~/programs/ConsSeq/tests")
    source("../R/consseg_r.R")
    library(Rcpp)
    sourceCpp("../src/consseg.cpp")
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
aeh <- function(L) {L^2/2}
## weights
w <- rep(1/M, M)


aehr <- list(negentropy=function(L,n) (L/n)*log(L/n),
             subqu=function(L,n) L^(3/2),
             quad=function(L,n) L^2,
             cub=function(L,n) L^3/3,
             quin=function(L,n) L^5/5,
             expo=function(L,n) exp(L/10)-1) #skip /10
aehc <- list(negentropy="(L/n)*log(L/n)",
             subqu="pow(L,1.5)",
             quad="L*L",
             cub="L*L*L/3",
             quin="L*L*L*L*L/5",
             expo="exp(L/10)-1.0") 
exprs <- list(negentropy=expression(italic(z/n)*log(italic(z/n))),
              subqu=expression(italic(z^(3/2))),
              quad=expression(italic(z^2)),
              cub=expression(italic(z^3/3)),
              quin=expression(italic(z^5/5)),
              expo=expression(exp(italic(z/2))-1))

bps <- rep(0,length(aehr))


png("consseg_potential_scan_r.png", width=2*3.5, height=2.5,
    units="in",res=200)
par(mfrow=c(2,3),mai=c(.35,.05,.05,.05), mgp=c(1.4,.3,0), tcl=-.25)
for ( i in 1:length(aehr) ) {
    cons <- consensus_r(b, n=n, w=w, e=aehr[[i]])

    plot_breaklist(b,axis1=FALSE, axis2=FALSE, col=NA)
    abline(v=cons$breakpoints, col="#0000FFCC", lwd=2)
    plot_breaklist(b,add=TRUE, col=1, lwd=1, length=.05,
                   axis1=FALSE, axis2=FALSE)
    axis(1)
    mtext(exprs[[names(aehr)[i]]],1,1.25*par("mgp")[1], cex=.9)
    cat(paste(paste(cons$breakpoints, collapse=", "), "\n"))
}
dev.off()

## the same, using compiler
## TODO: this does not work properly, perhaps
## the functions need proper casting?
## TODO:  n ignored

fnct <- "long double my_aeh(int l, int n) { long double L=1.0*l;long double N=1.0*n; return("

png("consseg_potential_scan_c.png", width=2*3.5, height=2.5,
    units="in",res=200)
par(mfrow=c(2,3),mai=c(.35,.05,.05,.05), mgp=c(1.4,.3,0), tcl=-.25)
for ( i in 1:length(aehc) ) {
    fnc <- paste(fnct,aehc[[i]],");}")
    cons <- consensus(b, n=n, w=w, e=fnc, verb=0)

    plot_breaklist(b,axis1=FALSE, axis2=FALSE, col=NA)
    abline(v=cons, col="#0000FFCC", lwd=2)
    plot_breaklist(b,add=TRUE, col=1, lwd=1, length=.05,
                   axis1=FALSE, axis2=FALSE)
    axis(1)
    mtext(exprs[[names(aehc)[i]]],1,1.25*par("mgp")[1], cex=.9)
    ## compare R and C output
    conr <- consensus_r(b, n=n, w=w, e=aehr[[i]])
    cat(paste(names(aehc)[i], "in C++ vs R:\n"))
    cat(paste("C:",paste(cons, collapse=", "), "\n"))
    cat(paste("R:",paste(conr$breakpoints, collapse=", "), "\n"))
}
dev.off()
