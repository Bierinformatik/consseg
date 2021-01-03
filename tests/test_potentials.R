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
n <- 50# 5000 #SEQUENCE LENGTH
M <- 10 # number segmentations (breakpoint lists)
l <- 4# # average number of segments
set.seed(1) # for constant results
b <- random_breakpoints(m=M,n=n,lambda=l)


## potential
aeh <- function(L) {L^2/2}
## weights
w <- rep(1/M, M)


aeh_r <- list(negentropy=function(L,n) (L/n)*log(L/n),
             subqu=function(L,n) L^(3/2),
             quad=function(L,n) L^2,
             cub=function(L,n) L^3/3,
             quin=function(L,n) L^5/5,
             expo=function(L,n) exp(L/10)-1) #skip /10
aeh_c <- list(negentropy="(L/n)*log(L/n)",
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

bps <- rep(0,length(aeh_r))





## COMPARE R AND C OUTPUT
## TODO: start testthat and/or microbenchmark here

for ( i in 1:length(aeh_c) ) {
    conc <- consensus(b, n=n, w=w, e=aeh_c[[i]], verb=0)
    conr <- consensus_r(b, n=n, w=w, e=aeh_r[[i]])$breakpoints
    if ( length(conc)!=length(conr) ) {
        cat(paste("WARNING: difference in R and Rcpp implementations",
                  " at ", i, names(aeh_c)[i], "\n"))
        cat(paste("Fnc:\t", names(aeh_c)[i], "\n"))
        cat(paste("C++:\t",paste(conc, collapse=", "), "\n"))
        cat(paste("R:\t",paste(conr, collapse=", "), "\n"))
    }
}

    
## SCAN OF POTENTIAL FUNCTIONS IN Rcpp
png("test_potentials_c.png", width=2*3.5, height=2.5,
    units="in",res=200)
par(mfrow=c(2,3),mai=c(.35,.05,.05,.05), mgp=c(1.4,.3,0), tcl=-.25)
for ( i in 1:length(aeh_c) ) {
    
    cons <- consensus(b, n=n, w=w, e=aeh_c[[i]], verb=0)

    plot_breaklist(b,axis1=FALSE, axis2=FALSE, col=NA)
    abline(v=cons, col="#0000FFCC", lwd=2)
    plot_breaklist(b,add=TRUE, col=1, lwd=1, length=.05,
                   axis1=FALSE, axis2=FALSE)
    axis(1)
    mtext(exprs[[names(aeh_c)[i]]],1,1.25*par("mgp")[1], cex=.9)
}
dev.off()


## SCAN OF POTENTIAL FUNCTIONS IN R
png("test_potentials_r.png", width=2*3.5, height=2.5,
    units="in",res=200)
par(mfrow=c(2,3),mai=c(.35,.05,.05,.05), mgp=c(1.4,.3,0), tcl=-.25)
for ( i in 1:length(aeh_r) ) {
    cons <- consensus_r(b, n=n, w=w, e=aeh_r[[i]])

    plot_breaklist(b,axis1=FALSE, axis2=FALSE, col=NA)
    abline(v=cons$breakpoints, col="#0000FFCC", lwd=2)
    plot_breaklist(b,add=TRUE, col=1, lwd=1, length=.05,
                   axis1=FALSE, axis2=FALSE)
    axis(1)
    mtext(exprs[[names(aeh_r)[i]]],1,1.25*par("mgp")[1], cex=.9)
    cat(paste(paste(cons$breakpoints, collapse=", "), "\n"))
}
dev.off()

## TEST USER-SUPPLIED PRE-COMPILATION
e <- "(L/n)*log(L/n)" # "L*L*L/3" #
ec <- compileEquation(e)
consensus(b, n=n, w=w, e=e, verb=0)
consensus(b, n=n, w=w, e=ec, verb=0)

## scan over L
res <- rep(NA, 500)
for ( L in 1:length(res) )
    res[L] <- evaluateEquation(e=ec, L=L, n=500)
png("test_potentials_evaluate.png", width=3.5, height=3.5,
    units="in",res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.4,.3,0), tcl=-.25)
plot(1:length(res),res, xlab="L", ylab=e)
dev.off()

