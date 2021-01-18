Sys.setenv("R_TESTS" = "")

library(testthat)
skip("This only for manual testing")

debug <- FALSE

if ( debug ) {
    setwd("~/programs/consseg/tests")
    source("../R/consseg_r.R")
    library(Rcpp)
    sourceCpp("../src/consseg.cpp")
    load("../data/primseg436_sset.rda")
} else {
    library(consseg)
    data(primseg436_sset)
}

## GENERATE RANDOM SEGMENTS
n <- 50 # SEQUENCE LENGTH
M <- 10 # number segmentations (breakpoint lists)
l <- 4 # average number of segments
set.seed(1) # for constant results across tests
b <- random_breakpoints(m=M,n=n,lambda=l)
bl <- lapply(b, function(x) sort(unique(c(1,x,n,n+1))))


#### TEST INTERNALS ####
## potential
aeh_r <- function(L) {L^2/2}
aeh_c <- compileEquation("L*L/2")

## weights
w <- rep(1/M, M)

## R implementation
cons_r <- consensus_r(bl, n=n, w=w, e=aeh_r, store=TRUE, test=TRUE)

## Rcpp implementation
cons_c <- consensus_c(bl, w=w, e=aeh_c, n=n,store=TRUE)


## plot results
png("test_internals.png", units="in", width=3.5, height=7, res=200)
par(mfcol=c(7,1),mai=c(.5,.5,.1,.1), mgp=c(1.4,.3,0), tcl=-.25)
par(mai=c(.1,.5,.1,.1))
plot_breaklist(b,lwd=1)
abline(v=cons_c$breakpoints, col="#FF000099",lwd=2)
abline(v=cons_r$breakpoints, col="#000000", lwd=1)
plot(cons_r$values$F,  type="p",xlab="sequence position k", ylab="F(k)")
lines(cons_c$values$F, col=2)
legend("topright",c("R/slow","Rcpp"),lty=c(NA,1),pch=c(1,NA),col=1:2)
plot(cons_r$values$ptr, type="p",
     ylab=expression("pointer[k]"))
lines(cons_c$values$ptr, col=2)
plot(cons_r$values$Dk, type="p",xlab="sequence position k",
     ylab=expression(min[j]~Delta(k)))
lines(cons_c$values$Dk, col=2)
plot(cons_r$values$dstar, ylab=expression(delta*"*"))
lines(cons_c$values$dstar, col=2)
plot(cons_r$values$dcd, ylab=expression(delta^intersect()),
     ylim=range(c(cons_r$values$dcd,cons_r$values$dcu),na.rm=TRUE))
points(cons_r$values$dcu, col=4)
lines(cons_c$values$dcd, col=2)
lines(cons_c$values$dcu, col=2)
legend("bottom",c(expression(delta["<"]^intersect()*(k)),
                  expression(delta[">"]^intersect()*(j+1))),
       bty="n",col=c(1,4),pch=1)
plot(cons_r$values$dsq, ylab=expression(delta),ylim=c(0,max(c(cons_r$values$dsq,cons_r$values$dsm))))
points(cons_r$values$dsm, col=4)
lines(cons_c$values$dsm, col=2)
lines(cons_c$values$dsq, col=2)
legend("bottomright",c(expression(delta[""<=""](j)),
                       expression(delta[""<""](k))),bty="n",col=c(1,4),pch=1)
dev.off()

#### TEST COMPILATION ####
## SLOW R IMPLEMENTATION
## TODO: compare incremental and slow R implementations
## TODO: relative error in test function, and count deviations
## potential


## print breakpoints
cons_r$breakpoints


### TEST MAIN WRAPPER
consensus(b,n=50,w=w)

## test compiling potential function
e <- "L*L*L/3"
consensus(b, n=50, w=w, e=e)

## test pre-compiled potential function
ec <- compileEquation(e)
consensus(b, n=50, w=w, e=ec)

## check with direct calculation above
cons_c$breakpoints
cons_r$breakpoints

#Plot compilation comparison
## plot results
png("test_compilation.png", units="in", width=3.5, height=7, res=200)
par(mfcol=c(7,1),mai=c(.5,.5,.1,.1), mgp=c(1.4,.3,0), tcl=-.25)
par(mai=c(.1,.5,.1,.1))
plot_breaklist(b,lwd=1)
abline(v=cons_c$breakpoints, col="#FF000099",lwd=2)
abline(v=cons_r$breakpoints, col="#000000", lwd=1)
plot(cons_r$values$F,  type="p",xlab="sequence position k", ylab="F(k)")
lines(cons_c$values$F, col=2)
legend("topright",c("R/fast","Rcpp"),lty=c(NA,1),pch=c(1,NA),col=1:2)
plot(cons_r$values$ptr, type="p",
     ylab=expression("pointer[k]"))
lines(cons_c$values$ptr, col=2)
plot(cons_r$values$Dk, type="p",xlab="sequence position k",
     ylab=expression(min[j]~Delta(k)))
lines(cons_c$values$Dk, col=2)
plot(cons_r$values$dstar, ylab=expression(delta*"*"))
lines(cons_c$values$dstar, col=2)
plot(cons_r$values$dcd, ylab=expression(delta^intersect()),
     ylim=range(c(cons_r$values$dcd,cons_r$values$dcu),na.rm=TRUE))
points(cons_r$values$dcu, col=4)
lines(cons_c$values$dcd, col=2)
lines(cons_c$values$dcu, col=2)
legend("bottom",c(expression(delta["<"]^intersect()*(k)),
                  expression(delta[">"]^intersect()*(j+1))),
       bty="n",col=c(1:4),pch=1)
plot(cons_r$values$dsq, ylab=expression(delta),ylim=c(0,max(c(cons_r$values$dsq,cons_r$values$dsm))))
points(cons_r$values$dsm, col=4)
lines(cons_c$values$dsm, col=2)
lines(cons_c$values$dsq, col=2)
legend("bottomright",c(expression(delta[""<=""](j)),
                       expression(delta[""<""](k))),bty="n",col=c(1:4),pch=1)
dev.off()


### TEST DATA ###
e <- compileEquation("(L/n)*log(L/n)")
csegs <- consensus(sset, w=c(1,1.01,1,1), e=e, return="breakpoints")

## plot all
png("test_data.png",width=7,height=3.5/2,res=300,units="in")
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



#### TEST POTENTIALS ####
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


### Benchmarking, R vs Rcpp implementation
### Currently not working
# Compare breakpoints of Rcpp and base R implementation
test_that("Breakpoints of consensus in Rcpp implementation are equal to base R implementation",{
    my_e <- "L*L/2"
    e <- compileEquation(my_e)

    expect_equal(cons_c$breakpoints, cons_r$breakpoints)
})


# Compare backtracing pointers between implementations
test_that("Backtracing pointer of consensus in Rcpp implementation are equal to base R implementation",{
    expect_equal(cons_c$values$ptr, cons_r$values$ptr)
})


# Compare Deltas between implementation
test_that("Delta_k of consensus in Rcpp implementation are equal to base R implementation",{
    expect_equal(cons_c$values$Dk, cons_r$values$Dk)
})


# Compare Deltas between implementation
test_that("Delta_star of consensus in Rcpp implementation are equal to base R implementation",{
    expect_equal(cons_c$values$dstar, cons_r$values$dstar)
})

