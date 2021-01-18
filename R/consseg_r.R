#' consseg: consensus segmentation from multiple input segmentations
#'@author Halima Saker, Rainer Machne \email{raim@tbi.univie.ac.at}, Jörg Fallmann \email{fall@bioinf.uni-leipzig.de}, Ahmad M. Shahin, Peter F. Stadler \email{studla@bioinf.uni-leipzig.de}
#'@docType package
#'@name consseg
#'@description Calculates consensus segmentation from cluster based segmentation
#'@section Dependencies: The package strictly depends on
#' \code{Rcpp} and \code{RcppXPtrUtils}.
#' All other dependencies are usually present in a
#' basic installation (\code{stats}, \code{graphics}, \code{grDevices})).
#' @references
#' Saker, Machne, Fallmann, Shahin & Stadler (2021) <>,
#' Machne, Murray & Stadler (2017) \doi{10.1038/s41598-017-12401-8},
#' @useDynLib consseg
NULL # this just ends the global package documentation

### DYNAMIC PROGRAMMING BASED CONSENSUS SEGMENTATION OF A CLUSTERING
### implemented by Rainer Machne & Jörg Fallmann

### FUNCTIONS

### MESSAGE UTILS
## nicer time-stamp
time <- function() format(Sys.time(), "%Y%m%d %H:%M:%S")
## messages
msg <- function(x) cat(x, file=stdout()) # until piping is implemented
## stored "warnings" (actually detailed messages)
warn <- function(w, warnings,verb=FALSE) {
    if (verb) cat(w)
    c(warnings,w)
}

#' pre-compile an potential function equation
#' @param e equation string using the interval length \code{L} and the total
#' sequence length \code{n} to calculate potentials during the \code{consseg}
#' recursion, eg. "(L/n)*log(L/n)" to use the negentropy.
# TODO: example using evaluateEquation to evaluate over L and n
#' @export
compileEquation <- function(e="L*L/2") {
    sign <- "long double my_aeh(double L, double n) { return("
    e <- paste(sign, e, ");}")
    e <- RcppXPtrUtils::cppXPtr(e)
    tmp <- try(RcppXPtrUtils::checkXPtr(e, type="long double",
                                        args=c("double","double")))
    if ( "try-error"%in%class(tmp) )
        stop("potential function can not be compile.")
    e
}

#' Calculate consensus segments from a list of segmentation breakpoints
#' @param b list of breakpoints of different segmentations or a
#' \code{data.frame} of segments with "start", "end" and "type columns",
#' or a \code{segmenTier} results object of class "segments".
#' @param n total sequence length (\code{max(b)} if not provided).
#' @param w weights vector, must sum up to 1 or will be normalized.
#' @param e potential function either a \code{XPtr} pointer to a
#' function pre-compiled with \code{\link{compileEquation}} or a string of
#' a (C++-compatible) mathematical equation using \code{L} for the
#' interval length and \code{n} for the total sequence length to
#' calculate potentials during the \code{consseg} recursion,
#' eg. "(L/n)*log(L/n)" to use the negentropy.
#' @param return return class, simple "breakpoints" or "segments", where
#' breakpoints are considered the start, and ends are cut one before.
#' @param store for debugging: store and return all internal vectors.
#' @param test for debugging.
#' @param verb verbosity level, 0 is silent.
## Rcpp bug requires importing it completely
#' @import Rcpp
#' @importFrom RcppXPtrUtils cppXPtr checkXPtr
#'@export
consensus <- function(b, n, w, e,
                      return="breakpoints", store=FALSE, test=FALSE, verb=1) {

    ## get pointer to potential function
    if ( missing(e) )
        e <- e_ptr() # default, aeh function in cpp file
    else if ( class(e)=="XPtr" ) { # test pre-compiled function
        tmp <- try(RcppXPtrUtils::checkXPtr(e, type="long double",
                                            args=c("double","double")))
        if ( "try-error"%in%class(tmp) )
            stop("wrong potential function signature, should be: ",
                 "`long double my_aeh(int L)`.")
    } else if ( class(e)=="character" ) { # compile from string
        if ( verb>0 )
            cat(paste("Compiling user supplied potential function.\n"))
        e <- compileEquation(e)
    }


    ## get class of breakpoints
    if ( "segments"%in%class(b) ) {
        n <- b$N
        blst <- split(b$segments, f=b$segments$type)
        b <- lapply(blst, function(x) c(x$start,x$end))
    } else if ( "data.frame"%in%class(b) ) { # start, end and type columns!
        blst <- split(b, f=b$type)
        b <- lapply(blst, function(x) c(x$start,x$end))
    }
    ## TODO:add iranges to non-executed example code
    ##} else if ( "IRanges"%in%class(b) ) {
    ##    b <- unique(c(IRanges::start(b), IRanges::end(b)))

    if ( sum(unlist(lapply(b, function(x) any(x<1))))>0 )
        stop("Breakpoints must be >0.")

    if ( missing(n) ) {
        n <- max(unlist(b))
        warning("total sequence length `n` is missing, ",
                "using the maximal breakpoint at ", n, ".")
    }

    ## prepare breakpoints such that they:
    ## * are sorted and unique,
    ## * start with 1 and end with n,
    ## and adding n+1 to for convenience in look-up table.
    b <- lapply(b, function(x) sort(unique(c(1,x,n,n+1))))
    M <- length(b)

    ## generate or normalize weight vector
    if ( missing(w) ) w <- rep(1/M, M)
    else if ( length(w)!=M )
        stop("Weight vector must be of the same length as number of",
             " input segmentations.")
    else if ( any(w<0) )
        stop("Weights can not be negative.")
    else if ( sum(w)!=1 ) {
        warning("Weight vector does not sum up to 1, normalizing.")
        w <- w/sum(w)
    }

    ## call recursion in C
    if ( verb>0 )
        cat(paste("Running Recursion.\n"))
    cons <- consensus_c(b, n=n, w=w, e=e, store=store)

    if ( return=="breakpoints" ) res <- cons$breakpoints
    else if ( return=="segments" ) res<-bp2seg(cons$breakpoints,trim.ends=TRUE)

    ## debug option: return all values
    if ( store ) return(cons)
    return(res)
}




#' Calculate consensus segments from a list of segmentation breakpoints
#'
#' This R implementation is maintained for debugging and testing
#' Rcpp code against expected results. If \code{test==TRUE} the
#' Delta of the recursion is additionally calculated directly
#' (very inefficiently) to compare results with an explicitly correct
#' calculation. This can also be useful when testing
#' potential functions that can yield large numbers and thereby cause
#' rounding errors.
#' @param b list of breakpoints of different segmentations.
#' @param n total sequence length (\code{max(b)} if not provided).
#' @param w weights vector, must sum up to 1 or will be normalized.
#' @param e potential function, taking two arguments: the length \code{L}
#' of the evaluated interval and the total sequence length \code{n}.
#' @param store for debugging: store and return all internal vectors.
#' @param test for debugging: calculate Delta directly without using
#' delta helpers. This is very slow but straightforward to implement,
#' the resulting Delta is tested against the Delta of the fast
#' implementation, and results can be used to test other implementations
#' (eg. in Rcpp).
#' @param rel.tol relative error tolerance to report and count differences
#' in Delta during test mode (\code{test==TRUE}).
#'@export
consensus_r <- function(b, n, w, e=function(L,n) L^2/2,
                        store=FALSE, test=FALSE, rel.tol=1e-8) {

    if ( missing(n) ) {
        n <- max(unlist(b))
        warning("total sequence length missing, ",
                "taking the maximal breakpoint at ", n, ".")
    }

    M <- length(b) # number of input segmentations


    ## prepare breakpoints such that they:
    ## * are sorted and unique,
    ## * start with 1 and end with n,
    ## and adding n+1 for convenience in look-up table.
    b <- lapply(b, function(x) sort(unique(c(1,x,n,n+1))))

    ## generate or normalize weight vector
    if ( missing(w) ) w <- rep(1/M, M)
    else if ( sum(w)!=1 ) {
        warning("Weight vector does not sum up to 1, normalizing.")
        w <- w/sum(w)
    }

    ##  FILL UP INTERVAL BORDER LOOKUP TABLES
    Blw <- Bup <- matrix(NA, nrow=n, ncol=M)
    for ( q in 1:M ) {
        i <- 1
        while ( b[[q]][i] <= n ) {
            lw = b[[q]][i]
            up = b[[q]][i+1]-1
            if ( up>n ) up = n
            for ( k in lw:up ) {
                Blw[k,q] = lw
                Bup[k,q] = up
            }
            i <- i+1
        }
    }

    ##  RECURSION

    ## initialize recursion vectors
    F <- rep(0, n)
    dsm <- dsq <-rep(0, n)
    dcd <- dcu <-rep(0, n)
    ptr <- rep(0, n)

    ## store \Delta and \delta^* for debugging
    if ( store )
        Dk <- Ds <- rep(NA,n)

    ## the recursion
    for ( k in 1:n ) {

        if ( k>1 ) {
            dsm[k] = dsm[k-1]
            dsq[k] = dsq[k-1]
        }
        dcd[k] = 0
        dcu[k] = 0

        for ( m in 1:M ) {
            if ( Bup[k,m] == k ) # \delta_<(k), start and end left of k
                dsm[k] = dsm[k] + w[m]*e(Bup[k,m]-Blw[k,m]+1, n)
            if ( Blw[k,m] == k ) # \delta_le(j), start left of j, to subtract
                dsq[k] = dsq[k] + w[m]*e(Bup[k,m]-Blw[k,m]+1, n)
            if ( Bup[k,m] > k ) # \delta^\cap_<(k), left end to k
                dcd[k] = dcd[k] + w[m]*e(k-Blw[k,m]+1, n)
            if ( Blw[k,m] < k ) # \delta^\cap_>(j+1), j+1 to right end
                dcu[k] = dcu[k] + w[m]*e(Bup[k,m]-k+1, n)
        }

        ## /* scan interval = [j+1,k] for minimum j */
        for ( j in 0:(k-1) ) {
            ## \delta^*: correct for segments that span [j+1,k]
            dstar = 0
            for ( m in (1:M) )
                if ( ( Blw[k,m] < j+1 ) & ( Bup[k,m] > k ) ) {
                    dtmp = e(Bup[k,m]-Blw[k,m]+1, n) + e(k-j, n) -
                        e(k-Blw[k,m]+1, n) - e(Bup[k,m]-j, n)
                    dstar = dstar + w[m]*dtmp
                }

            Dtmp = dsm[k] + dcd[k] + dcu[j+1] + dstar
            if ( j>0 ) Dtmp = Dtmp - dsq[j]

            D = e(k-j, n) - 2*Dtmp

            ## straightforward slow implementation for debugging
            ## this should deliver the correct D
            if ( test ) {
                Dslow = e(k-j)    #/* e(A) */
                for( m  in 1:M ) {  #/* loop ueber die input segmenierungen */
                    summe = 0
                    lw = j+1
                    up = Bup[lw,m]
                    while (up < k) {
                        summe = summe + e(up-lw+1,n)
                        lw = up+1
                        up = Bup[lw,m]
                    }
                    summe = summe + e(k-lw+1,n)
                    Dslow = Dslow - 2*w[m]*summe
                }
                if ( abs(1-D/Dslow) > rel.tol ) # test for relative error
                    cat(paste("DIFFERENCE Delta_tested/Delta_correct",
                              k, j, signif(abs(1-D/Dslow),4), "\n"))
                D <- Dslow
            }


            ## find F[k] = min Delta(j+1,k) + F(j)
            ## and store the j that delivered it
            if ( j==0 ) {
                F[k] = D
            } else if ( F[j]+D < F[k] ) {
                F[k] = F[j]+D
                ptr[k] = j
                if ( store ) { # store for debugging or plots
                    Dk[k] <- D
                    Ds[k] <- dstar
                }
            }
        }
    }

    ## BACKTRACE
    bp <- backtrace_r(ptr)
    ## remove helper n+1, which should always be the last
    bp <- bp[1:(length(bp)-1)]

    ## results
    results <- list(breakpoints=bp)
    if ( store )
        results <- append(results,
                          list(values=list(ptr=ptr, F=F,
                                           dsm=dsm, dsq=dsq,
                                           dcu=dcu, dcd=dcd,
                                           dstar=Ds,
                                           Bup=Bup, Blw=Blw,
                                           Dk=Dk)))
    return(results)

}


#' backtrace function for consseg
#' @param ptr pointer of segment ends.
#' @return breakpoints
#' @importFrom utils tail
#' @export
backtrace_r <- function(ptr) {

    ends <- numeric()
    ends <- c(ends, length(ptr))
    k <- tail(ptr, n=1)
    while(!is.na(k) & k>0){
        ends <- c(ends, k)
        k <- ptr[k]
    }

    ends <- rev(ends)
    ends <- ends +1 # j -> j+1
    ends <- unique(c(1,ends)) # ensure that first position is a breapkpoint

    return(ends)
}



#### SOME UTILS

## TODO: do this for a list of breakpoints, as the input for consensus
#' convert a vector of breakpoints incl. 1 and n into a list of
#' adjacent segments, starting at the breakpoints and
#' ending 1 before the next start/breakpoint
#' @param bp vector of breakpoints.
## @param start optional value to add to starts, eg. +1, if breakpoints
## are segment ends.
## @param end optional value to add to ends, eg. -1.
#' @param trim.ends assume that breakpoints are starts of segments,
#' and the ends are the \code{next breakpoint -1}, EXCEPT for
#' the last breakpoint, which corresponds to the total sequence length
#' and is a true end.
#' @return a data.frame with start and end positions of each segment
#' @export
bp2seg <- function(bp, trim.ends=FALSE) {
    df <- data.frame(start=bp[2:length(bp)-1],
                     end=bp[2:length(bp)], type="consensus")
    ## add end (-1), except for LAST
    if ( trim.ends )
        df$end[-length(df$end)] <- df$end[-length(df$end)] -1
    df
}

#' random segmentations
#' @param m number of segmentations.
#' @param n total length of sequence.
#' @param lambda mean of poisson distribution to select number of segments.
#' @importFrom stats rpois rnorm
#' @export
random_breakpoints <- function(m=10,n=50,lambda=5) {

    ## poisson distributed number of segments
    sgnums <- rpois(n=m, lambda=lambda)
    ## normally distributed segment lengths
    blst <- lapply(sgnums, function(x) cumsum(abs(rnorm(x))))
    names(blst) <- paste0("S",1:m)
    ## scale to total length n, round and add 1
    blst <- lapply(blst, function(x) unique(c(1,1+round((n-1)*x/max(x)))))
    blst
}

## TODO: allow breakpoint list (use bp2seg) and table
#' Simple plot function for a list of segmentation
#' breakpoints or segment tables.
#' @param blst a vector of breakpoints (output from \code{\link{consensus}}),
#' a named list of breakpoint vectors, a named list of \code{data.frame}s
#' with start and end columns, or a segmenTier results object
#' (\code{class(blst)=="segments"}).
#' @param n total sequence length (\code{max(unlist(blst))} if not provided).
#' @param add add to existing plot.
#' @param length length of the edges of the arrow head (in inches).
#' @param angle angle from the shaft of the arrow to the edge of the arrow
#' head.
#' @param code integer code, determining _kind_ of arrows to be drawn.
#' @param col arrow colors, a single color or a vector for each segmentation.
#' @param lwd arrow line width.
#' @param axis1 draw x-axis.
#' @param axis2 draw y-axis.
#' @param xlab x-axis label.
#' @param ... further arguments to \code{\link{arrows}}.
#' @importFrom graphics arrows axis mtext par plot
#' @export
plot_breaklist <- function(blst, n, add=FALSE,
                           length=.05, angle=45, code=3, col=1, lwd=2,
                           xlab="position",
                           axis1=TRUE, axis2=TRUE, ...) {

    ## convert list of breakpoints to table
    if ( class(blst)=="numeric" ) { # breakpoint vector
        blst <- list(segments=bp2seg(blst))
    } else if ( "segments"%in%class(blst) )  { # segmenTier class segments
        n <- blst$N
        blst <- split(blst$segments, f=blst$segments$type)
    } else { # list of breakpoint vectors
        if ( class(blst[[1]])!="data.frame" )
            blst <- lapply(blst, bp2seg)
    }
    if ( class(blst)!="list" )
        stop("blst should be a vector of breakpoints, ",
             "a list of breakpoint vectors, ",
             "a segmenTier results object (class=='segments'), ",
             "or a list of segment tables with start, end and type columns.")

    M <- length(blst)
    if ( missing(n) )
        n <- max(unlist(lapply(blst, function(x) max(c(x$start,x$end)))))

    if ( !add ) {
        plot(1:n,col=NA,ylim=c(.5,M+.5),ylab=NA,xlab=NA,
             axes=FALSE)
        if ( axis1 ) {
            axis(1)
            mtext(xlab,1,par("mgp")[1])
        }
        if ( axis2 )
            axis(2, at=1:M, labels=rev(names(blst)), las=2)
    }
    if ( length(col)==1 ) col <- rep(col,length(M))
    for ( i in seq_len(M) ) {
        arrows(x0=blst[[i]]$start, x1=blst[[i]]$end, y0=M-i+1,
               angle=angle,length=length, code=code, col=col[i], lwd=lwd, ...)
    }
}

