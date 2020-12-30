


#' Calculate consensus segments from a list of segmentation breakpoints
#' @param b list of breakpoints of different segmentations
#' @param n total sequence length (\code{max(b)} if not provided)
#' @param w weight function, taking one argument: the index \code{m} of
#' the respective segmentation in the breakpoint list \code{b}
#' @param e potential function, taking one argument: the length \code{L}
#' of the evaluated interval
#' @param store.matrix for debugging: store and return all internal vectors
#' @param test.slow for debugging: compare the incrementally calculated
#' Delta with the very slow direct calculation 
#'@export
consensus_r <- function(b, n, w=function(m) return(1/m), aeh=function(L) L^5/5,
                      store.matrix=FALSE, test.slow=FALSE) {

    if ( missing(n) ) {
        n <- max(unlist(b))
        warning("total sequence length missing, ",
                "taking the maximal breakpoint at ", n)
    }

    M <- length(b) # number of input segmentations

    ## prepare breakpoints such that they:
    ## * are sorted and unique,
    ## * start with 1 and end with n, 
    ## and adding n+1 to for convenience in look-up table.
    b <- lapply(b, function(x) sort(unique(c(1,x,n,n+1))))
    
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
    if ( store.matrix )
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
                dsm[k] = dsm[k] + w(m)*aeh(Bup[k,m]-Blw[k,m]+1)
            if ( Blw[k,m] == k ) # \delta_le(j), start left of j, to subtract
                dsq[k] = dsq[k] + w(m)*aeh(Bup[k,m]-Blw[k,m]+1)
            if ( Bup[k,m] > k ) # \delta^\cap_<(k), left end to k
                dcd[k] = dcd[k] + w(m)*aeh(k-Blw[k,m]+1)
            if ( Blw[k,m] < k ) # \delta^\cap_>(j+1), j+1 to right end
                dcu[k] = dcu[k] + w(m)*aeh(Bup[k,m]-k+1) 
        }
        
        ## /* scan interval = [j+1,k] for minimum j */
        for ( j in 0:(k-1) ) { 
            ## \delta^*: correct for segments that span [j+1,k]
            dstar = 0
            for ( m in (1:M) ) 
                if ( ( Blw[k,m] < j+1 ) & ( Bup[k,m] > k ) ) {
                    dtmp = aeh(Bup[k,m]-Blw[k,m]+1) + aeh(k-j) -
                        aeh(k-Blw[k,m]+1) - aeh(Bup[k,m]-j)
                    dstar = dstar + w(m)*dtmp
                }
        
            Dtmp = dsm[k] + dcd[k] + dcu[j+1] + dstar
            if ( j>0 ) Dtmp = Dtmp - dsq[j] 
        
            D = aeh(k-j) - 2*Dtmp

            ## straightforward slow implementation for debugging
            ## this should deliver the correct D
            if ( test.slow ) {
                Dslow = aeh(k-j)    #/* e(A) */
                for( m  in 1:M ) {  #/* loop ueber die input segmenierungen */
                    summe = 0
                    lw = j+1 
                    up = Bup[lw,m]
                    while (up < k) { 
                        summe = summe + aeh(up-lw+1)
                        lw = up+1 
                        up = Bup[lw,m]
                    } 
                    summe = summe + aeh(k-lw+1)
                    Dslow = Dslow - 2*w(m)*summe 
                }
                if ( abs(D-Dslow)>0.000000000001 ) # TODO: relative error
                    cat(paste("DIFFERENCE", k, j, D, Dslow, "\n"))
                D <- Dslow
            }
            
            
            ## find F[k] = min Delta(j+1,k) + F(j)
            ## and store the j that delivered it
            if ( j==0 ) {
                F[k] = D
            } else if ( F[j]+D < F[k] ) {
                F[k] = F[j]+D
                ptr[k] = j
                if ( store.matrix ) { # store for debugging or plots
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
    if ( store.matrix )
        results <- append(results,
                          list(F=F, Dk=Dk, dstar=Ds,
                               dsm=dsm, dsq=dsq, dcd=dcd, dcu=dcu))
    return(results)
                                        
}


#' @param ptr pointer of segment ends
#' @return breakpoints
#' @export
backtrace_r <- function(ptr) {
    
    ends <- numeric()
    ends <- c(ends, length(ptr))
    k <- tail(ptr, n=1)
    while(!is.na(k) & k>1){
        ends <- c(ends, k)
        k <- ptr[k]
    }

    ends <- rev(ends)
    ends <- ends +1 # j -> j+1
    ends <- unique(c(1,ends)) # ensure that first position is a breapkpoint

    return(ends)
}

## TODO: do this for a list of breakpoints, as the input for consensus
#' convert a vector of breakpoints incl. 1 and n into a list of
#' adjacent segments, starting at the breakpoints and
#' ending 1 before the next start/breakpoint
#' @param bp vector of breakpoints 
#' @return a data.frame with start and end positions of each segment
bp2seg <- function(bp) 
    data.frame(start=bp[2:length(bp)-1],
               end=bp[2:length(bp)], type="consensus")


#' Simulate IRanges of segments
#' @param l length of sequence
#' @param n average number of segments per segmentation
#' @param s upper bound for number of segmentations per sequence
#' @param r repeat first segmentation n times
#' @param df return dataframe
#' @importFrom IRanges disjoin
#' @return An IRanges object of simulated segments
#' @export
simulate_ranges_r <- function(l, n, s, r, df=FALSE){

    total <- IRanges::IRanges(start=1,end=l)
    ranges <- NULL

    if (r){
        nr_of_segments <- sample(s, size = 1)
        if (nr_of_segments %% 2){ #need start/end pairs so %2 == 0
            nr_of_segments = nr_of_segments + 1
        }
        print(paste0("Simulating ", n,
                     " sequences of length ", l, " with ", nr_of_segments,
                     " segments"))
        segmentations <- sort(sample(l, size = nr_of_segments, replace = FALSE))
        starts <- segmentations[c(TRUE, FALSE)]
        ends <- segmentations[c(FALSE, TRUE)]
        subr <- c(IRanges::IRanges(start=starts, end=ends), total)
        ranges <- IRanges::disjoin(subr)
        ranges <- rep(ranges,each=n)
    }
    else{
        print(paste0("Simulating ", n, " sequences of length ", l))
        for (i in 1:n){
            nr_of_segments <- sample(s, size = 1)
            if (nr_of_segments %% 2){ #need start/end pairs so %2 == 0
                nr_of_segments = nr_of_segments + 1
            }
            print(paste0("Sequence ", i, " with ", nr_of_segments, " segments"))
            segmentations <- sort(sample(l, size = nr_of_segments,
                                         replace = FALSE))
            starts <- segmentations[c(TRUE, FALSE)]
            ends <- segmentations[c(FALSE, TRUE)]
            subr <- c(IRanges::IRanges(start=starts, end=ends), total)
            ranges <- c(IRanges::disjoin(subr), ranges)
        }
    }

    if(df){
        ranges <- as.data.frame(ranges)
        counter <- 0
        ranges$type <- sapply(ranges$start, function(s) if(s == 1){paste0("segmentation",counter<<-counter+1)} else paste0("segmentation",counter))
    }

    return(ranges)
}


#' simple plot function for a list of segmentation tables
plot.breaklist <- function(blst, n, add=FALSE,
                           length=.1, angle=45, code=3, col=1, lwd=2,
                           axis1=TRUE, axis2=TRUE, ...) {

    M <- length(blst)
    if ( missing(n) )
        n <- max(unlist(lapply(blst, function(x) max(c(x$start,x$end)))))
    
    if ( !add ) {
        plot(1:n,col=NA,ylim=c(1,M),ylab=NA,xlab=NA,
             axes=FALSE)
        if ( axis1 ) 
            axis(1)
        if ( axis2 ) {
            axis(2, at=1:M, labels=rev(names(blst)), las=2)
            mtext("sequence position",1,par("mgp")[1])
        }
    }
    for ( i in seq_len(M) ) {
        arrows(x0=blst[[i]]$start, x1=blst[[i]]$end, y0=M-i+1,
               angle=angle,length=length, code=code, col=col, lwd=lwd, ...)
    }
}
