#' consseg: consensus segment for segmenTier cluster-based segmentation
#' from a sequential clustering
#'@author Halima Saker, Rainer Machne \email{raim@tbi.univie.ac.at}, Jörg Fallmann \email{fall@bioinf.uni-leipzig.de},
#' Ahmad M. Shahin, Peter F. Stadler \email{studla@bioinf.uni-leipzig.de}
#'@docType package
#'@name consseg
#'@description Calculates consensus segmentation from cluster based segmentation
#'@section Dependencies: The package strictly depends on
#' \code{segmenTier} and thus on \code{Rcpp}.
#' All other dependencies are usually present in a
#' basic installation (\code{stats}, \code{graphics}, \code{grDevices})).
#' @references
#' Saker, Machne, Fallmann, Shahin & Stadler (2021) <>,
#' Machne, Murray & Stadler (2017) <doi:10.1038/s41598-017-12401-8>,
#' Machne & Murray (2012) <doi:10.1371/journal.pone.0037906>, and
#' Lehmann et al. (2013) <doi:10.1186/1471-2105-14-133>
##'@import segmenTier
##'@import IRanges
NULL # this just ends the global package documentation

### DYNAMIC PROGRAMMING BASED CONSENSUS SEGMENTATION OF A CLUSTERING
### implemented by Jörg Fallmann & Rainer Machne

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

###TODO Check what is really needed
### We define
### Scoring Segments d[i,j], E .. potential,
### e..potential function, w_q .. weight of segment(boundary), kommt vom user
## We need
## Gleichung 9 und 10 vorher, dann 8, und damit i das DELTA bekomm fuer 8 muss i aus den deltas die Gleichung 12 ausrechnen und dafuer on demand delta* entweder ausrechnen wenn segment das intervall schneidet, sonst is das eh 0
## Eventuell inverser index mit intervallgrenzen der mir sagt gibts fuer mein i-k intervall irgendwas wo ein intervall drueber haengt in O(1) -> das brauch ich fuer delta* in Gl 12
## Inverser index fuer jede segmentierung damit ma obere und untere intervallgrenze hat
## 9&10 muss ma nur einmal befuellen weil die ueber alle segmente zaehlen
## und dazu muss ma nur was aendern wenn ma eine der segmentgrenzen ueberschreitet und net fuer jede position
##


###RECURSION
#' Calculate consensus segements C from segmenTier segments
#' @param RS, raw segmentation object from segmenTier
#' @param w weight function or 1
#' @param e potential function
#' @importFrom methods is
#' @return The consesus segmentation as IRanges object
#' @export
consensus <- function(RS, w, e) {

    segs <- extract_ranges(RS) #list of segment ranges
    m <- length(RS$ids) #number of segmentations
    n <- RS$N #length of segmented sequence
    l <- length(RS$segments$ID)#number of segments
    S <- NULL

    if (!is.function(w)){ #if not function we replace with dummy, assuming a sensible weight functions needs i, j and m for normalization, we now assume the user inputs some value and we strictly normalize such that sum(w) = 1
        w <- function(m){return(1/m)} #For now we only normalize by nr of segments, each segment has same weight
    }

    if (!is.function(e)){ #if not function we replace with dummy. Assuming a sensible potential function is based on segment length we just say lenght^2/2. Could be entropie function of whatever
        e <- function(width){ return(width^2/2)}
                                        # Gleichung 8/12, e gibt laenge des segments zurueck, e is wieder user definiert, sowas wie laenge Intervall ^ 2 /2 oder entropie
                                        # [j +1, k] intervall der einen interessiert
                                        # Max so viele Eintraege wie inputsegmentierungen
    }

    #Precalculate the potentials for dl, dle, dlov, drov over all i's

    dl <- list()
    dle <- list()
    dlov <- list()
    drov <- list()

    for (i in 1:n){
        dl[i] <- calc_dl(i, w, e, segs)
        dle[i] <- calc_dle(i, w, e, segs)
        dlov[i] <- calc_dlov(i, w, e, segs)
        drov[i] <- calc_drov(i, w, e, segs)
    }

    dstar = 0
    F <- list()
    ptr <- list()

    for (k in 1:n){
        F[k] <- .Machine$integer.max #close to infinite
        ptr[k] <- k
        if (k == 1){
            next
        }
        for(j in 1:(k-1)) { #j+1
            dstar <- calc_ds(j, k, w, e, segs)
            Dtmp = scoref(j, k, dl, dle, dlov, drov, dstar)
            D = e(k-j) - 2*Dtmp
            ##we don't want to store D, so we keep the pointer
            if( (F[[j]] + D) < F[[k]] ) {
                F[[k]] = F[[j]] + D
                ptr[k] = j
            }
        }
    }

    ##Backtrace Kette von letztem k nach 0
    for (k in n:1){
        if(ptr[[k]] != k){
            print(ptr[[k]])
            k = ptr[[k]]
        }
    }

}

# EXTRACT SEGMENTS
#' EXTRACT segements from segmenTier to ordered data.frame SO
#' @param S list of segmentations
#' @return A list of segment starts and ends (breakpoints)
#' @export
extract_segments <- function(S){

    #we need to transform this into a list of sequences that contains a list of starts,ends per segment of that sequence
    SO = subset(S$segments, select = c("ID","type","CL","start","end"))
    SO$width <- SO$end-SO$start+1

    startlist <- split(SO$width, SO$start)
    endlist <- split(SO$width, SO$end)
    returnlist <- list("starts" = startlist, "ends" = endlist)

    return(returnlist)

}

#### NOTE: Create IRanges object for easy intersection of position
#' EXTRACT segements from segmenTier to IRanges
#' @param S list of segmentations
#' @importFrom IRanges disjoin
#' @return An IRanges object of segement blocks
#' @export
extract_ranges <- function(S){

    SO = subset(S$segments, select = c("ID","type","CL","start","end"))
    total <- IRanges(start=1,end=S$N)
    ranges <- NULL

    for (t in unique(SO$type)){
        sub <- subset(SO, SO$type==t)
        subr <- c(IRanges(start=sub$start, end=sub$end),total)
        ranges <- c(disjoin(subr), ranges)
    }

    return(ranges)
}

## \delta_{<}(i)
## all segments that start and end left of i (used for right border, k)
#' Calculates delta(<i)
#' @param i current position
#' @param w weight function
#' @param e potential function
#' @param segs starts/ends vector of segments returning segment width
#' @return delta(<i)
#' @export
calc_dl <- function(i, w, e, segs) {

    bps <- subset(segs, end <= i) #<i or <=i ?
    potential = 0

    if (length(bps) < 1){
        return(0)
    }

    for (len in width(bps)){
                potential = potential + w(length(bps)) * e(len)
    }

    return(potential)

}

## \delta_{\le}(i)
## all segments that start left of i (used for left border, j+1)
#' Calculates delta(<=i)
#' @param i current position
#' @param w weight function
#' @param e potential function
#' @param segs starts/ends vector of segments returning segment width
#' @return delta(<=i)
#' @export
calc_dle <- function(i, w, e, segs) {

    bps <- subset(segs, start < i) #<i or <= 1 ?
    potential = 0

    if (length(bps) < 1){
        return(0)
    }

    for (len in width(bps)){
        potential = potential + w(length(bps)) * e(len)
    }

    return(potential)

}

## \delta^{\cap}_{<}(i)
## all segments that span i, count left of i (used for right border, k)
#' Calculates delta^(<i)
#' @param i current position
#' @param w weight function
#' @param e potential function
#' @param segs starts/ends vector of segments returning segment width
#' @importFrom IRanges restrict
#' @return delta^(<i)
#' @export
calc_dlov <- function(i, w, e, segs) {

    breakpoints <- subset(segs, start < i & end > i) # < > or <= >=?
    diff <- restrict(breakpoints, end = (i-1)) # i-1?

    if (length(diff) < 1){
        return(0)
    }

    potential = 0

    for (len in width(diff)){
        potential = potential + w(length(diff)) * e(len)
    }

    return(potential)

}

## \delta^{\cap}_{>}(i)
## all segments that span i, count right of i (used left border, j+1)
#' Calculates delta^(>i)
#' @param i current position
#' @param w weight function
#' @param e potential function
#' @param segs starts/ends vector of segments returning segment width
#' @importFrom IRanges restrict
#' @return delta^(>i)
#' @export
calc_drov <- function(i, w, e, segs) {

    breakpoints <- subset(segs, start < i & end > i) # < > or <= >= ?
    diff <- restrict(breakpoints, start = (i+1)) # i+1?

    if (length(diff) < 1){
        return(0)
    }

    potential = 0

    for (len in width(diff)){
        potential = potential + w(length(diff)) * e(len)
    }

    return(potential)

}

## \delta^*(i',i'')
## all segments that span j+1/k (current left and right border)
#' Calculates delta*(i',i'')
#' @param j current span start
#' @param k current span end
#' @param w weight function
#' @param e potential function
#' @param segs starts/ends vector of segments returning segment width
#' @import IRanges
#' @import S4Vectors
#' @return delta*(i',i'')
#' @export
calc_ds <- function(j, k, w, e, segs) {

    look <- IRanges(start=(j+1), end=k)  #Check to make sure we got boundaries right, here we have direct overlap
    ov <- segs[subjectHits(findOverlaps(look, segs, type="within"))]

    if (length(ov) < 1){
        return(0)
    }

    diffj <- restrict(ov, end = (j)) # -1?
    diffk <- restrict(ov, start = (k)) # +1?

    potential = 0

    for (len in width(diffj)){
        potential = potential - e(len)
    }

    for (len in width(diffk)){
        potential = potential - e(len)
    }

    for (len in width(ov)){
        potential = potential + e(len) + e(k-j)
    }

    potential <- potential * w(length(ov))

    return(potential)

}


##  \Delta([j+1,k]) \ref{eq:Delta}
## score function
#' Calculates Delta(j+1,k)
#' @param j current span start
#' @param k current span end
#' @param dl dl list
#' @param dle dle list
#' @param dlov dlov list
#' @param drov drov list
#' @param dstar delta star
#' @return Delta(j+1,k)
#' @export
scoref <- function(j, k, dl, dle, dlov, drov, dstar){
    return(dl[[k]] - dle[[j]] + dlov[[k]] + drov[[(j+1)]] + dstar)
}

#' Simulate IRanges of segments
#' @param l length of sequence
#' @param n number of sequences
#' @param s upper bound for number of segmentations per sequence
#' @param r repeat first segmentation n times
#' @importFrom IRanges disjoin
#' @return An IRanges object of simulated segments
#' @export
simulate_ranges <- function(l, n, s, r){

    total <- IRanges(start=1,end=l)
    ranges <- NULL

    if (r){
        nr_of_segments <- sample(s, size = 1)
        if (nr_of_segments %% 2){
            nr_of_segments = nr_of_segments + 1
        }
        print(paste0("Simulating ", n, " sequences of length ", l, " with ", nr_of_segments, " segments"))
        segmentations <- sort(sample(l, size = nr_of_segments, replace = FALSE))
        starts <- segmentations[c(TRUE, FALSE)]
        ends <- segmentations[c(FALSE, TRUE)]
        subr <- c(IRanges(start=starts, end=ends), total)
        ranges <- disjoin(subr)
        ranges <- rep(ranges,each=n)
    }
    else{
        print(paste0("Simulating ", n, " sequences of length ", l))
        for (i in 1:n){
            nr_of_segments <- sample(s, size = 1)
            if (nr_of_segments %% 2){
                nr_of_segments = nr_of_segments + 1
            }
            print(paste0("Sequence ", i, " with ", nr_of_segments, " segments"))
            segmentations <- sort(sample(l, size = nr_of_segments, replace = FALSE))
            starts <- segmentations[c(TRUE, FALSE)]
            ends <- segmentations[c(FALSE, TRUE)]
            subr <- c(IRanges(start=starts, end=ends), total)
            ranges <- c(disjoin(subr), ranges)
        }
    }

    return(ranges)
}


#' Plot Ranges of segments
#' Adopted from [IRangesOverview](https://www.bioconductor.org/packages/devel/bioc/vignettes/IRanges/inst/doc/IRangesOverview.pdf)
#' @param x ranges
#' @param xlim x-axis limit
#' @param main plot main title
#' @param col color for plot
#' @param sep separator for consecutive segment blocks
#' @return A plot of segments
#' @import graphics
#' @importFrom methods is
#' @export
plotRanges <- function(x, xlim=x, main=deparse(substitute(x)), col="black", sep=0.5){

    height <- 1
    if (is(xlim, "IntegerRanges")){
        xlim <- c(min(start(xlim)), max(end(xlim)))
    }
    bins <- disjointBins(IRanges(start(x), end(x) + 1))
    plot.new()
    plot.window(xlim, c(0, max(bins)*(height + sep)))
    ybottom <- bins*(sep + height) - height
    rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col=col)
    title(main)
    axis(1)

    return(plot)
}
