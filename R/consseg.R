#' consseg: consensus segment for segmenTier cluster-based segmentation
#' from a sequential clustering
#'@author Halima Saker, Rainer Machne \email{raim@tbi.univie.ac.at}, Jörg Fallmann \email{fall@bioinf.uni-leipzig.de},
#' Ahmad M. Shahin, Peter F. Stadler \email{studla@bioinf.uni-leipzig.de}
#'@docType package
#'@name consseg
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
    F <- list() #

    for (i in 1:n){
        dl[i] <- calc_dl(i, w, e, segs)
        dle[i] <- calc_dle(i, w, e, segs)
        dlov[i] <- calc_dlov(i, w, e, segs)
        drov[i] <- calc_drov(i, w, e, segs)
        F[i] = Inf
    }

    dstar = 0
    ptr <- list()

    for (k in 1:n){
        ptr[k] = k
        for(j in 1:k) { #i+1
            dstar <- calc_ds(j+1, k, e, segs)
            Dtmp = dl[l] - dle[j+1] + dlov[k] + drov[j+1] + dstar
            D = e(k-j+1) - 2*Dtmp
            ##we don't want to store D, so we keep the pointer
            if( F[j] + D < F[k] ) {
                F[k] = F[j] + D
                ptr[k] = j
            }
        }
    }
}

# EXTRACT SEGMENTS
#' EXTRACT segements from segmenTier to ordered data.frame SO
#' @param S list of segmentations
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
#' @export
extract_ranges <- function(S){

    #we need to transform this into a list of sequences that contains a list of starts,ends per segment of that sequence
    SO = subset(S$segments, select = c("ID","type","CL","start","end"))
    total <- IRanges(start=1,end=S$N)
    ranges <- NULL

    for (t in unique(SO$type)){
        sub <- subset(SO, type==t)
        subr <- c(IRanges(start=sub$start, end=sub$end),total)
        ranges <- c(disjoin(subr), ranges)
    }

    return(ranges)
}

## \delta_{<}(i)
## all segments that start and end left of i (used for right border, k)
#' Calculates \delta_{<}(i)
#' @param i current position
#' @param w weight function
#' @param e potential function
#' @param segs starts/ends vector of segments returning segment width
#' @export
calc_dl <- function(i, w, e, segs) {

    breakpoints <- subset(segs, end < i)
    potential = 0

    for (len in width(breakpoints)){
                potential = potential + w(length(breakpoints)) * e(len)
    }

    return(potential)

}

## \delta_{\le}(i)
## all segments that start left of i (used for left border, j+1)
#' Calculates \delta_{\le}(i)
#' @param i current position
#' @param w weight function
#' @param e potential function
#' @param segs starts/ends vector of segments returning segment width

#' @export
calc_dle <- function(i, w, e, segs) {

    breakpoints <- subset(segs, end < i)
    potential = 0

    for (len in width(breakpoints)){
        potential = potential + w(length(breakpoints)) * e(len)
    }

    return(potential)

}

## \delta^{\cap}_{<}(i)
## all segments that span i, count left of i (used for right border, k)
#' Calculates \delta^{\cap}_{<}(i)
#' @param i current position
#' @param w weight function
#' @param e potential function
#' @param segs starts/ends vector of segments returning segment width
#' @export
calc_dlov <- function(i, w, e, segs) {

    breakpoints <- subset(segs, start <= i & end >= i)
    diff <- restrict(breakpoints, end = i)

    potential = 0

    for (len in width(diff)){
        potential = potential + w(length(diff)) * e(len)
    }

    return(potential)

}

## \delta^{\cap}_{>}(i)
## all segments that span i, count right of i (used left border, j+1)
#' Calculates \delta^{\cap}_{>}(i)
#' @param i current position
#' @param w weight function
#' @param e potential function
#' @param segs starts/ends vector of segments returning segment width
#' @export
calc_drov <- function(i, w, e, segs) {

    breakpoints <- subset(segs, start <= i & end >= i)
    diff <- restrict(breakpoints, start = i)

    potential = 0

    for (len in width(diff)){
        potential = potential + w(length(diff)) * e(len)
    }

    return(potential)

}

## \delta^*(i',i'')
## all segments that span j+1/k (current left and right border)
#' Calculates \delta^*(i',i'')
#' @param j1 current span start
#' @param k current span end
#' @param w weight function
#' @param e potential function
#' @param segs starts/ends vector of segments returning segment width
#' @export
calc_ds <- function(j1, k, e, segs) {

    look <- IRanges(start=j, end=k)  #Check to make sure we got boundaried right, here we have direct overlap
    ov <- segs[subjectHits(findOverlaps(look, segs, type="within"))]
    diffj <- restrict(ov, end = j1)
    diffk <- restrict(ov, start = k)

    potential = 0

    for (len in width(ov)){
        potential = potential + e(len) + e(k-j1+1) - e(diffj) -  e(diffk)
    }

    potential = potential * w(length(ov))

    return(potential)

}


##  \Delta([j+1,k]) \ref{eq:Delta}
## score function
#' Calculates \Delta([j+1,k])
#' @param e potential function
#' @param j1 current span start
#' @param k current span end
#' @param dl dl list
#' @param dle dle list
#' @param dlov dlov list
#' @param drov drov list
#' @param ds delta star
#' @export
scoref <- function(e, j1, k, dl, dle, dlov, drov, ds){
    return (e(j1,k) -2*(dl(k) - dle(j1) + dlov(k) + drov(j1) + ds(j1,k)))
}


### DATA SET DOC

#' Segements from transcriptome time-series from budding yeast.
#'
#' segmenTier processed transcriptome time-series data from a region encompassing
#' four genes and a regulatory upstream non-coding RNA in budding yeast.
#' The data set is described in more detail in the publication
#' Machne, Murray & Stadler (2017) <doi:10.1038/s41598-017-12401-8>.
#'
#' @name primseg436_sset
#' @docType data
NULL

