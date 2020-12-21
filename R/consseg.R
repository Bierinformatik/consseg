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

    segs <- extract_segments(RS) #list of segment start/end coordinates
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

# Generate lookup for interval overlapping segments
#### HIER WEITER ####
#### for(k=1;k<=n;k++) {
### dsm[k] = dsm[k-1]
### dsq[k] = dsm[k-1]
### dcd[k] = 0
### dcu[k] = 0
### for (m=0, m<M, m++) {
###   if (Bup[k,m]== k) dsm[k] += w(m)*aeh(Bup[k,m]-Bdw[k,m]+1);
###   if (Bdw[k,m]== k) dsq[k] += w(m)*aeh(Bup[k,m]-Bdw[k,m]+1);
###   dcd[k] += w(m)*aeh(k-vor[k,m]+1);
###   dcu[k] += w(m)*aeh(Bup[k,m]-k+1);
### }
### dstar = 0
### F[k]   = infty
### for(j=0;j<k;j++) {
###   /* interval = [j+1,k] */
###     for (m=0, m<M, m++) {
###       if ( (Blw[k]<j+1) && (Bup[k]>k) ) {
###         dtmp = aeh(Bup[k]-Blw[k]+1) + aeh(k-j)
###         - aeh(Blw[k]-j) - aeh(Bup[k]-k+1)
###         dstar += w(m)*dtmp
###       }
###     }
###   Dtmp = dsm[l] - dsq[j+1] + dcd[k] + dcu[j+1] + dstar
###   D = aeh(k-j) - 2*Dtmp
###   /* we don't want to store D, so we keep the pointer */
###     if( F[j] + D < F[k] ) {
###         F[k] = F[j]+D
###         prt[k] = j
###     }
###   }
### }
####
####
#    lookup <- list()
#    lookup_segments(j, k, segs)
#
#    F <- rep(NA, m) ## recursion vector, F[k] = min(scoref(j+1,k) + F[j])
#    jmin <- F       ## backtracing vector: store position j in recursion
#
#
#    for ( j in 1:n) {
#        for ( k in j+1:n ) {
#        }
#        ## F[k] = min(scoref(j+1,k) + F[j])
#        ## j1min ## store position of min
#        score <= scoref(e, j1, k, dl, dle, dlov, drov, ds)
#    }
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
extract_ranges <- function(S){

    #we need to transform this into a list of sequences that contains a list of starts,ends per segment of that sequence
    SO = subset(S$segments, select = c("ID","type","CL","start","end"))
    ranges <- IRanges(start=SO$start, end=SO$end)

    return(ranges)
}

## inverse_lookup
## all segments that overlap my j,k interval
#' @param j interval start
#' @param k interval end
#' @param segs starts/ends vector of segments returning segment width
#' @export
lookup_overlap <- function(j, k, segs) {

    overlap_start <- segs$starts[ i>= which(as.integer(names(segs$starts)) <= k)]
    potential = 0
    for (points in breakpoints){
        for (point in points){
            potential = potential + e(as.integer(point))
        }
    }

    return(potential)

}


## \delta_{<}(i)
## all segments that start and end left of i (used for right border, k)
#' @param i current position
#' @param w weight function
#' @param e potential function
#' @param segs starts/ends vector of segments returning segment width
#' @export
calc_dl <- function(i, w, e, segs) {

    breakpoints <- segs$ends[which(as.integer(names(segs$ends)) < i)]
    potential = 0
    m = lenght(breakpoints)

    for (points in breakpoints){
        m = m + lenght(points) #Sum up the number of breakpoint to make sum weight eq 1
    }
    for (points in breakpoints){
        for (point in points){
            potential = potential + w(m) * e(as.integer(point))
        }
    }

    return(potential)

}

## \delta_{\le}(i)
## all segments that start left of i (used for left border, j+1)
#' @param i current position
#' @param w weight function
#' @param e potential function
#' @param segs starts/ends vector of segments returning segment width

#' @export
calc_dle <- function(i, w, e, segs) {

    breakpoints <- segs$starts[which(as.integer(names(segs$starts)) < i)]
    potential = 0
    m = lenght(breakpoints)

    for (points in breakpoints){
        m = m + lenght(points) #Sum up the number of breakpoint to make sum weight eq 1
    }
    for (points in breakpoints){
        for (point in points){
            potential = potential + w(m) * e(as.integer(point))
        }
    }

    return(potential)

}

## \delta^{\cap}_{<}(i)
## all segments that span i, count left of i (used for right border, k)
#' @param i current position
#' @param w weight function
#' @param e potential function
#' @param segs starts/ends vector of segments returning segment width
#' @export
calc_dlov <- function(i, w, e, segs) {

    breakpoints <- segs$starts[which(as.integer(names(segs$starts)) <= i)]
    #segs$ends[which(as.integer(names(segs$ends)) >= i)]
    potential = 0
    m = 0

    for (seg in segs$starts){
        for (start in seg){
            end = start+-i
            if (diff >0){
               m = lenght(points) #Sum up the number of segments that span to make sum weight eq 1
            }
        }
    }
    for (points in breakpoints){
        m = lenght(points)
        for (point in points){
            diff = points+point-i
            if (diff > 0){
                potential = potential + w(m) * e(as.integer(diff)) #i is the start and we want the cyan part that spans into the j1,k interval
            }
        }
    }

    return(potential)

}

## \delta^{\cap}_{>}(i)
## all segments that span i, count right of i (used left border, j+1)
#' @param i current position
#' @param w weight function
#' @param e potential function
#' @param segs starts/ends vector of segments returning segment width
#' @export
calc_drov <- function(i, w, e, segs) {

    breakpoints <- segs$starts[which(as.integer(names(segs$starts)) <= i)]
    potential = 0
    m = 0

    for (points in breakpoints){
        for (point in points){
            if (points + point > i){
                m = m + lenght(points) #Sum up the number of breakpoints to make sum weight eq 1
            }
        }
    }
    for (points in breakpoints){
        m = lenght(points)
        for (point in points){
            diff = points-(point+i)
            if (diff > 0){
               potential = potential + w(m) * e(as.integer(point-(points-i-1)+1)) # i is the end and we want the magenta part that spans into j1,k
            }
        }
    }

    return(potential)

}

## \delta^*(i',i'')
## all segments that span j+1/k (current left and right border)
#' @param j1 current span start
#' @param k current span end
#' @param w weight function
#' @param e potential function
#' @param segs starts/ends vector of segments returning segment width
#' @export
calc_ds <- function(j1, k, e, segs) {

}


##  \Delta([j+1,k]) \ref{eq:Delta}
## score function
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


##### TRASH #####
##### ##### #####

                                        # CREATE BREAKPOINT REVERSE LOOKUP
#####FUTURE Sparsify from S to SM with j = 6 , instead of SM with j = Nr. of segments
## Along sequence n..m all segments
## starting and ending at interval n..i-1
## all segments starting and ending at i and
## all segements starting and ending in interval i+1..m
#####For now:
#' Matrix from segements S from segmenTier to SM matrix of segmentations,
#' where breakpoints are indicated by a weight \eqn{w\in 0,{0,\dots,1}}
#' where weight is 1/total_nr_of_breakpoints or to be determined
#' @param S list of segmentations (breakpoints)
#' @export
lookup_segments <- function(S){

    segs <- extract_segments(S)
    m <- length(segs)/2
    n <- S$N #Total length of input sequences

    SM <- matrix(0, nrow = m, ncol = n)
                                        # Statt matrix die start/end listen
                                        # Zu jeder Position i den intervall wo bin ich fuer jede segmentierung
                                        # Auf der einen Seite das interval of interest j+1,k
                                        # auf der anderen Seite die Breakpoints
                                        # D(k,l) rekursiv indem
                                        #
                                        # Array der delta erst berechnen als Summe der breakpoints
    for (j in 1:n){
        for (i in 1:m){
            if ( any(sapply(segs[paste0("start.",i)],function(x) x==j)) | any(sapply(segs[paste0("end.",i)],function(x) x==j)) ){
                SM[i:i,j:j] <- w(i,j,m)  # w marks existence of boundary and adds weight of that boundary
                                        # Hier net die matrix sondern die deltas anfuellen
                                        # delta_i durchschnitt ding is immer e(von der intervallcoverage)
            }
        }
    }

    ret <- list(SM, segs, m, n)
    names(ret) <- c("SM", "segs", "m", "n")

    return(ret)
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

