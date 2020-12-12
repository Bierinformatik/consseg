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
##'@importFrom segmenTier segments
##'@importFrom Rcpp evalCpp
##'@importFrom stats qt sd var BIC AIC
##'@importFrom graphics image axis par plot matplot points lines legend arrows strheight strwidth text mtext abline polygon
##'@importFrom grDevices png dev.off rainbow gray xy.coords
NULL # this just ends the global package documentation

### DYNAMIC PROGRAMMING BASED CONSENSUS SEGMENTATION OF A CLUSTERING
### implemented by Jörg Fallmann, based
### on Rainer Machne's segmenTier implementation


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
### e..potential function, w_q .. weight of segment(boundary),
## We need
## dli = sum(1..m){w_q}*sum(seg_end<i){e()}


###RECURSION
#' @param S, segmentation object from segmenTier
consensus <- function(S) {

  SM = collapse_segments(S)
  F <- rep(NA, n) ## recursion vector, F[k] = min(scoref(j+1,k) + F[j])
  jmin <- F       ## backtracing vector: store position j in recursion

  for ( k in 1:n) {
    for ( q in 1:m ) {
    }
    ## F[k] = min(scoref(j+1,k) + F[j])
    ## j1min ## store position of min
  }
}

# COLLAPSE SEGMENTS
#####FUTURE Sparsify from S to SM with j = 6 , instead of SM with j = Nr. of segments
## Along sequence n..m all segments
## starting and ending at interval n..i-1
## all segments starting and ending at i and
## all segements starting and ending in interval i+1..m
#####For now:
#' Collapse segements S from segmenTier to SM matrix of segmentations,
#' where breakpoints are indicated by a weight \eqn{w\in 0,{0,\dots,1}}
#' where weight is 1/total_nr_of_breakpoints or to be determined
#' @param S list of segmentations (breakpoints)
#' @param w weightfunction or 1
collapse_segments <- function(S, w){

    segs = subset(S$segments, select = c("CL","start","end"))
    segstart = segs[order(segs[,"start"], segs[,"end"]), ][1,"start"]
    segend = segs[order(segs[,"start"], segs[,"end"]), ][nrow(segs),"end"]
    segnr = nrow(segs)
    seglen = segend-segstart+1

    if (!is.function(w) & w == 1){
      w = 1/segnr
    }

    SM <- matrix(0,nrow = segnr, ncol = seglen)

    for (i in segstart:segend){
      for (j in 1:segnr){
        SM[i,j] =
      }
    }

  return(SM)
}


####From RAIM
## \eqn{\e}
## function that evaluates an individual segment
sval <- function() {}

##  \Delta([j+1,k]) \ref{eq:Delta}
## score function,
scoref <- function(j1,k)
  sval(j1,k) -2*(dl(k) - dle(j1) + dlov(k) + drov(j1) + ds(j1,k))

## \delta_{<}(i)
## all segments that start and end left of i (used for right border, k)
dl <- function(i) {}

## \delta_{\le}(i)
## all segments that start left of i (used for left border, j+1)
dle <- function(i) {}

## \delta^{\cap}_{<}(i)
## all segments that span i, count left of i (used for right border, k)
dlov <- function(i) {}

## \delta^{\cap}_{>}(i)
## all segments that span i, count right of i (used left border, j+1)
drov <- function(i) {}

## \delta^*(i',i'')
## all segments that span j+1/k (current left and right border)
ds <- function(i) {}

