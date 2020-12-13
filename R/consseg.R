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
#' Calculate consensus segements C from segmenTier segments
#' @param RS, raw segmentation object from segmenTier
#' @export
consensus <- function(RS) {

  ##collapse will create extremely sparse matrix of nxm dimensions
  #S <- collapse_segments(RS)
  ##extract will retain a data.frame of segments ordered by start/end
  S <- extract_segments(RS)
  m <- nrow(S)
  start <- S[order(S[,"start"], S[,"end"]), ][1,"start"]
  end <- S[order(S[,"start"], S[,"end"]), ][m,"end"]
  n <- end-start#+1

  F <- rep(NA, m) ## recursion vector, F[k] = min(scoref(j+1,k) + F[j])
  jmin <- F       ## backtracing vector: store position j in recursion
  SM <- matrix(0, nrow = m, ncol = n)

  for ( k in 1:n) {
    for ( q in 1:m ) {
    }
    ## F[k] = min(scoref(j+1,k) + F[j])
    ## j1min ## store position of min
  }
}


# EXTRACT SEGMENTS
#' EXTRACT segements from segmenTier to ordered data.frame SO
#' @param S list of segmentations (breakpoints)
#' @export
extract_segments <- function(S){

  #we need to transform this into a list of sequences that contains a list of starts,ends per segment of that sequence
  SO = subset(S$segments, select = c("ID","type","CL","start","end"))
  #SO = SO[order(SO[,"start"], SO[,"end"]), ] # no longer need this

  startlist <- split(SO$start, SO$type)
  names(startlist) <- paste0("start.",1:length(startlist))
  endlist <- split(SO$end, SO$type)
  names(endlist) <- paste0("end.",1:length(endlist))
  returnlist <- c(startlist, endlist)

  return(returnlist)
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
#' @export
collapse_segments <- function(S, w){

    segs <- extract_segments(S)

    segstart <- segs[order(segs[,"start"], segs[,"end"]), ][1,"start"] #order not needed if segs already ordered
    segend <- segs[order(segs[,"start"], segs[,"end"]), ][nrow(segs),"end"] #order not needed if segs already ordered
    m <- nrow(segs)
    #n <- segend-segstart+1 # or do we need from 1-end of real sequence?
    n <- segend # seems so

    if (!is.function(w) & w == 1){
      w = 1/segnr #For now we only normalize by nr of segments, each segment has same weight
    }

    SM <- matrix(0, nrow = m, ncol = segend)

    for (j in 1:n){
      for (i in 1:m){
        if (i %in% segs["CL"] & (j %in% segs["start"] | j %in% segs["end"])){
          SM[i,j] = w  # w marks existence of boundary and adds weight of that boundary
        }
      }
    }

  return(SM)
}


####From RAIM
## \eqn{\e}
## function that evaluates an individual segment
#' @export
sval <- function() {}

##  \Delta([j+1,k]) \ref{eq:Delta}
## score function,
#' @export
scoref <- function(j1,k)
  sval(j1,k) -2*(dl(k) - dle(j1) + dlov(k) + drov(j1) + ds(j1,k))

## \delta_{<}(i)
## all segments that start and end left of i (used for right border, k)
#' @export
dl <- function(i) {

}

## \delta_{\le}(i)
## all segments that start left of i (used for left border, j+1)
#' @export
dle <- function(i) {

}

## \delta^{\cap}_{<}(i)
## all segments that span i, count left of i (used for right border, k)
#' @export
dlov <- function(i) {

}

## \delta^{\cap}_{>}(i)
## all segments that span i, count right of i (used left border, j+1)
#' @export
drov <- function(i) {

}

## \delta^*(i',i'')
## all segments that span j+1/k (current left and right border)
#' @export
ds <- function(i) {}


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

