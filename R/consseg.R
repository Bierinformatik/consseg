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
### e..potential function, w_q .. weight of segment(boundary), kommt vom user
## We need
## Gleichung 9 und 10 vorher, dann 8, und damit i das DELTA bekomm fuer 8 muss i aus den deltas die Gleichung 12 ausrechnen und dafuer on demand delta* entweder ausrechnen wenn segment das intervall schneidet, sonst is das eh 0
## Eventuell inverser index mit intervallgrenzen der mir sagt gibts fuer mein i-k intervall irgendwas wo ein intervall drueber haengt in O(1) -> das brauch ich fuer delta* in Gl 12
## Inverser index fuer jede segmentierung damit ma obere und untere intervallgrenze hat
## 8&9 muss ma nur einmal befuellen weil die ueber alle segmente zaehlen
## und dazu muss ma nur was aendern wenn ma eine der segmentgrenzen ueberschreitet und net fuer jede position
##


###RECURSION
#' Calculate consensus segements C from segmenTier segments
#' @param RS, raw segmentation object from segmenTier
#' @export
consensus <- function(RS) {

  ##matrixfy will create extremely sparse matrix of nxm dimensions, mostly made up of 0's
  ##and a list of segment start/end coordinates
  S <- matrixfy_segments(RS)

  SM <- S$SM
  segs <- S$segs
  m <- S$m
  n <- S$n
  S <- NULL

  F <- rep(NA, m) ## recursion vector, F[k] = min(scoref(j+1,k) + F[j])
  jmin <- F       ## backtracing vector: store position j in recursion


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

  startlist <- split(SO$start, SO$type)
  names(startlist) <- paste0("start.",1:length(startlist))
  endlist <- split(SO$end, SO$type)
  names(endlist) <- paste0("end.",1:length(endlist))
  returnlist <- c(startlist, endlist)

  return(returnlist)
}



# CREATE SEGMENTS MATRIX
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
#' @param w weightfunction or 1
#' @export
matrixfy_segments <- function(S, w){

  segs <- extract_segments(S)
  m <- length(segs)/2
  n <- S$N #Total length of input sequences

  if (!is.function(w)){ #if not function we replace with dummy, assuming a sensible weight functions needs i, j and m for normalization
      w <- function(i,j,m){ i=0; j=0; return(1/m)} #For now we only normalize by nr of segments, each segment has same weight
      # Der user gibt irgendwas an und wir scalieren des so das die Summer immer 1 is
  }

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


####From RAIM
## \eqn{\e}
## function that evaluates an individual segment
#' @export
sval <- function() {
  # Gleichung 8/12, e gibt laenge des segments zurueck
  # [j +1, k] intervall der einen interessiert
  # Max so viele Eintraege wie inputsegmentierungen
  # e is wieder user definiert, sowas wie laenge Intervall ^ 2 /2 oder entropie
}

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
ds <- function(i) {

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

