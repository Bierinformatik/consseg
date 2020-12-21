###TEST AREA
library(ConsSeg)
library(IRanges)
#Load Data
data(primseg436_sset)

#Extract starts and end of segments
segs <- extract_ranges(sset)

# test potential functions
m <- length(sset$ids) #number of segmentations
#n <- sset$N #length of segmented sequence
n <- 1000 # for testing
l <- length(sset$segments$ID)#number of segments

w <- function(m){return(1/m)}
e <- function(width){ return(width^2/2)}

dl <- list()
dle <- list()
dlov <- list()
drov <- list()
F <- list()

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
    ptr[k] <- k
    if (k == 1){
        next
    }
    print(paste0("k ",k))
    for(j in 1:(k-1)) { #j+1
        print(paste0("j ",j))
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


##Backtrace
#for (k in 1:n){
#    if(ptr[[k]] != k){
#        print(F[[k]])
#    }
#}
