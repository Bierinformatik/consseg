###TEST AREA
library(ConsSeg)
#Load Data
data(primseg436_sset)

#Extract starts and end of segments
segs <- extract_ranges(sset)

# test potential functions
n <- length(sset$ids) #number of segmentations
l <- sset$N #length of segmented sequence
c <- length(sset$segments$ID)#number of segments

w <- function(m){return(1/m)}
e <- function(width){ return(width^2/2)}

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

head(dl)
head(dle)
head(dlov)
head(drov)

