###TEST AREA
library(ConsSeg)
#Load Data
data(primseg436_sset)

#Extract starts and end of segments
segs <- extract_segments(sset)

# test potential functions
m <- length(sset$ids) #number of segmentations
n <- sset$N #length of segmented sequence
l <- length(sset$segments$ID)#number of segments

w <- function(i,j,m){ i=0; j=0; return(1/m)}
e <- function(width){ return(width^2/2)}

dl <- list()
dle <- list()
dlov <- list()
drov <- list()

for (i in 1:n){
    dl[i] <- calc_dl(i, e, segs)
    dle[i] <- calc_dle(i, e, segs)
    dlov[i] <- calc_dlov(i, e, segs)
    drov[i] <- calc_drov(i, e, segs)
}



