###TEST AREA
library(ConsSeg)
library(IRanges)
#Load Data
data(primseg436_sset)

#Extract starts and end of segments
segs <- simulate_ranges(200,10,12,TRUE)

# test potential functions
#m <- length(sset$ids) #number of segmentations
#n <- sset$N #length of segmented sequence
#l <- length(sset$segments$ID)#number of segments
m <- length(segs)/10
n <- 200 # for testing
l <- 10

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
    F[i] = .Machine$integer.max #as close to infinite as is comes
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
        dstar <- calc_ds(j, k, w, e, segs)
        Dtmp = scoref(j, k, dl, dle, dlov, drov, dstar)
        D = e(k-j) - 2*Dtmp
        ##we don't want to store D, so we keep the pointer
        if( (F[[j]] + D) < F[[k]] ) { # Da hats was, F[k] is hier ja noch Inf, also immer groesser als F[j] + D
            print(paste0("FOUND NEW MINIMA ",(F[[j]] + D)," for k ",k," and j ",j))
            F[[k]] = F[[j]] + D
            ptr[k] = j
        }
    }
}

##Backtrace Kette nach von letztem j nach 0
for (k in n:1){
    print(paste0("k ",k," : ",ptr[[k]]))
    if(ptr[[k]] != k){
        print(ptr[[k]])
        #k = ptr[[k]]
    }
}


##TODO
##Simulieren 1000 lang, 30-40 Segmente, zufaellig segmentieren
##Potentialfunktion x^3/x^4 ob da was passiert wos nimma passt
##Peter meint des stimmt wohl immer, aber mal schaun
##
##Raim soll mal schauen nach Beispiel, vielleicht aus dem anderen Beispiel
##aus dem segmentier Papierl
##Eventuell no ein Rcpp draus basteln
##
##Test drei oder vier mal die gleiche segmentierung,
##dann sieht ma flott ob da was net passt
