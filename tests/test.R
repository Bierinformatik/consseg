###TEST AREA
library(ConsSeg)
library(IRanges)
#Load Data
data(primseg436_sset)

#Extract starts and end of segments
#segs <- simulate_ranges(200,10,12,TRUE) #repeat
segs <- simulate_ranges(200,3,5,FALSE)
plotRanges(segs)

# test potential functions
#m <- length(sset$ids) #number of segmentations
#n <- sset$N #length of segmented sequence
#l <- length(sset$segments$ID)#number of segments
n <- 200 # for testing

w <- function(m){return(1/m)}
e <- function(width){ return(width^2/2)}

dl <- numeric()
dle <- numeric()
dlov <- numeric()
drov <- numeric()
F <- numeric()
ptr <- numeric(1)

for (k in 1:n){
    dl[k] <- calc_dl(k, w, e, segs)
    dle[k] <- calc_dle(k, w, e, segs)
    dlov[k] <- calc_dlov(k, w, e, segs)
    drov[k] <- calc_drov(k, w, e, segs)
    F[k] <- .Machine$integer.max #as close to infinite as is comes
    dstar = 0
    if (k == 1){
        next #just need to initialize F[k] and pointer here
    }
    print(paste0("k ",k))
    for(j in 1:(k-1)) { #j+1
        dstar <- calc_ds(j, k, w, e, segs)
        Dtmp = scoref(j, k, dl, dle, dlov, drov, dstar)
        D = e(k-j) - 2*Dtmp
        ##we don't want to store D, so we keep the pointer
        if( (F[j] + D) < F[k] ) { # Da hats was, F[k] is hier fast immer groesser als F[j] + D
            #print(paste0("FOUND NEW MINIMA ",(F[j] + D)," lower than ",F[k]," for k ",k," and j ",j))
            F[k] = F[j] + D
            ptr[k] = j
        }
    }
}

##Backtrace Kette nach von letztem j nach 0
for (k in n:2){
    if(ptr[k]){
        print(paste0("k ",k," : ",ptr[k]))
        #print(ptr[[k]])
        k = ptr[k]
    }
}


####MANUAL


dl <- numeric()
dle <- numeric()
dlov <- numeric()
drov <- numeric()
F <- numeric()
ptr <- numeric(1)


for (k in 1:n){
    dl[k] <- calc_dl(k, w, e, segs)
    dle[k] <- calc_dle(k, w, e, segs)
    dlov[k] <- calc_dlov(k, w, e, segs)
    drov[k] <- calc_drov(k, w, e, segs)
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
