###TEST AREA
library(ConsSeg)
library(IRanges)
#Load Data
data(primseg436_sset)

#Extract starts and end of segments
#segs <- simulate_ranges(200,10,12,TRUE) #repeat
#plotRanges(segs,rep=TRUE)
segs <- simulate_ranges(200,10,12,FALSE)
plot_Ranges(segs)
segs

# test potential functions
#m <- length(sset$ids) #number of segmentations
#n <- sset$N #length of segmented sequence
#l <- length(sset$segments$ID)#number of segments
n <- max(end(segs)) # for testing
M <- length(subset(start(segs), start(segs) == 1))
w <- function(m){return(1/m)}
e <- function(width){ return(width^2/2)}
w <- w(M)

dl <- numeric()
dle <- numeric()
dlov <- numeric()
drov <- numeric()
F <- numeric()
ptr <- numeric(1)

bla <- character()
blubb <- character()
blibb <- character()

for (k in 1:n){
    dl[k] <- calc_dl(k, w, e, segs)
    dle[k] <- calc_dle(k, w, e, segs)
    dlov[k] <- calc_dlov(k, w, e, segs)
    drov[k] <- calc_drov(k, w, e, segs)
    F[k] <- .Machine$integer.max #as close to infinite as is comes
    if (k == 1){
        next #just need to initialize F[k] and pointer here
    }
    print(paste0("k ",k))
    for(j in 1:(k-1)) { #j+1
        dstar <- calc_ds((j+1), k, w, e, segs)
        Dtmp <- scoref(j, k, dl, dle, dlov, drov, dstar)
        D <- e(k-j) - 2*Dtmp

        bla <- c(bla, paste0("j:",j," k:",k," dstar:",dstar))
        blubb <- c(blubb, paste0("j:",j," k:",k," Dtmp:",Dtmp))
        if (D<0){
            blibb <- c(blibb, paste0("j:",j," k:",k," D:",D))
        }

        ##we don't want to store D, so we keep the pointer
        if( (F[j] + D) < F[k] ) { # Da hats was, F[k] is hier fast immer groesser als F[j] + D
            #print(paste0("FOUND NEW MINIMA ",(F[j] + D)," lower than ",F[k]," for k ",k," and j ",j))
            F[k] = F[j] + D
            ptr[k] = j
        }
    }
}

subset(bla, !grepl('dstar:0',bla))
subset(blubb, !grepl('Dtmp:0',blubb))
subset(blibb, !grepl('D:0',blibb))
##Backtrace Kette nach von letztem j nach 0
##
ptr
for (k in n:2){
    if(!is.na(ptr[k])){
        print(paste0("k ",k," : ",ptr[k]))
        #print(ptr[[k]])
        k = ptr[k]
    }
}
segs

#plot all
#
library(scales)
height <- max(dl,dle,dlov,drov)
Fs <- rescale(F, to = c(0,height))
segs
xlim <- segs
sep <- 0.5
col="lightgrey"
border="black"
if (is(xlim, "IntegerRanges")){
    xlim <- c(min(start(xlim)), max(end(xlim)))
}
#bins <- disjointBins(IRanges(start(x), end(x) + 1))
plot.new()
segnr <- length(subset(start(segs), start(segs) == 1))
plot.window(xlim, c(0, (segnr)*(height + sep)))
nr <- 0
#par(mfrow=c(2,1))
for (i in 1:length(start(segs))){
    if (start(segs)[i] == 1){
        nr <- nr + 1
    }
    ybottom <- nr*(sep + height) - height
    rect(start(segs)[i], ybottom, end(segs)[i], ybottom + height, col=col, border=border)
}
title("compare")
axis(1)

lines(dl, type = "o", col = "blue")
lines(dle, type = "o", col = "red")
lines(dlov, type = "o", col = "magenta")
lines(drov, type = "o", col = "green")
lines(Fs,type = "o", col = "black")




####MANUAL D berechnen

segs <- simulate_ranges(10,1,3,FALSE)
n <- max(end(segs)) # for testing

w <- function(m){return(1/m)}
e <- function(width){ return(width^2/2)}

dl <- numeric()
dle <- numeric()
dlov <- numeric()
drov <- numeric()

for (k in 1:n){
    dl[k] <- calc_dl(k, w, e, segs)
    dle[k] <- calc_dle(k, w, e, segs)
    dlov[k] <- calc_dlov(k, w, e, segs)
    drov[k] <- calc_drov(k, w, e, segs)
}

D <- matrix(nrow = n, ncol = n)

for (k in 1:n){
    if (k == 1){
        D[,k] <- 1
        next #just need to initialize D[,k] here
    }
    print(paste0("k ",k))

    for(j in 1:(k-1)) { #i+1
        D[j,k] <- calc_DELTA(j,k,w,e,segs)
    }
}


#### VS AUTO
####
dl <- numeric()
dle <- numeric()
dlov <- numeric()
drov <- numeric()
F <- numeric()
ptr <- numeric(1)

bla <- character()
blubb <- character()
blibb <- character()

Dm <- matrix(nrow = n, ncol = n)

for (k in 1:n){
    dl[k] <- calc_dl(k, w, e, segs)
    dle[k] <- calc_dle(k, w, e, segs)
    dlov[k] <- calc_dlov(k, w, e, segs)
    drov[k] <- calc_drov(k, w, e, segs)
    F[k] <- .Machine$integer.max #as close to infinite as is comes
    if (k == 1){
        next #just need to initialize F[k] and pointer here
    }
    print(paste0("k ",k))
    for(j in 1:(k-1)) { #j+1
        dstar <- calc_ds((j+1), k, w, e, segs)
        Dtmp <- scoref(j, k, dl, dle, dlov, drov, dstar)
        Dm[j,k] <- e(k-j) - 2*Dtmp
    }
}


subset(bla, !grepl('dstar:0',bla))
subset(blubb, !grepl('Dtmp:0',blubb))
subset(blibb, !grepl('D:0',blibb))
##Backtrace Kette nach von letztem j nach 0
for (k in n:2){
    if(!is.na(ptr[k])){
        print(paste0("k ",k," : ",ptr[k]))
        #print(ptr[[k]])
        #k = ptr[k]
    }
}
segs

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
