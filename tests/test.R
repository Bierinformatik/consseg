###TEST AREA
library(ConsSeg)
data(primseg436_sset)
test <- sset
plotSegmentation(tset, cset, test, cex=.5, lwd=2)

SO <- extract_segments(test)
S <- collapse_segments(test,1)


for (j in 7429:7462){
  for (i in 1:3){
    if (j %in% segs["start"] | j in segs["end"]){
      SM[i,j] = w  # w marks existence of boundary and adds weight of that boundary
    }
  }
}

