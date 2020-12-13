###TEST AREA
library(ConsSeg)
#Load Data
data(primseg436_sset)

#Extract starts and end of segments
SO <- extract_segments(sset)
#Create segment matrix for dyn-prog approach
SM <- matrixfy_segments(sset, 1)


###Manual
segs <- extract_segments(sset)
m <- length(segs)/2
n <- sset$N #Total length of input sequences

