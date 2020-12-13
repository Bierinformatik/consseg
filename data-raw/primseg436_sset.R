## code to prepare `primseg436_sset` dataset goes here
library(segmenTier)
### DYNAMIC PROGRAMMING BASED SEGMENTATION OF A CLUSTERING
### implemented by Rainer Machne, hopefully
### as conceived by Peter F. Stadler
data(primseq436)
trafo <- "raw"     # transformation of the raw data
use.fft <- TRUE    # cluster discrete Fourier transform of data?
use.snr <- TRUE    # use DFT scaling (SNR is described as
# relative amplitude scaling in Machne&Murray, PLoS ONE 2012)
dft.range <- 1:7   # range of DFT to use for clustering
dc.trafo <- "ash"  # transformation of the first (DC) component of the DFT
# NOTE: add component 1 (DC) to DFT range to use
low.thresh <- -Inf # minimal total signal (DC component of DFT if use.fft)
### CLUSTERING PARAMETERS
K <- c(12)         # cluster number K; multiple allowed; specifically, note
# that k-means has a random effect at initialization
# and replicates of the same K can  yield different
# results for otherwise
nui.thresh <- 0.6  # threshold of position-cluster correlation below which
# the position will be assigned to the nuisance cluster
## k-means initialization
iter.max <- 100000 # max. iterations in kmeans
nstart <- 100      # number of initial configurations tested in kmeans
### SEGMENTATION PARAMETERS
## segmenTier parameters are handled via the settings function,
## where all parameters can be passed as vectors.
vary <- setVarySettings(
  E=c(1,3), # scale exponent of similarity matrices csim
  S="icor", # SCORING FUNCTIONS
  M=c(150), # scoring function minimal length penalty
  Mn=100,   # M for nuisance clusters
  nui=c(1,3)#-/+ correlation of nuisance cluster with others and itself
)
## PRE-PROCESS TIME SERIES FOR CLUSTERING
## take DFT and scale amplitudes, and
## select components of DFT
tset <- processTimeseries(ts=tsd, na2zero=TRUE,
                          trafo=trafo, dc.trafo=dc.trafo,
                          use.fft=use.fft, dft.range=dft.range,
                          use.snr=use.snr, low.thresh=low.thresh)
## CLUSTER PRE-PROCESSED TIME SERIES
set.seed(15) # stable kmeans clustering
cset <- clusterTimeseries(tset, K=K, iter.max=iter.max, nstart=nstart,
                          nui.thresh=nui.thresh)
### CALCULATE SEGMENTS FOR ALL CLUSTERINGS and
### FOR CHOSEN SEGMENTATION PARAMETERS
## NOTE, that the function optionally also stores the
## the scoring and backtracing matrices (see demo/segment_test.R)
sset <- segmentCluster.batch(cset, varySettings=vary,
                             id="mysegments",
                             type.name=c("E","M","nui"), # segment type names
                             verb=1, save.matrix=FALSE)

usethis::use_data(primseg436_sset, overwrite = TRUE)
