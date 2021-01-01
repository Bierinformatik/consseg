// signature of potential functions
// TODO: rcpp export causes compilation error in sourcCpp,
// but missing causes missing typedef error in R CMD INSTALL!
typedef long double (*funcPtr)(int);

