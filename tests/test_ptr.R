

## PASS COMIPLED FUNCTION POINTER TO Rcpp FUNCTION
## see https://gallery.rcpp.org/articles/passing-cpp-function-pointers/
## also see https://gallery.rcpp.org/articles/passing-cpp-function-pointers-rcppxptrutils/


library(Rcpp)
library(RcppXPtrUtils)
sourceCpp("test_ptr.cpp")

## user-supplied function
## NOTE: 3.0 is required, otherwise 21 is returned!
cpp <- "long double my_aeh(int L) { return (L*L*L)/3.0; }"
ptr <- cppXPtr(cpp)
checkXPtr(ptr, type="long double", args=c("int")) # returns silently
callViaXPtr(4, ptr)

cpp <- "long double my_aeh(int L) { return std::pow(L,3)/3; }"
ptr <- cppXPtr(cpp)
checkXPtr(ptr, type="long double", args=c("int")) # returns silently
callViaXPtr(4, ptr)

cpp <- "long double my_aeh(int L) { return exp(L/2)-1; }"
ptr <- cppXPtr(cpp)
checkXPtr(ptr, type="long double", args=c("int")) # returns silently
callViaXPtr(4, ptr)

## pre-compiled function
ptr <- putFunPtrInXPtr("aeh")
callViaXPtr(4, ptr)

## error check
cpp <- "double my_aeh(int L) { return (exp(L/2)-1); }"
ptr <- cppXPtr(cpp)
tmp <- try(checkXPtr(ptr, type="long double", args=c("int")))
if ( !"try-error"%in%class(tmp) ) {
    callViaXPtr(4, ptr)
}

cpp <- "long double my_aeh(int L) { return (exp(L/2)-1); }"
ptr <- cppXPtr(cpp)
tmp <- try(checkXPtr(ptr, type="long double", args=c("int")))
if ( !"try-error"%in%class(tmp) ) {
    callViaXPtr(4, ptr)
}


