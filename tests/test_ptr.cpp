#include <Rcpp.h>
using namespace Rcpp;

typedef long double (*funcPtr)(int L);

// [[Rcpp::export]]
long double aeh(int L) {
  return (L*L)/2;
}

// [[Rcpp::export]]
double callViaXPtr(int L, SEXP xpsexp) {
    XPtr<funcPtr> xpfun(xpsexp);
    funcPtr fun = *xpfun;
    long double y = fun(L);
    return (y);
}

// [[Rcpp::export]]
XPtr<funcPtr> putFunPtrInXPtr(std::string fstr) {
  return(XPtr<funcPtr>(new funcPtr(&aeh)));
}
