#include "ConsSeg_types.h"
#include <Rcpp.h>
using namespace Rcpp;

//' default potential function, \code{L^2/2}
//' @param L interval length
//' @param n total sequence length, required for function signature,
//' but not used in this example
// [[Rcpp::export]]
long double aeh(int L, int n) {
  long double e =  L*1.0; // implicit cast
  return(e*e/2);
}


// TODO: use string to switch between internal functions?
//' get internal potential function
// [[Rcpp::export]]
XPtr<funcPtr> e_ptr() {
  return(XPtr<funcPtr>(new funcPtr(&aeh)));
}


// backtrace function
// [[Rcpp::export]]
NumericVector backtrace_c(NumericVector imax) {

  int end = imax.size();
  NumericVector ends (end);
  ends[0] = end;
  int cnt = 1; 
  while ( end>1 ) {
    end = imax[end-1];
    ends[cnt] = end;
    cnt++;
  }
  // cut and reverse
  NumericVector res (cnt);
  for ( int i=0; i<cnt; i++ ) 
    res[cnt-i-1] = ends[i];
  return res;

}


// TODO: don't return all internal values, and use more specific
// types for them (long double, int)
//' Calculate consensus segments from a list of segmentation breakpoints
//' @param b list of breakpoints of different segmentations
//' @param n total sequence length 
//' @param w weights vector, must sum up to 1
//' @param e potential function, taking one argument: the length \code{L}
//' of the evaluated interval
//' @param store for debugging: store and return all internal vectors
//'@export
// [[Rcpp::export]]
List consensus_c(List b, int n, NumericVector w, SEXP e,
		 bool store=0) {

  // get potential function
  XPtr<funcPtr> xpfun(e);
  funcPtr aeh = *xpfun;

  int M = b.length(); // number of input segmentations

  
  //  FILL UP INTERVAL BORDER LOOKUP TABLES

  NumericMatrix Blw(n, M); // todo: int matrix
  NumericMatrix Bup(n, M);
  NumericVector bq;
  int up;
  int lw;
  for ( int q=0; q<M; q++ ) {
    bq = b[q];
    int i = 0;
    while ( bq[i] <= n ) {
      lw = bq[i];
      up = bq[i+1]-1;
      if ( up>n ) up = n ;
      for ( int k=lw-1; k<up; k++) { // k in lw:up ) { 
	Blw(k,q) = lw;
	Bup(k,q) = up;
      }
      i++; 
    } 
  }

    
  //  RECURSION
    
  // initialize recursion vectors
  NumericVector F(n+1);   // TODO: C++ long double
  NumericVector dsm(n+1);
  NumericVector dsq(n+1);
  NumericVector dcd(n+1);
  NumericVector dcu(n+1);
  NumericVector ptr(n+1); // TODO: integer
    
  // store \Delta and \delta^* for debugging
  if ( store ) {
    NumericVector Dk(n+1);
    NumericVector Ds(n+1);
  }

  // initialize vectors
  F[0] = dsm[0] = dsq[0] = dcd[0] = dsq[0] = 0.0;
  std::fill( ptr.begin(), ptr.end(), 0.0 );

  // helper variables and Delta
  long double D = 0.0; 
  long double Dtmp = 0.0;
  long double dtmp = 0.0;
  long double dstar = 0.0;

  // the recursion
  // NOTE: using n+1/F(0) et al. for convenience
  // and comparability with R:
  // counter = position index = vector index
  for ( int k=1; k<=n; k++ ) { // in 1:n ) { 
    
    dsm[k] = dsm[k-1];
    dsq[k] = dsq[k-1];
    dcd[k] = 0;
    dcu[k] = 0; 
        
    for ( int m=0; m<M; m++ ) {
      if ( Bup(k-1,m) == k ) // \delta_<(k), start and end left of k 
	dsm[k] += w[m]*aeh(Bup(k-1,m)-Blw(k-1,m)+1, n);
      if ( Blw(k-1,m) == k ) // \delta_le(j), start left of j, to subtract
	dsq[k] += w[m]*aeh(Bup(k-1,m)-Blw(k-1,m)+1, n);
      if ( Bup(k-1,m) > k ) // \delta^\cap_<(k), left end to k
	dcd[k] += w[m]*aeh(k-Blw(k-1,m)+1, n);
      if ( Blw(k-1,m) < k ) // \delta^\cap_>(j+1), j+1 to right end
	dcu[k] += w[m]*aeh(Bup(k-1,m)-k+1, n);
    }

    /* scan interval = [j+1,k] for minimum j */
    for ( int j=0; j<k; j++ ) { //in 0:(k-1) ) { 
      // \delta^*: correct for segments that span [j+1,k]
      dstar = 0.0;
      for ( int m=0; m<M; m++ ) {
	if ( ( Blw(k-1,m) < j+1 ) && ( Bup(k-1,m) > k ) ) {
	  dtmp = aeh(Bup(k-1,m)-Blw(k-1,m)+1, n) + aeh(k-j, n) -
	    aeh(k-Blw(k-1,m)+1, n) - aeh(Bup(k-1,m)-j, n);
	  dstar += w[m]*dtmp;
	}
      }

      // calculate \Delta(j+1,k)
      Dtmp = dsm[k] - dsq[j] + dcd[k] + dcu[j+1] + dstar;
      D = aeh(k-j, n) - 2*Dtmp;

      // find F[k] = min \Delta(j+1,k) + F(j)
      // and store the j that delivered it
      if ( j==0 ) {
	F[k] = D;
      } else if ( F[j]+D < F[k] ) {
	F[k] = F[j]+D;
	ptr[k] = j;
	//if ( store ) { // store for debugging or plots
	//  Dk[k] = D;
	//  Ds[k] = dstar;
	//}
      }
    }
  }
    
  // BACKTRACE
  NumericVector bp;
  bp = backtrace_c(ptr+1);
  bp = bp[Rcpp::Range(0, bp.length()-2)];
  // TODO: only add  matrices with option
  return List::create(Named("breakpoints") = bp,
		      Named("ptr") = ptr,
		      Named("F") = F[Rcpp::Range(1,n)],
		      Named("dsm") = dsm[Rcpp::Range(1,n)],
		      Named("dsq") = dsq[Rcpp::Range(1,n)],
		      Named("dcu") = dcu[Rcpp::Range(1,n)],
		      Named("dcd") = dcd[Rcpp::Range(1,n)],
		      Named("Bup") = Bup,
		      Named("Blw") = Blw);  
}

