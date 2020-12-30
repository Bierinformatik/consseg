


//' Calculate consensus segments from a list of segmentation breakpoints
//' @param b list of breakpoints of different segmentations
//' @param n total sequence length (\code{max(b)} if not provided)
//' @param w weight function, taking one argument: the index \code{m} of
//' the respective segmentation in the breakpoint list \code{b}
//' @param e potential function, taking one argument: the length \code{L}
//' of the evaluated interval
//' @param store.matrix for debugging: store and return all internal vectors
//'@export
// [[Rcpp::export]]
List consensus_c(List b, int n, std::string w, std::string aeh,
		 bool store_matrix=0) {

  // TODO: outside in R, add argument n if missing, add 1, n  to b list
  
  M = b.length(); // number of input segmentations

  
  //  FILL UP INTERVAL BORDER LOOKUP TABLES
  // TODO: make sure that this is size m and
  // position index is index-1
  NumericMatrix Blw(n, M);
  NumericMatrix Bup(n, M);
  for ( int q=0; q<M; q++ ) {
    int i = 1;
    while ( b[[q]][i] <= n ) {
      lw = b[[q]][i];
      up = b[[q]][i+1]-1;
      if ( up>n ) up = n ;
      for ( int k=lw-1; k<up; k++) { // k in lw:up ) { 
	Blw[k,q] = lw;
	Bup[k,q] = up;
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
  if ( store_matrix ) {
    NumericVector Dk(n+1);
    NumericVector Ds(n+1);
  }

  // initialize vectors
  F[0] = dsm[0] = dsq[0] = dcd[0] = dsq[0] = 0;
  
  // the recursion
  // NOTE: using n+1/F(0) et al. for convenience
  // and comparability with R:
  // counter = position index = vector index
  for ( int k=1; k<=n; n++ ) { // in 1:n ) { 
    
    dsm[k] = dsm[k-1];
    dsq[k] = dsq[k-1];
    dcd[k] = 0;
    dcu[k] = 0; 
        
    for ( int m=0; m<M; m++ ) {
      if ( Bup[k-1,m] == k ) // \delta_<(k), start and end left of k 
	dsm[k] = dsm[k] + w(m)*aeh(Bup[k-1,m]-Blw[k-1,m]+1);
      if ( Blw[k-1,m] == k ) // \delta_le(j), start left of j, to subtract
	dsq[k] = dsq[k] + w(m)*aeh(Bup[k-1,m]-Blw[k-1,m]+1);
      if ( Bup[k-1,m] > k ) // \delta^\cap_<(k), left end to k
	dcd[k] = dcd[k] + w(m)*aeh(k-Blw[k-1,m]+1);
      if ( Blw[k-1,m] < k ) // \delta^\cap_>(j+1), j+1 to right end
	dcu[k] = dcu[k] + w(m)*aeh(Bup[k-1,m]-k+1);
    }
        
    /* scan interval = [j+1,k] for minimum j */
    for ( int j=0; j<k-1; j++ ) { //in 0:(k-1) ) { 
      // \delta^*: correct for segments that span [j+1,k]
      long double dstar = 0;
      for ( int m=0; m<M; m++ ) 
	if ( ( Blw[k-1,m] < j+1 ) && ( Bup[k-1,m] > k ) ) {
	  dtmp = aeh(Bup[k-1,m]-Blw[k-1,m]+1) + aeh(k-j) -
	    aeh(k-Blw[k-1,m]+1) - aeh(Bup[k-1,m]-j);
	  dstar = dstar + w(m)*dtmp;
	}
      
      Dtmp = dsm[k] + dcd[k] + dcu[j+1] + dstar;
      if ( j>0 ) Dtmp = Dtmp - dsq[j];
        
      D = aeh(k-j) - 2*Dtmp;

      // find F[k] = min Delta(j+1,k) + F(j)
      // and store the j that delivered it
      if ( j==0 ) {
	F[k] = D;
      } else if ( F[j]+D < F[k] ) {
	F[k] = F[j]+D;
	ptr[k] = j;
	if ( store_matrix ) { // store for debugging or plots
	  Dk[k] = D;
	  Ds[k] = dstar;
	}
      }
    }
  }
    
  // TODO: BACKTRACE
  // TODO: add generated matrices
  return ptr;
}

