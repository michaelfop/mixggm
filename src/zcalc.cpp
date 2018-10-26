#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
Rcpp::List zcalc(arma::mat denspro)
{
  int N = denspro.n_rows;
  int K = denspro.n_cols;
  arma::uvec kappa = linspace<uvec>(0, K-1, K);
  double llk = 0.0;
  arma::mat z(N,K);

  for ( int i = 0; i < N; i++ ) {
    double m = max( denspro.row(i) );
    double loghood = m + std::log( sum( exp(denspro.row(i) - m) ) );
    llk += loghood;
    z.row(i) = exp( denspro.row(i) - loghood );
  }

  return Rcpp::List::create( Named("z") = z,
                             Named("loglik") = llk );
}
