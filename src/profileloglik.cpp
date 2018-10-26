#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List profileloglik(arma::mat Sigma, arma::mat S, int N)
{
  int V = Sigma.n_rows;

  arma::mat inv = inv_sympd(Sigma);
  double val;
  double sign;
  log_det(val, sign, Sigma);
  val = -N/2.0 * val - N/2.0 * trace( inv*S ) ;

  return Rcpp::List::create( Named("loglik") = val,
                             Named("omega") =  inv);
}
