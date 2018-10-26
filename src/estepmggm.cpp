#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List estepmggm(arma::mat x,
                     arma::mat mu,
                     arma::cube sigma,
                     arma::vec pro)
{
  int N = x.n_rows;
  int K = mu.n_rows;
  int V = mu.n_cols;
  double llk, loghood, m;

  arma::mat dens(N,K), denspro(N,K), z(N,K);
  for ( int k = 0; k < K; k++ ) {
    arma::mat ss = sigma.slice(k);
    arma::mat rooti = arma::trans( arma::inv( trimatu( arma::chol(sigma.slice(k)) ) ) );
    double rootisum = arma::sum( log(rooti.diag()) );
    double constants = -(static_cast<double>(V)/2.0) * log2pi;

    for ( int i = 0; i < N; i++ ) {
      arma::vec d = rooti * arma::trans( x.row(i) - mu.row(k)) ;
      dens(i,k) = constants - 0.5 * arma::sum(d%d) + rootisum;
      denspro(i,k) = dens(i,k) + std::log(pro(k));

      if ( k == (K-1) ) {
        double m = max( denspro.row(i) );
        loghood = m + std::log( sum( exp(denspro.row(i) - m) ) );
        llk += loghood;
        z.row(i) = exp( denspro.row(i) - loghood );
      }
    }
  }

  return Rcpp::List::create( Named("dens") = dens,
                             Named("denspro") = denspro,
                             Named("z") = z,
                             Named("loglik") = llk
  );
}
