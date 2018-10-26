#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
using namespace std;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List conggm(arma::mat S,
                  arma::umat graph,
                  int N,
                  double tol,
                  int maxit,
                  bool traceOut)
{
  int V = S.n_rows;
  arma::mat omega(V,V), sigma(V,V);
  arma::uvec var = linspace<uvec>(0, V-1, V);
  arma::mat W = S;          // W is Sigma with this notation
  arma::mat Wnoj(V-1,V-1);
  arma::mat Snoj(V-1,1);
  arma::vec w(V-1);
  Rcpp::List neigh(V);

  for ( int j = 0; j < V; j++ ) {
    arma::uvec neighj = find( graph.col(j) );
    arma::uvec set = find( neighj > j );
    neighj.elem(set) = neighj.elem(set) - 1;
    neigh(j) = neighj;
  }

  // algorithm for estimation ---------------------------------------
  double loglik;
  bool crit = true;
  double err;
  int it = 0;
  double val;
  double sign;
  arma::mat WPrev = W;
  log_det(val, sign, WPrev);
  arma::mat invPrev = inv_sympd(WPrev);
  double llkPrev = -N/2.0 * val - N/2.0 * trace( invPrev*S );

  while ( crit ) {

      for ( int j = 0; j < V; j++ ) {
      arma::uvec noj = find(var != j);          // set [-j]
      arma::uvec sij = find(var == j);          // set [j]
      arma::uvec neighj = neigh(j);             // neighbours of node j
      arma::vec beta = zeros<vec>(V-1);

      if ( neighj.size() < 1 ) {
        w = beta;
      } else {
        Wnoj = W(noj,noj);
        Snoj = S(noj,sij);
        beta(neighj) = solve(Wnoj(neighj,neighj), Snoj.rows(neighj));
        w = Wnoj*beta;
      }
      W(noj, sij) = w;
      W(sij, noj) = w.t();
    }

    // loglikelihood ------------------------------------------------
    // -N/2*log( det(Sigma) ) - N/2*sum( diag(solve(Sigma) %*% S) )
    log_det(val, sign, W);
    omega = inv_sympd(W);
    loglik = -N/2.0 * val - N/2.0 * trace( omega*S );

    // check
    err = std::abs(loglik - llkPrev) / (1 + std::abs(loglik));
    llkPrev = loglik;
    WPrev = W;
    it++;
    crit = ( (err > tol) & (it < maxit) );
    if ( traceOut ) Rcpp::Rcout << loglik << std::endl;
  }

  arma::vec d = omega.diag();
  omega = omega % graph;
  omega.diag() = d;

  return Rcpp::List::create( Named("sigma") = W,
                             Named("omega") = omega,
                             Named("loglik") = loglik,
                             Named("it") = it,
                             Named("err") = err );
}
