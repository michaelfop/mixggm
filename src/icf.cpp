#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List icf ( arma::mat Sigma,
                 arma::mat S,
                 arma::umat graph,
                 double N,
                 double tol,
                 int maxit,
                 bool traceOut,
                 bool regularize,
                 double psi )
{
  int V = S.n_rows;
  arma::uvec var = linspace<uvec>(0, V-1, V);
  arma::mat omega = Sigma;
  Rcpp::List spo(V), spo2(V), nspo(V);
  arma::uvec numspo(V);

  // find spouses and other stuff -----------------------------------
  for ( int j = 0; j < V; j++ ) {
    arma::uvec set = find( graph.col(j) );
    spo(j) = set;
    numspo(j) = set.n_rows;
    arma::uvec set2 = find( set > j );
    set.elem(set2) = set.elem(set2) - 1;
    spo2(j) = set;
    set = find( graph.col(j) == 0 );
    set2 = find(set != j);
    set = set.elem(set2);
    nspo(j) = set;
  }
  arma::uvec nontriv = find(numspo != 0);


  // icf main -------------------------------------------------------
  // outputs and initialization
  double loglik, llkPrev;
  bool crit = true;
  double err;
  int it = 0;
  double val;
  double sign;
  arma::mat SigmaPrev = Sigma;
  log_det(val, sign, SigmaPrev);
  arma::mat invPrev = inv_sympd(SigmaPrev);
  if ( regularize ) {
    llkPrev = -(N + V + psi + 2.0)/2.0 * val -(N + V + psi + 1.0)/2.0 * trace( invPrev*S ) ;
  } else {
    llkPrev = -N/2.0 * val - N/2.0 * trace( invPrev*S );
  }

  while ( crit ) {
    for ( int i = 0; i < nontriv.size(); ++i ) {

      int j = nontriv(i);

      arma::uvec sp = spo(j);               // spouses
      arma::uvec ns = nspo(j);              // no spouses
      arma::uvec sp2 = spo2(j);
      arma::uvec sub = find(var != j);        // set [-j]
      arma::uvec h = find(var == j);          // set [j]

      arma::mat Sginv = inv_sympd( Sigma(sub,sub) );
      arma::mat YZT = S.cols(sub);              // S[,-j]
      YZT = YZT.row(j) * Sginv.cols(sp2);       // S[j,-j] %*% inv[,s]
      arma::mat ZZT = Sginv.rows(sp2) * S(sub,sub) * Sginv.cols(sp2);   // inv[s,] %*% S[-i,-i] %*% inv[,s]
      Sigma(h,sp) = YZT * inv_sympd(ZZT);
      Sigma(sp,h) = Sigma(h,sp).t();
      Sigma(h,ns).fill(0.0);
      Sigma(ns,h).fill(0.0);
      arma::mat lambda = S(h,h) - Sigma(h,sp)*YZT.t();
      Sigma(h,h) = lambda + Sigma(h,sp) * Sginv(sp2,sp2) * Sigma(sp,h);
    }
    //------------------------------------------------------------------------

    // loglikelihood---------------------------------------------------------
    // -N/2*log( det(Sigma) ) - N/2*sum( diag(solve(Sigma) %*% S) )
    omega = inv_sympd(Sigma);
    log_det(val, sign, Sigma);
    if ( regularize ) {
      loglik = -(N + V + psi + 2.0)/2.0 * val -(N + V + psi + 1.0)/2.0 * trace( omega*S );
    } else {
      loglik = -N/2.0 * val - N/2.0 * trace( omega*S );
    }
    //-----------------------------------------------------------------------

    // check
    err = std::abs(loglik - llkPrev) / (1 + std::abs(loglik));
    llkPrev = loglik;
    SigmaPrev = Sigma;
    it++;
    crit = ( (err > tol) & (it < maxit) );
    if ( traceOut ) Rcpp::Rcout << loglik << std::endl;
  }

  return Rcpp::List::create( Named("sigma") = Sigma,
                             Named("omega") = omega,
                             Named("loglik") = loglik,
                             Named("it") = it,
                             Named("err") = err );
}
