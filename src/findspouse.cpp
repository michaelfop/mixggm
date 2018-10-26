#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
Rcpp::List findspouse(arma::umat graph)
{
  int V = graph.n_rows;
  Rcpp::List spo(V), spo2(V), nspo(V);
  arma::uvec numspo(V);

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
  arma::uvec nontrivial = find(numspo != 0);

  return Rcpp::List::create( Named("spo") = spo,
                             Named("spo2") = spo2,
                             Named("nspo") = nspo,
                             Named("numspo") = numspo,
                             Named("nontrivial") = nontrivial );
}
