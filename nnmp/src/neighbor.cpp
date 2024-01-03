#include <iostream>
#include <RcppArmadillo.h>
#include <math.h>

using namespace std;
using namespace Rcpp;

//////////////////////////////////////////////////////////////////////////
// Get Euclidean distance between two locations
//////////////////////////////////////////////////////////////////////////
double getED(const arma::rowvec& loc1, 
             const arma::rowvec& loc2) {
  
  double dd = sqrt(sum(arma::pow(loc1 - loc2, 2)));
  return dd;
  
}

//////////////////////////////////////////////////////////////////////////
// Get neighbor information of the reference locations
/////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List refNe(const int& nne, 
           const arma::mat& coords, 
           const arma::colvec& obs) {
  
  int nn = obs.n_rows;
  arma::mat ne_dist(nn, nne), ne_idx(nn, nne), ne_obs(nn, nne);
  arma::rowvec iloc(2), jloc(2), idist;
  arma::uvec ine_idx;
  ne_dist.fill(NA_REAL);
  ne_idx.fill(NA_REAL);
  ne_obs.fill(NA_REAL);
  
  for (int i = 1; i < nne; ++i) {
    
    iloc = coords.row(i);
    
    idist.set_size(i);
    
    for (int j = 0; j < i; ++j) {
      
      jloc = coords.row(j);
      idist(j) = getED(iloc, jloc);
      
    }
    
    ine_idx = arma::sort_index(idist);
    ne_dist(i, arma::span(0, i-1)) = idist(ine_idx).t();
    ne_idx(i, arma::span(0, i-1)) = arma::conv_to<arma::rowvec>::from(ine_idx) + 1;
    ne_obs(i, arma::span(0, i-1)) = obs(ine_idx).t();
    
  }
  
  for (int i = nne; i < nn; ++i) {
    
    iloc = coords.row(i);
    idist.set_size(i);
    
    for (int j = 0; j < i; ++j) {
      
      jloc = coords.row(j);
      idist(j) = getED(iloc, jloc);
      
    }
    
    ine_idx = arma::sort_index(idist);
    ne_dist.row(i) = idist(ine_idx.head(nne)).t();
    ne_idx.row(i) = arma::conv_to<arma::rowvec>::from(ine_idx.head(nne)) + 1;
    ne_obs.row(i) = obs(ine_idx.head(nne)).t();
    
  }
  
  return List::create(Named("ne_dist") = ne_dist, 
                      Named("ne_index") = ne_idx, 
                      Named("ne_obs") = ne_obs);
  
}



//////////////////////////////////////////////////////////////////////////
// Get neighbor information of the non-reference locations
/////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List gridNe(const int& nne, const arma::mat& coords, const arma::mat& grid, const arma::colvec& obs) {
  
  int nn = obs.n_rows;
  int mm = grid.n_rows;
  arma::mat ne_dist(mm, nne), ne_idx(mm, nne), ne_obs(mm, nne);
  arma::rowvec iloc(2), jloc(2), idist(nn);
  arma::uvec ine_idx(nn);
  
  for (int i = 0; i < mm; ++i) {
    
    iloc = grid.row(i);
    
    for (int j = 0; j < nn; ++j) {
      
      jloc = coords.row(j);
      idist(j) = getED(iloc, jloc);
      
    }
    
    ine_idx = arma::sort_index(idist);
    ne_dist.row(i) = idist(ine_idx.head(nne)).t();
    ne_idx.row(i) = arma::conv_to<arma::rowvec>::from(ine_idx.head(nne)) + 1;
    ne_obs.row(i) = obs(ine_idx.head(nne)).t();
    
  }
  
  return List::create(Named("ne_dist") = ne_dist, 
                      Named("ne_index") = ne_idx, 
                      Named("ne_obs") = ne_obs);
  
}


