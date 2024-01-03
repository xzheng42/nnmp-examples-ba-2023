#include <iostream>
#include <RcppArmadillo.h>
#include <math.h>
#include "utils.h"

using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//////////////////////////////////////////////////////////////////////////
// Logit Gaussian Weights for full likelihood inference
/////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List logitGausWeight1_cpp(const int& nne,
                          const arma::mat& dist,
                          const double& zeta,
                          const arma::colvec& mu,
                          const arma::colvec& ka){
  
  int data_size = dist.n_rows;
  arma::colvec icutoff, ilogit_cutoff, idist, irho, ilog_rho, iweight, iprobs;
  arma::mat weights(data_size, nne);
  arma::mat cutoff(data_size, nne + 1);
  weights.fill(NA_REAL);
  cutoff.fill(NA_REAL);
  
  for (int i = 2; i < nne; ++i) {
    idist = arma::trans(dist(i, arma::span(0, i - 1)));
    ilog_rho = -idist / zeta;
    irho = arma::exp(ilog_rho - ilog_rho.max());
    icutoff.set_size(i + 1);
    icutoff.fill(0);
    icutoff.tail(i) = cumsum_cpp(irho / sum(irho));
    ilogit_cutoff = arma::log(icutoff / (1.0 - icutoff));
    iprobs = arma::normcdf(ilogit_cutoff, mu(i-2), ka(i-2));
    iweight = arma::diff(iprobs);
    weights(i, arma::span(0, i - 1)) = iweight.t();
    cutoff(i, arma::span(0, i)) = icutoff.t();
  }
  
  icutoff.set_size(nne + 1);
  for (int i = nne; i < data_size; ++i) {
    idist = dist.row(i).t();
    ilog_rho = -idist / zeta;
    ilog_rho -= max(ilog_rho);
    irho = arma::exp(ilog_rho);
    icutoff.fill(0);
    icutoff.tail(nne) = cumsum_cpp(irho / sum(irho));
    ilogit_cutoff = arma::log(icutoff / (1.0 - icutoff));
    iprobs = arma::normcdf(ilogit_cutoff, mu(i-2), ka(i-2));
    iweight = arma::diff(iprobs);
    weights.row(i) = iweight.t();
    cutoff.row(i) = icutoff.t();
  }
  
  return List::create(Named("weights") = weights, Named("cutoff") = cutoff);
  
}

//////////////////////////////////////////////////////////////////////////
// Logit Gaussian Weights for conditional likelihood inference 
/////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List logitGausWeight2_cpp(const int& nne,
                          const arma::mat& dist,
                          const double& zeta,
                          const arma::colvec& mu,
                          const arma::colvec& ka){
  
  int data_size = dist.n_rows;
  arma::colvec icutoff(nne + 1), ilogit_cutoff, idist, irho, ilog_rho, iweight, iprobs;
  arma::mat weights(data_size, nne);
  arma::mat cutoff(data_size, nne + 1);
  
  for (int i = 0; i < data_size; ++i) {
    idist = dist.row(i).t();
    ilog_rho = -idist / zeta;
    irho = arma::exp(ilog_rho - ilog_rho.max());
    icutoff.fill(0);
    icutoff.tail(nne) = cumsum_cpp(irho / sum(irho));
    ilogit_cutoff = arma::log(icutoff / (1.0 - icutoff));
    iprobs = arma::normcdf(ilogit_cutoff, mu(i), ka(i));
    iweight = arma::diff(iprobs);
    weights.row(i) = iweight.t();
    cutoff.row(i) = icutoff.t();
  }
  
  return List::create(Named("weights") = weights, Named("cutoff") = cutoff);
  
}