#include <iostream>
#include <RcppArmadillo.h>
#include <math.h>
#include <truncnorm.h>
#include "utils.h"
#include "copulas.h"

using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

//////////////////////////////////////////////////////////////////////////
// Compute copula NNMP conditional cdf
//////////////////////////////////////////////////////////////////////////
typedef double (*margfunc)(double, double, double);

arma::cube conCopCDF(const arma::mat& grid_ne_obs,
                     const arma::mat& grid_ne_dist,
                     const int& cop_family,
                     const arma::cube& cop_param,
                     const arma::mat& marg_param,
                     const arma::colvec& zeta,
                     const arma::mat& ga,
                     const arma::colvec kasq,
                     const arma::mat& DD,
                     const arma::colvec& probs,
                     margfunc pMarg) {

  int grid_size = grid_ne_obs.n_rows;
  int sample_size = zeta.n_rows;
  int nne = grid_ne_obs.n_cols;
  int probs_size = probs.n_rows;

  arma::colvec ine_dist, ine_obs, iter_log_cutoff_rho, iter_cutoff_rho, iter_weight, 
               iter_logit_cutoff, iter_cdf;
  arma::colvec iter_log_comp_concdf(nne);
  arma::colvec iter_cutoff(nne + 1);
  arma::rowvec iDD;

  double iter_marg_param1, iter_marg_param2, iter_cop_param, 
         iter_zeta, iter_logit_mu, iter_ka, pp1, pp2_k, ine_k;
  double pmax = 1 - 1e-10;
  double pmin = 1e-10;

  arma::cube concdf(grid_size, probs_size, sample_size);

  for (int i = 0; i < grid_size; ++i) {

    ine_dist = grid_ne_dist.row(i).t();
    ine_obs = grid_ne_obs.row(i).t();
    iDD = DD.row(i);

    for (int iter = 0; iter < sample_size; ++iter) {

      iter_zeta = zeta(iter);
      iter_log_cutoff_rho = -ine_dist / iter_zeta;
      iter_cutoff_rho = arma::exp(iter_log_cutoff_rho - iter_log_cutoff_rho.max());
      iter_cutoff.fill(0);
      iter_cutoff.tail(nne) = cumsum_cpp(iter_cutoff_rho / sum(iter_cutoff_rho));
      iter_logit_cutoff = arma::log(iter_cutoff / (1.0 - iter_cutoff));
      iter_logit_mu = arma::as_scalar(iDD * ga.col(iter));
      iter_ka = sqrt(kasq(iter));
      iter_cdf = arma::normcdf(iter_logit_cutoff, iter_logit_mu, iter_ka);
      iter_weight = arma::diff(iter_cdf);

      iter_marg_param1 = marg_param(iter, 0);
      iter_marg_param2 = marg_param(iter, 1);

      for (int j = 0; j < probs_size; ++j) {

        pp1 = probs(j);

        for (int k = 0; k < nne; ++k) {

          ine_k = ine_obs(k);
          pp2_k = pMarg(ine_k, iter_marg_param1, iter_marg_param2);
          if (pp2_k < pmin) {
            pp2_k = pmin;
          } else if (pp2_k > pmax) {
            pp2_k = pmax;
          }
          
          iter_cop_param = cop_param(i, k, iter);
          iter_log_comp_concdf(k) = pConCop(pp1, pp2_k, cop_family, iter_cop_param, true);

        }

        concdf(i, j, iter) = sum(arma::exp(arma::log(iter_weight) + iter_log_comp_concdf));

      }
    }
  }

  return concdf;

}

// [[Rcpp::export]]
arma::cube cdfCopNNMP_simple(const arma::mat& grid_ne_obs,
                             const arma::mat& grid_ne_dist,
                             const int& cop_family,
                             const arma::cube& cop_param,
                             const int& marg_family,
                             const arma::mat& marg_param,
                             const arma::colvec& zeta,
                             const arma::mat& ga,
                             const arma::colvec kasq,
                             const arma::mat& DD,
                             const arma::colvec& probs){

  arma::cube cdf_out;
  
  if (marg_family == 1)
  {
    cdf_out = conCopCDF(grid_ne_obs, grid_ne_dist,
                        cop_family, cop_param, marg_param,
                        zeta, ga, kasq, DD, probs, pbeta);
  } 
  else if (marg_family == 2)
  {
    cdf_out = conCopCDF(grid_ne_obs, grid_ne_dist,
                        cop_family, cop_param, marg_param,
                        zeta, ga, kasq, DD, probs, pgamma);
  } 
  else {
    stop("error: this family is currently not supported.");
  }

  return cdf_out;

}