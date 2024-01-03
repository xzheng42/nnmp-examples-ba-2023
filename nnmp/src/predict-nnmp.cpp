#include <iostream>
#include <RcppArmadillo.h>
#include <math.h>
#include <RcppTN.h>
#include "utils.h"
#include "copulas.h"

using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo, RcppTN)]]

//////////////////////////////////////////////////////////////////////////
// Prediction function for Gaussian NNMPs
//////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List predGausNNMP(const int& nne, 
                  const arma::mat& XX,
                  const arma::mat& DD,
                  const arma::mat& bb,
                  const arma::mat& zz,
                  const arma::colvec& sigmasq, 
                  const arma::colvec& tausq, 
                  const arma::colvec& phi,
                  const arma::colvec& zeta,
                  const arma::mat& ga,
                  const arma::colvec& kasq,
                  const arma::mat& grid_ne_idx,
                  const arma::mat& grid_ne_dist,
                  const arma::colvec& probs,
                  const bool& verbose,
                  const int& nreport,
                  const bool sam = false){
  
  int grid_size = grid_ne_dist.n_rows;
  int sample_size = phi.n_rows;
  int probs_size = probs.n_rows;
  int label, iter_zz_ne_label;
  
  double iter_rho, iter_ine_dist, iter_zz_ne, iter_sigmasq, iter_phi, iter_zeta, 
         iter_mu, iter_logit_mu, iter_ka;
  
  arma::colvec label_vals = arma::linspace<arma::colvec>(0, nne-1, nne); 
  arma::colvec iter_bb, ine_dist, iter_cutoff_rho, iter_log_cutoff_rho, iter_logit_cutoff, iter_cdf, iter_weight;
  arma::colvec iter_cutoff(nne + 1);
  arma::rowvec iXX, iDD;
  
  arma::colvec ipred_zz(sample_size);
  arma::colvec pred_zz_mean(grid_size);
  arma::mat pred_zz_qq(probs_size, grid_size);
  
  arma::colvec ipred_obs(sample_size);
  arma::colvec pred_obs_mean(grid_size);
  arma::mat pred_obs_qq(probs_size, grid_size);
  
  arma::mat pred_zz_sam(grid_size, sample_size);
  arma::mat pred_obs_sam(grid_size, sample_size);
  
  Rprintf("----------------------------------------\n");
  Rprintf("\t  Predicting\n");
  Rprintf("----------------------------------------\n");
  
  int icount = 0;
  
  for (int i = 0; i < grid_size; ++i) {
    
    iXX = XX.row(i);
    iDD = DD.row(i);
    ine_dist = grid_ne_dist.row(i).t();
    
    for (int iter = 0; iter < sample_size; ++iter) {
      
      //Compute the weights
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
      
      // Sample the configuration variable
      label = arma::as_scalar(armaIntSample(label_vals, 1, 1, iter_weight));
      
      // Generate the latent z
      iter_ine_dist = ine_dist(label);
      iter_phi = phi(iter);
      iter_rho = exp(-iter_ine_dist / iter_phi);
      iter_zz_ne_label = grid_ne_idx(i, label);
      iter_zz_ne = zz(iter_zz_ne_label - 1, iter);
      iter_sigmasq = sigmasq(iter);
      ipred_zz(iter) = R::rnorm(iter_rho * iter_zz_ne, sqrt(iter_sigmasq * (1 - pow(iter_rho, 2))));
      
      // Generate the outcome y
      iter_bb = bb.col(iter);
      iter_mu = arma::as_scalar(iXX * iter_bb) + ipred_zz(iter);
      ipred_obs(iter) = R::rnorm(iter_mu, sqrt(tausq(iter)));
      
      if (sam) {
        pred_zz_sam(i, iter) = ipred_zz(iter);
        pred_obs_sam(i, iter) = ipred_obs(iter);
      }
      
    }
    
    // Compute mean and quantiles
    pred_zz_mean(i) = mean(ipred_zz);
    pred_zz_qq.col(i) = arma::quantile(ipred_zz, probs);
    pred_obs_mean(i) = mean(ipred_obs);
    pred_obs_qq.col(i) = arma::quantile(ipred_obs, probs);
    
    icount++;
    if (verbose) {
      if (icount == nreport) {
        Rprintf("  Locations: %i/%i, %3.2f%%\n", i+1, grid_size, 100.0*(i+1)/grid_size);
        icount = 0;
      }
    }
    
    
  }
  
  
  if (sam) {
    return(List::create(Named("zz_mu") = pred_zz_mean, 
                        Named("zz_qq") = pred_zz_qq,
                        Named("obs_mu") = pred_obs_mean, 
                        Named("obs_qq") = pred_obs_qq,
                        Named("zz_sam") = pred_zz_sam,
                        Named("obs_sam") = pred_obs_sam));
  } else {
    return(List::create(Named("zz_mu") = pred_zz_mean, 
                        Named("zz_qq") = pred_zz_qq, 
                        Named("obs_mu") = pred_obs_mean, 
                        Named("obs_qq") = pred_obs_qq));
  }  
  
  
}

//////////////////////////////////////////////////////////////////////////
// Prediction function for skew-Gaussian NNMPs
//////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List predSkewGausNNMP_simple(const arma::colvec& obs,
                             const int& nne, 
                             const arma::mat& DD,
                             const arma::colvec& la,
                             const arma::colvec& sigmasq, 
                             const arma::colvec phi,
                             const arma::colvec zeta,
                             const arma::mat& ga,
                             const arma::colvec kasq,
                             const arma::mat grid_ne_idx, 
                             const arma::mat grid_ne_dist, 
                             const arma::colvec& probs,
                             const bool& verbose,
                             const int& nreport,
                             const bool sam = false){
  
  
  int grid_size = grid_ne_dist.n_rows;
  int sample_size = phi.n_rows;
  int probs_size = probs.n_rows;
  int label, iter_ne_label;
  
  double iter_ine_dist, iter_rho, iter_ne, iter_sigmasq, iter_phi, iter_zeta, iter_tau2;
  double iter_xi, iter_la, z0, iter_mu, iter_logit_mu, iter_ka;  
  
  arma::colvec label_vals = arma::linspace<arma::colvec>(0, nne-1, nne); 
  arma::colvec ine_dist, iter_log_cutoff_rho, iter_cutoff_rho, iter_logit_cutoff, iter_cdf, iter_weight;
  arma::mat pred_obs(probs_size, grid_size);
  arma::colvec ipred_obs(sample_size);
  arma::colvec pred_obs_mean(grid_size);
  arma::colvec iter_cutoff(nne + 1);
  arma::rowvec iDD;
  arma::mat pred_sam(grid_size, sample_size);
  
  Rprintf("----------------------------------------\n");
  Rprintf("\t  Predicting\n");
  Rprintf("----------------------------------------\n");
  
  int icount = 0;
  
  for (int ii = 0; ii < grid_size; ++ii) {
    
    iDD = DD.row(ii);
    ine_dist = grid_ne_dist.row(ii).t();
    
    for (int iter = 0; iter < sample_size; ++iter) {
      
      // Compuate weights and generate a configuration variable
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
      label = arma::as_scalar(armaIntSample(label_vals, 1, 1, iter_weight));
      iter_ine_dist = ine_dist(label);
      
      // Generate y0
      iter_phi = phi(iter);
      iter_rho = exp(-iter_ine_dist / iter_phi);
      iter_sigmasq = sigmasq(iter);
      iter_la = la(iter);
      iter_tau2 = iter_sigmasq / (iter_sigmasq + pow(iter_la, 2));
      iter_ne_label = grid_ne_idx(ii, label);
      iter_ne = obs(iter_ne_label - 1);
      iter_xi = iter_tau2 * iter_ne * iter_la / iter_sigmasq;
      z0 = RcppTN::rtn1(iter_xi, sqrt(iter_tau2), 0, arma::datum::inf);  
      iter_mu = (1.0 - iter_rho) * iter_la * z0 + iter_rho * iter_ne;
      ipred_obs(iter) = R::rnorm(iter_mu, sqrt(iter_sigmasq * (1 - iter_rho * iter_rho)));
      
      if (sam) {
        pred_sam(ii, iter) = ipred_obs(iter);
      }    
      
    }
    
    pred_obs.col(ii) = arma::quantile(ipred_obs, probs);
    pred_obs_mean(ii) = mean(ipred_obs);
    
    icount++;
    if (verbose) {
      if (icount == nreport) {
        Rprintf("  Locations: %i/%i, %3.2f%%\n", ii+1, grid_size, 100.0*(ii+1)/grid_size);
        icount = 0;
      }
    }
    
    
    
  }
  
  if (sam) {
    return(List::create(Named("obs_mu") = pred_obs_mean, 
                        Named("obs_qq") = pred_obs,
                        Named("sam") = pred_sam));
  } else {
    return(List::create(Named("obs_mu") = pred_obs_mean, 
                        Named("obs_qq") = pred_obs));
  }  
  
  
}


// [[Rcpp::export]]
List predSkewGausNNMP_covars_lapar(const arma::colvec& obs,
                                   const int& nne, 
                                   const arma::mat& ref_XX,
                                   const arma::mat& nonref_XX,
                                   const arma::colvec& ref_par_lables,
                                   const arma::colvec& nonref_par_labels,
                                   const arma::mat& DD,
                                   const arma::mat& bb,
                                   const arma::mat& la,
                                   const arma::colvec& sigmasq, 
                                   const arma::colvec& phi,
                                   const arma::colvec& zeta,
                                   const arma::mat& ga,
                                   const arma::colvec& kasq,
                                   const arma::mat& grid_ne_idx,
                                   const arma::mat& grid_ne_dist, 
                                   const arma::colvec& probs,
                                   const bool& verbose,
                                   const int& nreport,
                                   const bool sam = false){
  
  
  int grid_size = grid_ne_dist.n_rows;
  int sample_size = phi.n_rows;
  int probs_size = probs.n_rows;
  int label, ipar_label, iter_ine_idx;
  
  double iter_ine_dist, iter_rho, iter_ne, iter_sigma2, iter_phi, iter_zeta, iter_tau2;
  double iter_xi, iter_la, z0, iter_mu, iter_logit_mu, iter_ka, iter_nu, iter_nu_ne, iter_la_ne;  
  
  arma::colvec label_vals = arma::linspace<arma::colvec>(0, nne-1, nne); 
  arma::colvec ine_dist, iter_log_cutoff_rho, iter_cutoff_rho, iter_logit_cutoff, iter_cdf, iter_weight;
  arma::mat pred_obs(probs_size, grid_size);
  arma::colvec ipred_obs(sample_size);
  arma::colvec pred_obs_mean(grid_size);
  arma::colvec iter_cutoff(nne + 1);
  arma::rowvec iDD;
  arma::colvec iter_bb, inu, iter_ref_nu, ipar_la;
  arma::mat pred_sam(grid_size, sample_size);
  
  Rprintf("----------------------------------------\n");
  Rprintf("\t  Predicting\n");
  Rprintf("----------------------------------------\n");
  
  int icount = 0;
  
  for (int ii = 0; ii < grid_size; ++ii) {
    
    iDD = DD.row(ii);
    ine_dist = grid_ne_dist.row(ii).t();
    inu = arma::trans(nonref_XX.row(ii) * bb);
    ipar_label = nonref_par_labels(ii) - 1;
    ipar_la = la.row(ipar_label).t();
    
    for (int iter = 0; iter < sample_size; ++iter) {

      // Compuate weights and generate a configuration variable
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
      label = arma::as_scalar(armaIntSample(label_vals, 1, 1, iter_weight));
      iter_ine_dist = ine_dist(label);
      iter_ine_idx = grid_ne_idx(ii, label) - 1;

      // Generate y0
      iter_phi = phi(iter);
      iter_rho = exp(-iter_ine_dist / iter_phi);
      iter_sigma2 = sigmasq(iter);
      iter_nu = inu(iter);
      iter_la = ipar_la(iter);
      iter_la_ne = la(ref_par_lables(iter_ine_idx) - 1, iter);
      iter_tau2 = iter_sigma2 / (iter_sigma2 + pow(iter_la_ne, 2));
      iter_ne = obs(iter_ine_idx);
      iter_bb = bb.col(iter);
      iter_ref_nu = ref_XX * iter_bb;
      iter_nu_ne = iter_ref_nu(iter_ine_idx);
      iter_xi = iter_tau2 * (iter_ne - iter_nu_ne) * iter_la_ne / iter_sigma2;
      z0 = RcppTN::rtn1(iter_xi, sqrt(iter_tau2), 0, arma::datum::inf);  
      iter_mu = iter_nu + iter_la * z0 + iter_rho * (iter_ne - iter_nu_ne - iter_la_ne * z0);
      ipred_obs(iter) = R::rnorm(iter_mu, sqrt(iter_sigma2 * (1 - iter_rho * iter_rho)));

      if (sam) {
        pred_sam(ii, iter) = ipred_obs(iter);
      }

    }
    
    pred_obs.col(ii) = arma::quantile(ipred_obs, probs);
    pred_obs_mean(ii) = mean(ipred_obs);
    
    icount++;
    if (verbose) {
      if (icount == nreport) {
        Rprintf("  Locations: %i/%i, %3.2f%%\n", ii+1, grid_size, 100.0*(ii+1)/grid_size);
        icount = 0;
      }
    }
    
  }
  
  if (sam) {
    return(List::create(Named("obs_mu") = pred_obs_mean, 
                        Named("obs_qq") = pred_obs,
                        Named("sam") = pred_sam));
  } else {
    return(List::create(Named("obs_mu") = pred_obs_mean, 
                        Named("obs_qq") = pred_obs));
  }  
  
}

//////////////////////////////////////////////////////////////////////////
// Prediction functions for copula NNMPs with continuous marginals
//////////////////////////////////////////////////////////////////////////

typedef double (*margfunc)(double, double, double);

List predCop_simple(const arma::colvec& obs,
                    const int& cop_family,
                    const arma::cube& cop_param,
                    const arma::mat marg_param,
                    const int& nne,
                    const arma::mat& DD,
                    const arma::colvec zeta,
                    const arma::mat& ga,
                    const arma::colvec kasq,
                    const arma::mat grid_ne_idx,
                    const arma::mat grid_ne_dist,
                    const arma::colvec& probs,
                    const bool sam,
                    margfunc pMarg,
                    margfunc qMarg,
                    const bool& verbose,
                    const int& nreport){

  int grid_size = grid_ne_dist.n_rows;
  int sample_size = zeta.n_rows;
  int probs_size = probs.n_rows;
  int label, iter_ne_label;

  double iter_ne, iter_zeta, iter_logit_mu, iter_ka,
         iter_cop_param, iter_marg_param1, iter_marg_param2,
         iter_vv, iter_uu;

  arma::colvec label_vals = arma::linspace<arma::colvec>(0, nne - 1, nne);
  arma::colvec ipred_obs(sample_size);
  arma::colvec pred_obs_mean(grid_size);
  arma::mat pred_obs(probs_size, grid_size);
  arma::colvec ine_dist, iter_log_cutoff_rho, iter_cutoff_rho,
  iter_logit_cutoff, iter_cdf, iter_weight;
  arma::colvec iter_cutoff(nne + 1);
  arma::rowvec iDD;
  arma::mat pred_sam(grid_size, sample_size);
  
  Rprintf("----------------------------------------\n");
  Rprintf("\t  Predicting\n");
  Rprintf("----------------------------------------\n");
  
  int icount = 0;

  for (int ii = 0; ii < grid_size; ++ii) {

    iDD = DD.row(ii);
    ine_dist = grid_ne_dist.row(ii).t();

    for (int iter = 0; iter < sample_size; ++iter) {
      
      // Compute weights and generate label
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
      label = arma::as_scalar(armaIntSample(label_vals, 1, 1, iter_weight));
      
      // Predic y0
      iter_ne_label = grid_ne_idx(ii, label);
      iter_ne = obs(iter_ne_label - 1);
      iter_cop_param = cop_param(ii, label, iter);
      iter_marg_param1 = marg_param(iter, 0);
      iter_marg_param2 = marg_param(iter, 1);
      iter_vv = pMarg(iter_ne, iter_marg_param1, iter_marg_param2);
      iter_uu = rConCop(iter_vv, cop_family, iter_cop_param);
      ipred_obs(iter) = qMarg(iter_uu, iter_marg_param1, iter_marg_param2);
 
      if (sam) {
        pred_sam(ii, iter) = ipred_obs(iter);
      }

    }

    pred_obs.col(ii) = arma::quantile(ipred_obs, probs);
    pred_obs_mean(ii) = mean(ipred_obs);
    
    icount++;
    if (verbose) {
      if (icount == nreport) {
        Rprintf("  Locations: %i/%i, %3.2f%%\n", ii+1, grid_size, 100.0*(ii+1)/grid_size);
        icount = 0;
      }
    }


  }

  if (sam) {
    return(List::create(Named("obs_mu") = pred_obs_mean, 
                        Named("obs_qq") = pred_obs,
                        Named("sam") = pred_sam));
  } else {
    return(List::create(Named("obs_mu") = pred_obs_mean, 
                        Named("obs_qq") = pred_obs));
  }


}


List predCop_covars(const arma::colvec& obs,
                    const int& cop_family,
                    const arma::cube& cop_param,
                    const arma::mat ref_marg_param1,
                    const arma::mat ref_marg_param2,
                    const arma::mat nonref_marg_param1,
                    const arma::mat nonref_marg_param2,
                    const int& nne,
                    const arma::mat& DD,
                    const arma::colvec zeta,
                    const arma::mat& ga,
                    const arma::colvec kasq,
                    const arma::mat grid_ne_idx,
                    const arma::mat grid_ne_dist,
                    const arma::colvec& probs,
                    const bool sam,
                    margfunc pMarg,
                    margfunc qMarg,
                    const bool& verbose,
                    const int& nreport){

  int grid_size = grid_ne_dist.n_rows;
  int sample_size = zeta.n_rows;
  int probs_size = probs.n_rows;
  int label, iter_ne_label;

  double iter_ne, iter_zeta, iter_logit_mu, iter_ka,
         iter_cop_param, iter_marg_param1, iter_marg_param2,
         iter_ne_marg_param1, iter_ne_marg_param2,
         iter_vv, iter_uu;

  arma::colvec label_vals = arma::linspace<arma::colvec>(0, nne - 1, nne);
  arma::colvec ipred_obs(sample_size);
  arma::colvec pred_obs_mean(grid_size);
  arma::mat pred_obs(probs_size, grid_size);
  arma::colvec ine_dist, iter_log_cutoff_rho, iter_cutoff_rho,
  iter_logit_cutoff, iter_cdf, iter_weight;
  arma::colvec iter_cutoff(nne + 1);
  arma::rowvec iDD;
  arma::mat pred_sam(grid_size, sample_size);
  
  Rprintf("----------------------------------------\n");
  Rprintf("\t  Predicting\n");
  Rprintf("----------------------------------------\n");
  
  int icount = 0;

  for (int ii = 0; ii < grid_size; ++ii) {

    iDD = DD.row(ii);
    ine_dist = grid_ne_dist.row(ii).t();

    for (int iter = 0; iter < sample_size; ++iter) {
      
      // Compute weights and generate label
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
      label = arma::as_scalar(armaIntSample(label_vals, 1, 1, iter_weight));
      
      // Predic y0
      iter_ne_label = grid_ne_idx(ii, label);
      iter_ne = obs(iter_ne_label - 1);
      iter_cop_param = cop_param(ii, label, iter);
      iter_marg_param1 = nonref_marg_param1(ii, iter);
      iter_marg_param2 = nonref_marg_param2(ii, iter);
      iter_ne_marg_param1 = ref_marg_param1(iter_ne_label - 1, iter);
      iter_ne_marg_param2 = ref_marg_param2(iter_ne_label - 1, iter);
      iter_vv = pMarg(iter_ne, iter_ne_marg_param1, iter_ne_marg_param2);
      iter_uu = rConCop(iter_vv, cop_family, iter_cop_param);
      ipred_obs(iter) = qMarg(iter_uu, iter_marg_param1, iter_marg_param2);
 
      if (sam) {
        pred_sam(ii, iter) = ipred_obs(iter);
      }

    }

    pred_obs.col(ii) = arma::quantile(ipred_obs, probs);
    pred_obs_mean(ii) = mean(ipred_obs);
    
    icount++;
    if (verbose) {
      if (icount == nreport) {
        Rprintf("  Locations: %i/%i, %3.2f%%\n", ii+1, grid_size, 100.0*(ii+1)/grid_size);
        icount = 0;
      }
    }


  }

  if (sam) {
    return(List::create(Named("obs_mu") = pred_obs_mean, 
                        Named("obs_qq") = pred_obs,
                        Named("sam") = pred_sam));
  } else {
    return(List::create(Named("obs_mu") = pred_obs_mean, 
                        Named("obs_qq") = pred_obs));
  }


}


// [[Rcpp::export]]
List predCopNNMP_simple(const arma::colvec& obs,
                        const int& cop_family,
                        const arma::cube& cop_param,
                        const int& marg_family,
                        const arma::mat marg_param,
                        const int& nne,
                        const arma::mat& DD,
                        const arma::colvec zeta,
                        const arma::mat& ga,
                        const arma::colvec kasq,
                        const arma::mat grid_ne_idx,
                        const arma::mat grid_ne_dist,
                        const arma::colvec& probs,
                        const bool& verbose,
                        const int& nreport,
                        const bool sam = false){

  List pred_out;
  
  if (marg_family == 1)
  {
    pred_out = predCop_simple(obs, cop_family, cop_param, marg_param, nne, DD,
                              zeta, ga, kasq, grid_ne_idx, grid_ne_dist, probs,
                              sam, pbeta, qbeta, verbose, nreport);
  } 
  else if (marg_family == 2)
  {
    pred_out = predCop_simple(obs, cop_family, cop_param, marg_param, nne, DD,
                              zeta, ga, kasq, grid_ne_idx, grid_ne_dist, probs,
                              sam, pgamma, qgamma, verbose, nreport);
  } 
  else {
    stop("error: this family is currently not supported.");
  }

  return pred_out;

}

// [[Rcpp::export]]
List predCopNNMP_covars(const arma::colvec& obs,
                        const int& cop_family,
                        const arma::cube& cop_param,
                        const int& marg_family,
                        const arma::mat ref_marg_param1,
                        const arma::mat ref_marg_param2,
                        const arma::mat nonref_marg_param1,
                        const arma::mat nonref_marg_param2,
                        const int& nne,
                        const arma::mat& DD,
                        const arma::colvec zeta,
                        const arma::mat& ga,
                        const arma::colvec kasq,
                        const arma::mat grid_ne_idx,
                        const arma::mat grid_ne_dist,
                        const arma::colvec& probs,
                        const bool& verbose,
                        const int& nreport,
                        const bool sam = false){

  List pred_out;
  
  if (marg_family == 1)
  {
    pred_out = predCop_covars(obs, cop_family, cop_param, 
                              ref_marg_param1, ref_marg_param2,
                              nonref_marg_param1, nonref_marg_param2, nne, DD,
                              zeta, ga, kasq, grid_ne_idx, grid_ne_dist, probs,
                              sam, pbeta, qbeta, verbose, nreport);
  } 
  else if (marg_family == 2)
  {
    pred_out = predCop_covars(obs, cop_family, cop_param, 
                              ref_marg_param1, ref_marg_param2,
                              nonref_marg_param1, nonref_marg_param2, nne, DD,
                              zeta, ga, kasq, grid_ne_idx, grid_ne_dist, probs,
                              sam, pgamma, qgamma, verbose, nreport);
  } 
  else {
    stop("error: this family is currently not supported.");
  }

  return pred_out;

}


