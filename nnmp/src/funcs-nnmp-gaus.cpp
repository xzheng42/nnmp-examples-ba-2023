#include <iostream>
#include <RcppArmadillo.h>
#include <math.h>
#include <truncnorm.h>
#include "utils.h"

using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

//////////////////////////////////////////////////////////////////////////
// Update the spatial random effects of the Guassian NNMP regression model
//////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::colvec updateGnnmpSPE(const arma::colvec& yy,
                            const arma::mat& XX,
                            const arma::colvec& bb,
                            const arma::colvec& sre,
                            const arma::colvec& sre_ne_label,
                            const arma::colvec& rho,
                            const double& tausq,
                            const double& sigmasq) {
  
  int nn = yy.n_rows, label;
  double ine, isigma2, imu;
  arma::uvec as_ne_idx;
  arma::colvec rho2 = arma::pow(rho, 2);
  arma::colvec new_sre = sre;
  arma::colvec resid = yy - XX * bb;
  arma::colvec rho_sigma2 = sigmasq * (1.0 - rho2);
  arma::colvec ipa, ipa_rho, is2;
  
  // When i = 1
  ine = 0;
  as_ne_idx = arma::find(sre_ne_label == 1);
  ipa = new_sre(as_ne_idx);
  ipa_rho = rho(as_ne_idx);
  is2 = rho_sigma2(as_ne_idx) / rho2(as_ne_idx);
  isigma2 = 1.0 / (1.0 / tausq + 1.0 / sigmasq + sum(1.0 / is2));
  imu = resid(0) / tausq + sum(ipa / (ipa_rho % is2));
  new_sre(0) = R::rnorm(imu * isigma2, sqrt(isigma2));
  
  // When i > 1
  for (int i = 1; i < nn; ++i) {
    
    label = sre_ne_label(i) - 1;
    ine = new_sre(label);        
    as_ne_idx = arma::find(sre_ne_label == i + 1); // find indices of z_j's whose neighbor is z_i
    ipa = new_sre(as_ne_idx); // find those z_j whose neighbor is z_i
    ipa_rho = rho(as_ne_idx); // find the corresponding rho_j
    is2 = rho_sigma2(as_ne_idx) / rho2(as_ne_idx); 
    isigma2 = 1.0 / (1.0 / tausq + 1.0 / rho_sigma2(i) + sum(1.0 / is2));
    imu = resid(i) / tausq + rho(i) * ine / rho_sigma2(i) + sum(ipa / (ipa_rho % is2));
    new_sre(i) = R::rnorm(imu * isigma2, sqrt(isigma2));
    
  }
  
  return(new_sre);
  
}


//////////////////////////////////////////////////////////////////////////
// Update the configuration variables of the Guassian NNMP regression model
//////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List updateGnnmpLabel(const arma::colvec& sre, 
                      const arma::mat& sre_ne, 
                      const int& nne, 
                      const arma::mat& weight, 
                      const arma::colvec& mu, 
                      const arma::colvec& sigmasq,
                      const arma::mat& rho,
                      const arma::colvec& logit_mu,
                      const arma::colvec& logit_ka,
                      const arma::mat& cutoff) {
  
  int data_size = sre.n_rows, label;
  arma::colvec label_vals = arma::linspace<arma::colvec>(1, nne, nne); 
  arma::Col<int> labels(data_size);
  arma::colvec tt(data_size);
  arma::colvec idistr, isre, irho, iweight, isre_ne, imu, isigma2, iloglik, ilabel_vals;  
  arma::colvec cmu, csigma2;
  double cutoff1, cutoff2, cutoffa, cutoffb;
  
  labels(0) = NA_REAL, labels(1) = 1;
  tt.fill(NA_REAL);
  
  for (int i = 2; i < nne; ++i) {
    // Update configuration variables
    isre.set_size(i);
    isre.fill(sre(i));
    cmu = mu.rows(0, i - 1);
    csigma2.set_size(i);
    csigma2 = sigmasq.rows(0, i - 1);
    ilabel_vals = arma::linspace<arma::colvec>(1, i, i);
    irho = arma::trans(rho(i, arma::span(0, i - 1)));
    iweight = arma::trans(weight(i, arma::span(0, i - 1)));
    isre_ne = arma::trans(sre_ne(i, arma::span(0, i - 1)));
    imu = cmu + irho % (isre_ne - cmu);
    isigma2 = csigma2 % (1.0 - arma::pow(irho, 2));
    iloglik = arma::log_normpdf(isre, imu, arma::sqrt(isigma2));
    idistr = log(iweight) + iloglik;
    idistr = exp(idistr - idistr.max());
    idistr = idistr / sum(idistr);
    label = arma::as_scalar(armaIntSample(ilabel_vals, 1, 1, idistr));
    labels(i) = label;
    // Update sre variables from the truncated Gaussian 
    cutoff1 = cutoff(i, label - 1);
    cutoff2 = cutoff(i, label);
    cutoffa = log(cutoff1 / (1.0 - cutoff1));
    cutoffb = log(cutoff2 / (1.0 - cutoff2));
    tt(i) = r_truncnorm(logit_mu(i-2), logit_ka(i-2), cutoffa, cutoffb);
  }
  
  isre.set_size(nne);
  for (int i = nne; i < data_size; ++i) {
    // Update configuration variables
    irho = rho.row(i).t();
    iweight = weight.row(i).t();
    isre.fill(sre(i));
    isre_ne = sre_ne.row(i).t();
    imu = mu + irho % (isre_ne - mu);
    isigma2 = sigmasq % (1.0 - arma::pow(irho, 2));
    iloglik = arma::log_normpdf(isre, imu, arma::sqrt(isigma2));
    idistr = log(iweight) + iloglik;
    idistr = exp(idistr - idistr.max());
    idistr = idistr / sum(idistr);
    label = arma::as_scalar(armaIntSample(label_vals, 1, 1, idistr));
    labels(i) = label;
    // Update sre variables from the truncated Gaussian 
    cutoff1 = cutoff(i, label - 1);
    cutoff2 = cutoff(i, label);
    cutoffa = log(cutoff1 / (1.0 - cutoff1));
    cutoffb = log(cutoff2 / (1.0 - cutoff2));
    tt(i) = r_truncnorm(logit_mu(i-2), logit_ka(i-2), cutoffa, cutoffb);
  }
  
  return List::create(Named("latent_label") = labels, Named("latent_t") = tt);
  
}


//////////////////////////////////////////////////////////////////////////
// Compute Gelfand and Ghosh PPLC for the Gaussian NNMP spatial regression
/////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::colvec pplcGausNNMP(const arma::colvec& yy,
                           const arma::mat& XX,
                           const arma::mat& bb,
                           const arma::mat& zz,
                           const arma::colvec& tausq, 
                           const int& weight = 1.0) {
  
  int data_size = yy.n_rows;
  int sample_size = tausq.n_rows;
  double ijmu, idiff;
  arma::colvec rep_data(sample_size);
  arma::colvec gscore(data_size);
  arma::colvec pscore(data_size);
  arma::rowvec iXX;
  arma::colvec gg(3);
  
  for (int i = 0; i < data_size; ++i) {
    iXX = XX.row(i);
    for (int j = 0; j < sample_size; ++j) {
      ijmu = arma::as_scalar(iXX * bb.col(j) + zz(i,j));
      rep_data(j) = R::rnorm(ijmu, sqrt(tausq(j)));
    }
    idiff = yy(i) - mean(rep_data);
    gscore(i) = pow(idiff, 2);
    pscore(i) = var(rep_data);
  }
  
  gg(0) = sum(gscore);
  gg(1) = sum(pscore);
  gg(2) = weight * gg(0) + gg(1);
  
  return gg;
  
}


//////////////////////////////////////////////////////////////////////////
// Compute DIC for the Gaussian NNMP spatial regression
/////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double dicGausNNMP(const arma::colvec& yy,
                      const arma::mat& XX,
                      const arma::mat& bb,
                      const arma::mat& zz,
                      const arma::colvec& tausq) {
  
  int data_size = yy.n_rows;
  int sample_size = tausq.n_rows;
  double imu;
  double ijmu;
  arma::colvec bb_mean = arma::mean(bb, 1);
  arma::colvec zz_mean = arma::mean(zz, 1);
  double tausq_mean = mean(tausq);
  double dic;
  
  arma::colvec loglik(data_size);
  arma::mat loglik2(data_size, sample_size);
  arma::rowvec iXX;
  
  for (int i = 0; i < data_size; ++i) {
    iXX = XX.row(i);
    imu = arma::as_scalar(iXX * bb_mean + zz_mean(i));
    loglik(i) = R::dnorm(yy(i), imu, sqrt(tausq_mean), true);
    for (int j = 0; j < sample_size; ++j) {
      ijmu = arma::as_scalar(iXX * bb.col(j) + zz(i,j));
      loglik2(i,j) = R::dnorm(yy(i), ijmu, sqrt(tausq(j)), true);
    }
  }
  
  dic = arma::as_scalar(2.0 * sum(loglik) - 4.0 * sum(sum(loglik2)) / sample_size);
  
  return dic;
  
}