#include <iostream>
#include <RcppArmadillo.h>
#include <math.h>
#include <truncnorm.h>
#include "utils.h"

using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

//////////////////////////////////////////////////////////////////////////
// Update skew-Gaussian NNMP's labels 
//////////////////////////////////////////////////////////////////////////

arma::colvec logC_cpp(const double& yy,
                      const arma::colvec& yy_ne,
                      const arma::colvec& rho,
                      const double& la,
                      const double& sigmasq) {
  
  int nne = yy_ne.n_rows;
  double s2, tt1, tt2;
  arma::colvec numer(nne), denom(nne), diff;
  
  tt2 = la / sqrt(sigmasq) / sqrt(pow(la, 2) + sigmasq);
  
  for (int j = 0; j < nne; ++j) {
    s2 = sigmasq * (1 + rho(j));
    tt1 = la / sqrt(s2) / sqrt(2.0 * pow(la, 2) + s2);
    numer(j) = R::pnorm(tt1 * (yy + yy_ne(j)), 0.0, 1.0, true, true);
    denom(j) = R::pnorm(tt2 * yy_ne(j), 0.0, 1.0, true, true);
  }
  
  diff = numer - denom;
  return diff;
  
}

arma::colvec logf_cpp(const double& yy,
                      const arma::colvec& yy_ne,
                      const arma::colvec& rho,
                      const double& la,
                      const double& sigmasq) {
  
  int nne = yy_ne.n_rows;
  double omega2 = sigmasq + pow(la, 2);
  double ss12, aa, tomega2;
  arma::colvec logdens(nne);
  
  for (int j = 0; j < nne; ++j) {
    ss12 = rho(j) * sigmasq + pow(la, 2);
    aa = ss12 / omega2;
    tomega2 = omega2 * (1.0 - pow(aa, 2));
    logdens(j) = R::dnorm(yy, aa * yy_ne(j), sqrt(tomega2), true);
  }
  
  return logdens;
  
}

// [[Rcpp::export]]
List updateSNnnmpLabel(const arma::colvec& yy, 
                       const arma::mat& yy_ne, 
                       const int& nne, 
                       const arma::mat& weight, 
                       const double& la,
                       const double& sigmasq,
                       const arma::mat& rho,
                       const arma::colvec& logit_mu,
                       const arma::colvec& logit_ka,
                       const arma::mat& cutoff) {
  
  int data_size = yy.n_rows;
  int label;
  arma::colvec label_vals = arma::linspace<arma::colvec>(1, nne, nne); 
  arma::Col<int> labels(data_size);
  arma::colvec tt(data_size), depNe(data_size);
  arma::colvec idistr, irho, iweight, iyy_ne, iloglik, ilogc, ilogf;
  double iyy, cutoff1, cutoff2, cutoffa, cutoffb;    
  
  for (int i = 0; i < data_size; ++i) {
    
    // Update configuration variables
    iyy = yy(i);
    iyy_ne = yy_ne.row(i).t();
    irho = rho.row(i).t();
    iweight = weight.row(i).t();
    ilogc = logC_cpp(iyy, iyy_ne, irho, la, sigmasq);
    ilogf = logf_cpp(iyy, iyy_ne, irho, la, sigmasq);
    iloglik = ilogc + ilogf;
    idistr = log(iweight) + iloglik;
    idistr = exp(idistr - idistr.max());
    idistr = idistr / sum(idistr);
    label = as_scalar(armaIntSample(label_vals, 1, 1, idistr));
    labels(i) = label;
    
    // Update sre variables from the truncated Gaussian 
    cutoff1 = cutoff(i, label - 1);
    cutoff2 = cutoff(i, label);
    cutoffa = log(cutoff1 / (1.0 - cutoff1));
    cutoffb = log(cutoff2 / (1.0 - cutoff2));
    tt(i) = r_truncnorm(logit_mu(i), logit_ka(i), cutoffa, cutoffb);
    
  }
  
  return List::create(Named("data_label") = labels, Named("latent_t") = tt);
  
}

arma::colvec logC_cpp_covars_lapar(const double& yy,
                                   const arma::colvec& yy_ne,
                                   const arma::colvec& rho,
                                   const double& nu,
                                   const arma::colvec& nu_ne,
                                   const double& la,
                                   const arma::colvec& la_ne,
                                   const double& sigma2) {
  
  int nne = yy_ne.n_rows;
  double s2, denom1, alp1, alp2, alp_ne, la_ne_j, la_ne_j2, omega2_ne;
  double la2 = pow(la, 2);
  arma::colvec numer(nne), denom(nne), diff;
  
  for (int j = 0; j < nne; ++j) {
    la_ne_j = la_ne(j);
    la_ne_j2 = pow(la_ne_j, 2);
    s2 = sigma2 * (1.0 + rho(j));
    denom1 = sqrt((1.0 - rho(j)) * s2) * sqrt((1.0 - rho(j)) * s2 + la2 + la_ne_j2 - 2.0 * rho(j) * la * la_ne_j);
    alp1 = (la - rho(j) * la_ne_j) / denom1;
    alp2 = (la_ne_j - rho(j) * la) / denom1;
    alp_ne = la_ne_j / sqrt(sigma2);
    omega2_ne = sigma2 + la_ne_j2;
    numer(j) = R::pnorm(alp1 * (yy - nu) + alp2 * (yy_ne(j) - nu_ne(j)), 0.0, 1.0, true, true);
    denom(j) = R::pnorm(alp_ne * (yy_ne(j) - nu_ne(j)) / sqrt(omega2_ne), 0.0, 1.0, true, true);
  }
  
  diff = numer - denom;
  return diff;
  
}

arma::colvec logf_cpp_covars_lapar(const double& yy,
                                   const arma::colvec& yy_ne,
                                   const arma::colvec& rho,
                                   const double& nu,
                                   const arma::colvec& nu_ne,
                                   const double& la,
                                   const arma::colvec& la_ne,
                                   const double& sigma2) {
  
  int nne = yy_ne.n_rows;
  double trho, aa, tomega2, omega2_ne, la_ne_j, la_ne_j2;
  double la2 = pow(la, 2);
  double omega2 = sigma2 + la2;
  arma::colvec logdens(nne);
  
  for (int j = 0; j < nne; ++j) {
    la_ne_j = la_ne(j);
    la_ne_j2 = pow(la_ne_j, 2);
    trho = rho(j) * sigma2 + la * la_ne_j;
    omega2_ne = sigma2 + la_ne_j2;
    aa = trho / omega2_ne;
    tomega2 = omega2 - pow(trho, 2) / omega2_ne;
    logdens(j) = R::dnorm(yy, nu +  aa * (yy_ne(j) - nu_ne(j)), sqrt(tomega2), true);
  }
  
  return logdens;
  
}

// [[Rcpp::export]]
List updateSNnnmpLabel_covars_lapar(const arma::colvec& yy, 
                                    const arma::mat& yy_ne, 
                                    const int& nne, 
                                    const arma::mat& weight, 
                                    const arma::colvec& nu,
                                    const arma::mat& nu_ne,
                                    const arma::colvec& la,
                                    const arma::mat& la_ne,
                                    const double& sigmasq,
                                    const arma::mat& rho,
                                    const arma::colvec& logit_mu,
                                    const arma::colvec& logit_ka,
                                    const arma::mat& cutoff) {
  
  int data_size = yy.n_rows;
  int label;
  arma::colvec label_vals = arma::linspace<arma::colvec>(1, nne, nne); 
  arma::Col<int> labels(data_size);
  arma::colvec tt(data_size), depNe(data_size);
  arma::colvec idistr, irho, iweight, iyy_ne, iloglik, ilogc, ilogf, inu_ne, ila_ne;
  double iyy, cutoff1, cutoff2, cutoffa, cutoffb, inu, ila;    
  
  for (int i = 0; i < data_size; ++i) {
    
    // Update configuration variables
    iyy = yy(i);
    inu = nu(i);
    ila = la(i);
    iyy_ne = yy_ne.row(i).t();
    inu_ne = nu_ne.row(i).t();
    ila_ne = la_ne.row(i).t();
    irho = rho.row(i).t();
    iweight = weight.row(i).t();
    ilogc = logC_cpp_covars_lapar(iyy, iyy_ne, irho, inu, inu_ne, ila, ila_ne, sigmasq);
    ilogf = logf_cpp_covars_lapar(iyy, iyy_ne, irho, inu, inu_ne, ila, ila_ne, sigmasq);
    iloglik = ilogc + ilogf;
    idistr = log(iweight) + iloglik;
    idistr = exp(idistr - idistr.max());
    idistr = idistr / sum(idistr);
    label = as_scalar(armaIntSample(label_vals, 1, 1, idistr));
    labels(i) = label;
    
    // Update sre variables from the truncated Gaussian 
    cutoff1 = cutoff(i, label - 1);
    cutoff2 = cutoff(i, label);
    cutoffa = log(cutoff1 / (1.0 - cutoff1));
    cutoffb = log(cutoff2 / (1.0 - cutoff2));
    tt(i) = r_truncnorm(logit_mu(i), logit_ka(i), cutoffa, cutoffb);
    
  }
  
  return List::create(Named("data_label") = labels, Named("latent_t") = tt);
  
}