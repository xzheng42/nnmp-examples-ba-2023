#include <iostream> 
#include <RcppArmadillo.h> 
#include <RcppDist.h> 
#include <RcppTN.h> 
#include <truncnorm.h> 
#include <math.h> 
#include "utils.h"  

using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo, RcppDist, RcppTN)]]

#define PMAX  1.0 - 1e-10
#define PMIN  1e-10

double dBiCop_cpp(const double& uu,
                  const double& vv,
                  const int& cop_family,
                  const double& cop_param,
                  const bool& logd = false) {
  
  double logdens;
  
  
  // Gaussian
  if (cop_family == 1) {
    double xx = R::qnorm(uu, 0.0, 1.0, true, false);
    double yy = R::qnorm(vv, 0.0, 1.0, true, false);
    double rho = cop_param;
    double rho2 = pow(rho, 2);
    double xx2 = pow(xx, 2);
    double yy2 = pow(yy, 2);
    double ee = rho2 * (xx2 + yy2) - 2.0 * rho * xx * yy;
    logdens = -.5 * ee / (1.0 - rho2) - .5 * log(1.0 - rho2);   
  } 
  // Gumbel
  else if (cop_family == 2) {
    double eta = cop_param;
    if (eta < 1) {
      stop("error: gumbel copula parameter must be greater than 1.");
    }
    if (eta > 50.0) {
      eta = 50.0;
    }
    double xx = -log(uu);
    double yy = -log(vv);
    double tt = pow(xx, eta) + pow(yy, eta);
    logdens = -pow(tt, 1.0 / eta) + log(pow(tt, 1.0 / eta) + eta - 1.0) + 
      (1.0 / eta - 2.0) * log(tt) + (eta - 1.0) * (log(xx) + log(yy)) + xx + yy;
  }
  // Clayton
  else if (cop_family == 3) {
    double delta = cop_param;
    if (delta > 98) {
      delta = 98;
    }
    double ee = pow(uu, -delta) + pow(vv, -delta) - 1.0;
    logdens = log(1.0 + delta) - (delta + 1.0) * (log(uu) + log(vv)) - (2.0 + 1.0 / delta) * log(ee);
  }
  else {
    stop("error: no families are found.");
  }
  
  if (logd) {
    return logdens;
  } else {
    return exp(logdens);
  }
  
}


arma::colvec dBiCop_cpp2(const arma::colvec& uu, 
                         const arma::colvec& vv,
                         const int& cop_family,
                         const arma::colvec& cop_param,
                         const bool& logd = false) {
  
  int nn = uu.n_rows;
  arma::colvec dens(nn);
  
  for (int i = 0; i < nn; ++i) {
    dens(i) = dBiCop_cpp(uu(i), vv(i), cop_family, cop_param(i), logd);
  } 
  
  return dens;
  
}

// [[Rcpp::export]]
double dBiCopMar_cpp(const double& xx,
                     const double& yy,
                     const int& cop_family,
                     const double& cop_param,
                     const int& marg_family,
                     const arma::colvec& mar_param,
                     const bool& ce = false,
                     const bool& logd = false) {
  
  double uu, vv;
  
  if (marg_family == 1) {
    uu = R::pbeta(xx, mar_param(0), mar_param(1), true, false);
    vv = R::pbeta(yy, mar_param(2), mar_param(3), true, false);
  }
  else if (marg_family == 2) {
    uu = R::pgamma(xx, mar_param(0), 1.0 / mar_param(1), true, false);
    vv = R::pgamma(yy, mar_param(2), 1.0 / mar_param(3), true, false);
  } 
  else {
    stop("error: no such a family is available.");
  }
  
  if (uu < PMIN) {
    uu = PMIN;
  } else if (uu > PMAX) {
    uu = PMAX;
  }
  if (vv < PMIN) {
    vv = PMIN;
  } else if (vv > PMAX) {
    vv = PMAX;
  }
  
  return dBiCop_cpp(uu, vv, cop_family, cop_param, logd);
  
  
}

// [[Rcpp::export]]
arma::colvec dBiCopMar_cpp2(const arma::colvec& xx,
                            const arma::colvec& yy,
                            const int& cop_family,
                            const arma::colvec cop_param,
                            const int& marg_family,
                            const arma::mat& mar_param,
                            const bool& ce = false,
                            const bool& logd = false) {
  
  int nn = xx.n_rows;    
  arma::colvec uu(nn), vv(nn);
  
  if (marg_family == 1) {
    for (int i = 0; i < nn; ++i) {
      uu(i) = R::pbeta(xx(i), mar_param(i,0), mar_param(i,1), true, false);
      vv(i) = R::pbeta(yy(i), mar_param(i,2), mar_param(i,3), true, false);
    }
  }
  else if (marg_family == 2) {
    for (int i = 0; i < nn; ++i) {
      uu(i) = R::pgamma(xx(i), mar_param(i,0), 1.0 / mar_param(i,1), true, false);
      vv(i) = R::pgamma(yy(i), mar_param(i,2), 1.0 / mar_param(i,3), true, false);
    }
  } 
  else {
    stop("error: no such a family is available.");
  }
  
  for (int i = 0; i < nn; ++i) {
    if (uu(i) < PMIN) {
      uu(i) = PMIN;
    } else if (uu(i) > PMAX) {
      uu(i) = PMAX;
    }
    if (vv(i) < PMIN) {
      vv(i) = PMIN;
    } else if (vv(i) > PMAX) {
      vv(i) = PMAX;
    }
  }
  
  return dBiCop_cpp2(uu, vv, cop_family, cop_param, logd);    
  
}


// [[Rcpp::export]]
List updateCopLabel(const arma::colvec& yy,
                    const arma::mat& yy_ne,
                    const int& nne,
                    const arma::mat& weight,
                    const arma::mat& rho,
                    const int& cop_family,
                    const int& marg_family,
                    const arma::colvec& yy_param1,
                    const arma::colvec& yy_param2,
                    const arma::mat& yy_ne_param1,
                    const arma::mat& yy_ne_param2,
                    const arma::colvec& logit_mu,
                    const arma::colvec& logit_ka,
                    const arma::mat& cutoff) {
  
  int data_size = yy.n_rows;
  int ilabel;
  arma::Col<int> labels(data_size);
  arma::colvec label_vals = arma::linspace<arma::colvec>(1, nne, nne);
  arma::colvec tt(data_size);
  arma::colvec ilabel_val, idistr, irho, iweight, iloglik, iyy_ne, iyy;
  arma::colvec iyy_param1, iyy_param2;
  arma::mat mar_param;
  double cutoff1, cutoff2, cutoffa, cutoffb;
  
  tt.fill(NA_REAL);
  
  iyy.set_size(nne);
  iyy_param1.set_size(nne);
  iyy_param2.set_size(nne);
  
  for (int i = 0; i < data_size; ++i) {
    
    iyy.fill(yy(i));
    iyy_ne = yy_ne.row(i).t();
    irho = rho.row(i).t();
    
    mar_param.set_size(nne, 4);
    iyy_param1.fill(yy_param1(i));
    iyy_param2.fill(yy_param2(i));
    mar_param.col(0) = iyy_param1;
    mar_param.col(1) = iyy_param2;
    mar_param.col(2) = yy_ne_param1.row(i).t();
    mar_param.col(3) = yy_ne_param2.row(i).t();
    
    iloglik = dBiCopMar_cpp2(iyy, iyy_ne, cop_family, irho, marg_family, mar_param, false, true);
    iweight = weight.row(i).t();
    idistr = log(iweight) + iloglik;
    idistr = exp(idistr - idistr.max());
    idistr = idistr / sum(idistr);
    ilabel = as_scalar(armaIntSample(label_vals, 1, 1, idistr));
    labels(i) = ilabel;
    
    cutoff1 = cutoff(i, ilabel - 1);
    cutoff2 = cutoff(i, ilabel);
    cutoffa = log(cutoff1 / (1.0 - cutoff1));
    cutoffb = log(cutoff2 / (1.0 - cutoff2));
    tt(i) = r_truncnorm(logit_mu(i), logit_ka(i), cutoffa, cutoffb);
    
  }
  
  return List::create(Named("data_label") = labels, Named("latent_t") = tt);
  
}


// [[Rcpp::export]]
List updateCeCopAux(const arma::colvec& yy,
                    const arma::colvec& oo,
                    const arma::colvec& rho,
                    const arma::colvec& ne_label,
                    const int& cop_family,
                    const int& mar_family,
                    const arma::colvec& yy_param1,
                    const arma::colvec& yy_param2) {
  
  int nn = yy.n_rows;
  int kk, ilabel;
  arma::colvec new_oo = oo;
  arma::uvec as_ne_idx;
  arma::colvec iyy_pa, ioo_pa, irho_pa, iprop_yy, icur_yy; 
  arma::colvec iparam1, iparam2, mar_param1;
  arma::mat mar_param, mar_param2;
  double prop_oo, cur_oo, prop_loglik, cur_loglik, diff_loglik, irho, iyy_ne, ioo_ne;
  arma::colvec oo_acct_save(nn);
  oo_acct_save.zeros();
  
  // When i = 1
  as_ne_idx = arma::find(ne_label == 1);
  kk = as_ne_idx.n_rows;
  irho_pa = rho(as_ne_idx);
  iyy_pa = yy(as_ne_idx) - new_oo(as_ne_idx);
  
  prop_oo = R::runif(0, 1);
  cur_oo = new_oo(0);  
  iprop_yy.set_size(kk);
  icur_yy.set_size(kk);
  iprop_yy.fill(yy(0) - prop_oo);
  icur_yy.fill(yy(0) - cur_oo);
  
  if (mar_family != 3) {
    iparam1.set_size(kk);
    iparam2.set_size(kk);
    iparam1.fill(yy_param1(0));
    iparam2.fill(yy_param2(0));
    mar_param.set_size(kk, 4);
    mar_param.col(0) = yy_param1(as_ne_idx);
    mar_param.col(1) = yy_param2(as_ne_idx);
    mar_param.col(2) = iparam1;
    mar_param.col(3) = iparam2;
  } else {
    iparam1.set_size(kk);
    iparam1.fill(yy_param1(0));
    mar_param.set_size(kk, 2);
    mar_param.col(0) = yy_param1(as_ne_idx);
    mar_param.col(1) = iparam1;
  }
  
  
  prop_loglik = sum(dBiCopMar_cpp2(iyy_pa, iprop_yy, cop_family, irho_pa, mar_family, mar_param, true, true));
  cur_loglik = sum(dBiCopMar_cpp2(iyy_pa, icur_yy, cop_family, irho_pa, mar_family, mar_param, true, true));
  diff_loglik = prop_loglik - cur_loglik;
  if (diff_loglik > log(R::runif(0, 1))) {
    new_oo(0) = prop_oo;
    oo_acct_save(0) = 1;
  }
  
  
  // when i > 1
  for (int i = 1; i < nn; ++i) {
    
    ilabel = ne_label(i) - 1;
    iyy_ne = yy(ilabel);
    ioo_ne = new_oo(ilabel);
    irho = rho(i);
    
    as_ne_idx = arma::find(ne_label == i + 1);
    kk = as_ne_idx.n_rows;
    iyy_pa = yy(as_ne_idx) - new_oo(as_ne_idx);
    irho_pa = rho(as_ne_idx);
    
    prop_oo = R::runif(0, 1);
    cur_oo = new_oo(i);
    iprop_yy.set_size(kk);
    icur_yy.set_size(kk);
    iprop_yy.fill(yy(i) - prop_oo);
    icur_yy.fill(yy(i) - cur_oo);
    
    if (mar_family != 3) {
      iparam1.set_size(kk);
      iparam2.set_size(kk);
      iparam1.fill(yy_param1(i));
      iparam2.fill(yy_param2(i));
      mar_param1 = { yy_param1(i), yy_param2(i), yy_param1(ilabel), yy_param2(ilabel) };
      mar_param2.set_size(kk, 4);
      mar_param2.col(0) = yy_param1(as_ne_idx);
      mar_param2.col(1) = yy_param2(as_ne_idx);
      mar_param2.col(2) = iparam1;
      mar_param2.col(3) = iparam2;
    } else {
      mar_param1 = { yy_param1(i), yy_param1(ilabel) };
      iparam1.set_size(kk);
      iparam1.fill(yy_param1(0));
      mar_param2.set_size(kk, 2);
      mar_param2.col(0) = yy_param1(as_ne_idx);
      mar_param2.col(1) = iparam1;    
    }
    
    prop_loglik = dBiCopMar_cpp(yy(i) - prop_oo, iyy_ne - ioo_ne, cop_family, irho, mar_family, mar_param1, true, true) +
      sum(dBiCopMar_cpp2(iyy_pa, iprop_yy, cop_family, irho_pa, mar_family, mar_param2, true, true));
    cur_loglik = dBiCopMar_cpp(yy(i) - cur_oo, iyy_ne - ioo_ne, cop_family, irho, mar_family, mar_param1, true, true) +
      sum(dBiCopMar_cpp2(iyy_pa, icur_yy, cop_family, irho_pa, mar_family, mar_param2, true, true));
    diff_loglik = prop_loglik - cur_loglik;
    if (diff_loglik > log(R::runif(0, 1))) {
      new_oo(i) = prop_oo;
      oo_acct_save(i) = 1;
    }
    
  }
  
  return List::create(Named("oo") = new_oo, Named("oo_acct_save") = oo_acct_save);
  
}


double rConGum_mnr(const double& vv,
                   const double& eta,
                   const double& tol,
                   const int& maxiter) {
  
  double ww = -log(vv);
  double pp = arma::randu();
  double upb = 100.0;
  double lob = ww;
  double zz0 = (lob + upb) / 2.0;
  
  double rr, tt, tt1, tt2, zz1, uu;
  
  for (int i = 0; i < maxiter; ++i) {    
    
    // Newton step
    tt = zz0 + (eta - 1.0) * log(zz0) - (ww + (eta - 1.0) * log(ww) - log(pp));
    tt1 = 1.0 + (eta - 1.0) / zz0;
    tt2 = (1.0 - eta) / pow(zz0, 2);
    zz1 = zz0 - tt * tt1 / (tt1 * tt1  - tt * tt2);
    
    // Bisection step
    if (zz1 < lob || zz1 > upb) {
      if (tt < 0) {
        lob = zz0;
      } else {
        upb = zz0;
      }
      zz1 = (lob + upb) / 2.0;
      
    } 
    
    rr = zz1 + (eta - 1.0) * log(zz1) - (ww + (eta - 1.0) * log(ww) - log(pp));
    
    if (abs(rr) <= tol) {
      break;
    }
    
    zz0 = zz1;
    
    if (isnan(rr)) {
      stop("Numerical issue occurred: the function outputs nan");
    }
    if (i == (maxiter - 1)){
      stop("Failed to converge!");
    }
  }
  
  uu = exp(-pow((pow(zz1, eta) - pow(ww, eta)), 1.0/eta));
  
  if (uu == 1){
    // cout << "eta = " << eta << endl;
    stop("Numerical issue occurred: uu = 1");
  }
  
  return uu;
  
}


double rConCop(const double& vv,
               const int& cop_family,
               const double& cop_param,
               const double& tol = 1e-5,
               const int& maxiter = 1000) {
  
  double zz = arma::randu();
  double uu;
  
  if (cop_family == 1) {
    double rho = cop_param;
    double qzz = R::qnorm(zz, 0.0, 1.0, true, false); 
    double qvv = R::qnorm(vv, 0.0, 1.0, true, false);
    uu = R::pnorm(sqrt(1.0 - rho * rho) * qzz + rho * qvv, 0.0, 1.0, true, false);
  }
  else if (cop_family == 2) {
    double eta = cop_param;
    if (eta > 50.0) {
      eta = 50.0;
    }
    uu = rConGum_mnr(vv, eta, tol, maxiter);
  }
  else if (cop_family == 3) {
    double delta = cop_param;
    if (delta > 98.0) {
      delta = 98.0;
    }
    double tt = exp(-delta / (1.0 + delta) * log(zz)) - 1.0;
    uu = pow(exp(log(tt) - delta * log(vv)) + 1.0, -1.0 / delta);
  }
  else {
    stop("error: no families are found.");
  }
  
  return uu;
  
}


double pConCop(const double& uu,
               const double& vv,
               const int& cop_family,
               const double& cop_param,
               const bool& logp = false) {
  
  if (cop_family == 1) {
    double rho = cop_param;
    double xx = R::qnorm(uu, 0.0, 1.0, true, false);
    double yy = R::qnorm(vv, 0.0, 1.0, true, false);
    double qq = (xx - rho * yy) / sqrt(1.0 - pow(rho, 2.0));
    double pp = R::pnorm(qq, 0.0, 1.0, true, false);
    if (logp) {
      return log(pp);
    } else {
      return pp;
    }
    
  } 
  else if (cop_family == 2) {
    
    double eta = cop_param;
    if (eta > 50.0) {
      eta = 50.0;
    }
    
    double xx = -log(uu);
    double yy = -log(vv);
    double logpp = -log(vv) - pow(pow(xx, eta) + pow(yy, eta), 1.0 / eta) + 
      (1.0 / eta - 1.0) * log(1.0 + pow(xx / yy, eta));
    
    if (logp) {
      return(logpp);
    } else {
      return exp(logpp);
    }
    
  } 
  else if (cop_family == 3) {
    
    double delta = cop_param;
    if (delta > 98.0) {
      delta = 98.0;
    }
    
    double ee = 1.0 + pow(vv, delta) * (pow(uu, -delta)  - 1);
    double logpp = -(1.0 + 1.0 / delta) * log(ee);
    
    if (logp) {
      return logpp;
    } else {
      return exp(logpp);
    }
    
  } 
  else {
    
    stop("error: no such a family is supported.");
    
  }
  
}