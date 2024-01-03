#include <iostream>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>

using namespace std;
using namespace Rcpp;

//////////////////////////////////////////////////////////////////////////
// RcppArmadillo's sample function
/////////////////////////////////////////////////////////////////////////
arma::colvec armaIntSample(const arma::colvec& x,
                           const int& size,
                           const bool& replace,
                           const arma::colvec& prob) {
  arma::colvec rsq = RcppArmadillo::sample(x, size, replace, prob);
  return(rsq);
}


//////////////////////////////////////////////////////////////////////////
// Compute cumulative sum
//////////////////////////////////////////////////////////////////////////
arma::colvec cumsum_cpp(const arma::colvec xx) {
  int kk = xx.n_rows;
  arma::colvec yy(kk);
  yy(0) = xx(0);
  for (int i = 1; i < kk; ++i) {
    yy(i) = sum(xx.rows(0, i));
    if (yy(i) > 1.0) {
      yy.rows(i, kk - 1).ones();
      break;
    }
  }
  return yy;
}


//////////////////////////////////////////////////////////////////////////
// Distribution functions
//////////////////////////////////////////////////////////////////////////
#define EPSILON DBL_EPSILON

double qbeta(double uu, double shape1, double shape2) {
  return R::qbeta(uu, shape1, shape2, true, false);
}

double pbeta(double uu, double shape1, double shape2) {
  return R::pbeta(uu, shape1, shape2, true, false);
}

double qgamma(double uu, double shape1, double shape2) {
  return R::qgamma(uu, shape1, 1.0 / shape2, true, false);
}

double pgamma(double uu, double shape1, double shape2) {
  return R::pgamma(uu, shape1, 1.0 / shape2, true, false);
}


