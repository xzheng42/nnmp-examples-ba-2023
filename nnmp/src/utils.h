#ifndef UTILS_H
#define UTILS_H
#include <RcppArmadillo.h>

arma::colvec armaIntSample(const arma::colvec& x, const int& size, const bool& replace, const arma::colvec& prob);

arma::colvec cumsum_cpp(const arma::colvec xx);

double qbeta(double uu, double shape1, double shape2);

double pbeta(double uu, double shape1, double shape2);

double qgamma(double uu, double shape1, double shape2);

double pgamma(double uu, double shape1, double shape2);
  
#endif 