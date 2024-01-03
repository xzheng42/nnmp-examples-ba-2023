#ifndef COPULAS_H
#define COPULAS_H

double rConCop(const double& vv, 
               const int& cop_family, 
               const double& cop_param, 
               const double& tol = 1e-5, 
               const int& maxiter = 1000);

double  pConCop(const double& uu,
                const double& vv,
                const int& cop_family,
                const double& cop_param,
                const bool& logp = false);

#endif 