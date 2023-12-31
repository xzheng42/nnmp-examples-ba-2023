// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// dBiCopMar_cpp
double dBiCopMar_cpp(const double& xx, const double& yy, const int& cop_family, const double& cop_param, const int& marg_family, const arma::colvec& mar_param, const bool& ce, const bool& logd);
RcppExport SEXP _nnmp_dBiCopMar_cpp(SEXP xxSEXP, SEXP yySEXP, SEXP cop_familySEXP, SEXP cop_paramSEXP, SEXP marg_familySEXP, SEXP mar_paramSEXP, SEXP ceSEXP, SEXP logdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< const double& >::type yy(yySEXP);
    Rcpp::traits::input_parameter< const int& >::type cop_family(cop_familySEXP);
    Rcpp::traits::input_parameter< const double& >::type cop_param(cop_paramSEXP);
    Rcpp::traits::input_parameter< const int& >::type marg_family(marg_familySEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type mar_param(mar_paramSEXP);
    Rcpp::traits::input_parameter< const bool& >::type ce(ceSEXP);
    Rcpp::traits::input_parameter< const bool& >::type logd(logdSEXP);
    rcpp_result_gen = Rcpp::wrap(dBiCopMar_cpp(xx, yy, cop_family, cop_param, marg_family, mar_param, ce, logd));
    return rcpp_result_gen;
END_RCPP
}
// dBiCopMar_cpp2
arma::colvec dBiCopMar_cpp2(const arma::colvec& xx, const arma::colvec& yy, const int& cop_family, const arma::colvec cop_param, const int& marg_family, const arma::mat& mar_param, const bool& ce, const bool& logd);
RcppExport SEXP _nnmp_dBiCopMar_cpp2(SEXP xxSEXP, SEXP yySEXP, SEXP cop_familySEXP, SEXP cop_paramSEXP, SEXP marg_familySEXP, SEXP mar_paramSEXP, SEXP ceSEXP, SEXP logdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type yy(yySEXP);
    Rcpp::traits::input_parameter< const int& >::type cop_family(cop_familySEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type cop_param(cop_paramSEXP);
    Rcpp::traits::input_parameter< const int& >::type marg_family(marg_familySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mar_param(mar_paramSEXP);
    Rcpp::traits::input_parameter< const bool& >::type ce(ceSEXP);
    Rcpp::traits::input_parameter< const bool& >::type logd(logdSEXP);
    rcpp_result_gen = Rcpp::wrap(dBiCopMar_cpp2(xx, yy, cop_family, cop_param, marg_family, mar_param, ce, logd));
    return rcpp_result_gen;
END_RCPP
}
// updateCopLabel
List updateCopLabel(const arma::colvec& yy, const arma::mat& yy_ne, const int& nne, const arma::mat& weight, const arma::mat& rho, const int& cop_family, const int& marg_family, const arma::colvec& yy_param1, const arma::colvec& yy_param2, const arma::mat& yy_ne_param1, const arma::mat& yy_ne_param2, const arma::colvec& logit_mu, const arma::colvec& logit_ka, const arma::mat& cutoff);
RcppExport SEXP _nnmp_updateCopLabel(SEXP yySEXP, SEXP yy_neSEXP, SEXP nneSEXP, SEXP weightSEXP, SEXP rhoSEXP, SEXP cop_familySEXP, SEXP marg_familySEXP, SEXP yy_param1SEXP, SEXP yy_param2SEXP, SEXP yy_ne_param1SEXP, SEXP yy_ne_param2SEXP, SEXP logit_muSEXP, SEXP logit_kaSEXP, SEXP cutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type yy(yySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type yy_ne(yy_neSEXP);
    Rcpp::traits::input_parameter< const int& >::type nne(nneSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int& >::type cop_family(cop_familySEXP);
    Rcpp::traits::input_parameter< const int& >::type marg_family(marg_familySEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type yy_param1(yy_param1SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type yy_param2(yy_param2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type yy_ne_param1(yy_ne_param1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type yy_ne_param2(yy_ne_param2SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type logit_mu(logit_muSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type logit_ka(logit_kaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type cutoff(cutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(updateCopLabel(yy, yy_ne, nne, weight, rho, cop_family, marg_family, yy_param1, yy_param2, yy_ne_param1, yy_ne_param2, logit_mu, logit_ka, cutoff));
    return rcpp_result_gen;
END_RCPP
}
// updateCeCopAux
List updateCeCopAux(const arma::colvec& yy, const arma::colvec& oo, const arma::colvec& rho, const arma::colvec& ne_label, const int& cop_family, const int& mar_family, const arma::colvec& yy_param1, const arma::colvec& yy_param2);
RcppExport SEXP _nnmp_updateCeCopAux(SEXP yySEXP, SEXP ooSEXP, SEXP rhoSEXP, SEXP ne_labelSEXP, SEXP cop_familySEXP, SEXP mar_familySEXP, SEXP yy_param1SEXP, SEXP yy_param2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type yy(yySEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type oo(ooSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type ne_label(ne_labelSEXP);
    Rcpp::traits::input_parameter< const int& >::type cop_family(cop_familySEXP);
    Rcpp::traits::input_parameter< const int& >::type mar_family(mar_familySEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type yy_param1(yy_param1SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type yy_param2(yy_param2SEXP);
    rcpp_result_gen = Rcpp::wrap(updateCeCopAux(yy, oo, rho, ne_label, cop_family, mar_family, yy_param1, yy_param2));
    return rcpp_result_gen;
END_RCPP
}
// cdfCopNNMP_simple
arma::cube cdfCopNNMP_simple(const arma::mat& grid_ne_obs, const arma::mat& grid_ne_dist, const int& cop_family, const arma::cube& cop_param, const int& marg_family, const arma::mat& marg_param, const arma::colvec& zeta, const arma::mat& ga, const arma::colvec kasq, const arma::mat& DD, const arma::colvec& probs);
RcppExport SEXP _nnmp_cdfCopNNMP_simple(SEXP grid_ne_obsSEXP, SEXP grid_ne_distSEXP, SEXP cop_familySEXP, SEXP cop_paramSEXP, SEXP marg_familySEXP, SEXP marg_paramSEXP, SEXP zetaSEXP, SEXP gaSEXP, SEXP kasqSEXP, SEXP DDSEXP, SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type grid_ne_obs(grid_ne_obsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type grid_ne_dist(grid_ne_distSEXP);
    Rcpp::traits::input_parameter< const int& >::type cop_family(cop_familySEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type cop_param(cop_paramSEXP);
    Rcpp::traits::input_parameter< const int& >::type marg_family(marg_familySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type marg_param(marg_paramSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ga(gaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type kasq(kasqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type DD(DDSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(cdfCopNNMP_simple(grid_ne_obs, grid_ne_dist, cop_family, cop_param, marg_family, marg_param, zeta, ga, kasq, DD, probs));
    return rcpp_result_gen;
END_RCPP
}
// updateGnnmpSPE
arma::colvec updateGnnmpSPE(const arma::colvec& yy, const arma::mat& XX, const arma::colvec& bb, const arma::colvec& sre, const arma::colvec& sre_ne_label, const arma::colvec& rho, const double& tausq, const double& sigmasq);
RcppExport SEXP _nnmp_updateGnnmpSPE(SEXP yySEXP, SEXP XXSEXP, SEXP bbSEXP, SEXP sreSEXP, SEXP sre_ne_labelSEXP, SEXP rhoSEXP, SEXP tausqSEXP, SEXP sigmasqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type yy(yySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type bb(bbSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type sre(sreSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type sre_ne_label(sre_ne_labelSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const double& >::type tausq(tausqSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigmasq(sigmasqSEXP);
    rcpp_result_gen = Rcpp::wrap(updateGnnmpSPE(yy, XX, bb, sre, sre_ne_label, rho, tausq, sigmasq));
    return rcpp_result_gen;
END_RCPP
}
// updateGnnmpLabel
List updateGnnmpLabel(const arma::colvec& sre, const arma::mat& sre_ne, const int& nne, const arma::mat& weight, const arma::colvec& mu, const arma::colvec& sigmasq, const arma::mat& rho, const arma::colvec& logit_mu, const arma::colvec& logit_ka, const arma::mat& cutoff);
RcppExport SEXP _nnmp_updateGnnmpLabel(SEXP sreSEXP, SEXP sre_neSEXP, SEXP nneSEXP, SEXP weightSEXP, SEXP muSEXP, SEXP sigmasqSEXP, SEXP rhoSEXP, SEXP logit_muSEXP, SEXP logit_kaSEXP, SEXP cutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type sre(sreSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sre_ne(sre_neSEXP);
    Rcpp::traits::input_parameter< const int& >::type nne(nneSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type sigmasq(sigmasqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type logit_mu(logit_muSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type logit_ka(logit_kaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type cutoff(cutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(updateGnnmpLabel(sre, sre_ne, nne, weight, mu, sigmasq, rho, logit_mu, logit_ka, cutoff));
    return rcpp_result_gen;
END_RCPP
}
// pplcGausNNMP
arma::colvec pplcGausNNMP(const arma::colvec& yy, const arma::mat& XX, const arma::mat& bb, const arma::mat& zz, const arma::colvec& tausq, const int& weight);
RcppExport SEXP _nnmp_pplcGausNNMP(SEXP yySEXP, SEXP XXSEXP, SEXP bbSEXP, SEXP zzSEXP, SEXP tausqSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type yy(yySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type bb(bbSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type zz(zzSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type tausq(tausqSEXP);
    Rcpp::traits::input_parameter< const int& >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(pplcGausNNMP(yy, XX, bb, zz, tausq, weight));
    return rcpp_result_gen;
END_RCPP
}
// dicGausNNMP
double dicGausNNMP(const arma::colvec& yy, const arma::mat& XX, const arma::mat& bb, const arma::mat& zz, const arma::colvec& tausq);
RcppExport SEXP _nnmp_dicGausNNMP(SEXP yySEXP, SEXP XXSEXP, SEXP bbSEXP, SEXP zzSEXP, SEXP tausqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type yy(yySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type bb(bbSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type zz(zzSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type tausq(tausqSEXP);
    rcpp_result_gen = Rcpp::wrap(dicGausNNMP(yy, XX, bb, zz, tausq));
    return rcpp_result_gen;
END_RCPP
}
// updateSNnnmpLabel
List updateSNnnmpLabel(const arma::colvec& yy, const arma::mat& yy_ne, const int& nne, const arma::mat& weight, const double& la, const double& sigmasq, const arma::mat& rho, const arma::colvec& logit_mu, const arma::colvec& logit_ka, const arma::mat& cutoff);
RcppExport SEXP _nnmp_updateSNnnmpLabel(SEXP yySEXP, SEXP yy_neSEXP, SEXP nneSEXP, SEXP weightSEXP, SEXP laSEXP, SEXP sigmasqSEXP, SEXP rhoSEXP, SEXP logit_muSEXP, SEXP logit_kaSEXP, SEXP cutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type yy(yySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type yy_ne(yy_neSEXP);
    Rcpp::traits::input_parameter< const int& >::type nne(nneSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const double& >::type la(laSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigmasq(sigmasqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type logit_mu(logit_muSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type logit_ka(logit_kaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type cutoff(cutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(updateSNnnmpLabel(yy, yy_ne, nne, weight, la, sigmasq, rho, logit_mu, logit_ka, cutoff));
    return rcpp_result_gen;
END_RCPP
}
// updateSNnnmpLabel_covars_lapar
List updateSNnnmpLabel_covars_lapar(const arma::colvec& yy, const arma::mat& yy_ne, const int& nne, const arma::mat& weight, const arma::colvec& nu, const arma::mat& nu_ne, const arma::colvec& la, const arma::mat& la_ne, const double& sigmasq, const arma::mat& rho, const arma::colvec& logit_mu, const arma::colvec& logit_ka, const arma::mat& cutoff);
RcppExport SEXP _nnmp_updateSNnnmpLabel_covars_lapar(SEXP yySEXP, SEXP yy_neSEXP, SEXP nneSEXP, SEXP weightSEXP, SEXP nuSEXP, SEXP nu_neSEXP, SEXP laSEXP, SEXP la_neSEXP, SEXP sigmasqSEXP, SEXP rhoSEXP, SEXP logit_muSEXP, SEXP logit_kaSEXP, SEXP cutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type yy(yySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type yy_ne(yy_neSEXP);
    Rcpp::traits::input_parameter< const int& >::type nne(nneSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type nu_ne(nu_neSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type la(laSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type la_ne(la_neSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigmasq(sigmasqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type logit_mu(logit_muSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type logit_ka(logit_kaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type cutoff(cutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(updateSNnnmpLabel_covars_lapar(yy, yy_ne, nne, weight, nu, nu_ne, la, la_ne, sigmasq, rho, logit_mu, logit_ka, cutoff));
    return rcpp_result_gen;
END_RCPP
}
// refNe
List refNe(const int& nne, const arma::mat& coords, const arma::colvec& obs);
RcppExport SEXP _nnmp_refNe(SEXP nneSEXP, SEXP coordsSEXP, SEXP obsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type nne(nneSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type obs(obsSEXP);
    rcpp_result_gen = Rcpp::wrap(refNe(nne, coords, obs));
    return rcpp_result_gen;
END_RCPP
}
// gridNe
List gridNe(const int& nne, const arma::mat& coords, const arma::mat& grid, const arma::colvec& obs);
RcppExport SEXP _nnmp_gridNe(SEXP nneSEXP, SEXP coordsSEXP, SEXP gridSEXP, SEXP obsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type nne(nneSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type obs(obsSEXP);
    rcpp_result_gen = Rcpp::wrap(gridNe(nne, coords, grid, obs));
    return rcpp_result_gen;
END_RCPP
}
// predGausNNMP
List predGausNNMP(const int& nne, const arma::mat& XX, const arma::mat& DD, const arma::mat& bb, const arma::mat& zz, const arma::colvec& sigmasq, const arma::colvec& tausq, const arma::colvec& phi, const arma::colvec& zeta, const arma::mat& ga, const arma::colvec& kasq, const arma::mat& grid_ne_idx, const arma::mat& grid_ne_dist, const arma::colvec& probs, const bool& verbose, const int& nreport, const bool sam);
RcppExport SEXP _nnmp_predGausNNMP(SEXP nneSEXP, SEXP XXSEXP, SEXP DDSEXP, SEXP bbSEXP, SEXP zzSEXP, SEXP sigmasqSEXP, SEXP tausqSEXP, SEXP phiSEXP, SEXP zetaSEXP, SEXP gaSEXP, SEXP kasqSEXP, SEXP grid_ne_idxSEXP, SEXP grid_ne_distSEXP, SEXP probsSEXP, SEXP verboseSEXP, SEXP nreportSEXP, SEXP samSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type nne(nneSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type DD(DDSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type bb(bbSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type zz(zzSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type sigmasq(sigmasqSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type tausq(tausqSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ga(gaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type kasq(kasqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type grid_ne_idx(grid_ne_idxSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type grid_ne_dist(grid_ne_distSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type probs(probsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int& >::type nreport(nreportSEXP);
    Rcpp::traits::input_parameter< const bool >::type sam(samSEXP);
    rcpp_result_gen = Rcpp::wrap(predGausNNMP(nne, XX, DD, bb, zz, sigmasq, tausq, phi, zeta, ga, kasq, grid_ne_idx, grid_ne_dist, probs, verbose, nreport, sam));
    return rcpp_result_gen;
END_RCPP
}
// predSkewGausNNMP_simple
List predSkewGausNNMP_simple(const arma::colvec& obs, const int& nne, const arma::mat& DD, const arma::colvec& la, const arma::colvec& sigmasq, const arma::colvec phi, const arma::colvec zeta, const arma::mat& ga, const arma::colvec kasq, const arma::mat grid_ne_idx, const arma::mat grid_ne_dist, const arma::colvec& probs, const bool& verbose, const int& nreport, const bool sam);
RcppExport SEXP _nnmp_predSkewGausNNMP_simple(SEXP obsSEXP, SEXP nneSEXP, SEXP DDSEXP, SEXP laSEXP, SEXP sigmasqSEXP, SEXP phiSEXP, SEXP zetaSEXP, SEXP gaSEXP, SEXP kasqSEXP, SEXP grid_ne_idxSEXP, SEXP grid_ne_distSEXP, SEXP probsSEXP, SEXP verboseSEXP, SEXP nreportSEXP, SEXP samSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< const int& >::type nne(nneSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type DD(DDSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type la(laSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type sigmasq(sigmasqSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ga(gaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type kasq(kasqSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type grid_ne_idx(grid_ne_idxSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type grid_ne_dist(grid_ne_distSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type probs(probsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int& >::type nreport(nreportSEXP);
    Rcpp::traits::input_parameter< const bool >::type sam(samSEXP);
    rcpp_result_gen = Rcpp::wrap(predSkewGausNNMP_simple(obs, nne, DD, la, sigmasq, phi, zeta, ga, kasq, grid_ne_idx, grid_ne_dist, probs, verbose, nreport, sam));
    return rcpp_result_gen;
END_RCPP
}
// predSkewGausNNMP_covars_lapar
List predSkewGausNNMP_covars_lapar(const arma::colvec& obs, const int& nne, const arma::mat& ref_XX, const arma::mat& nonref_XX, const arma::colvec& ref_par_lables, const arma::colvec& nonref_par_labels, const arma::mat& DD, const arma::mat& bb, const arma::mat& la, const arma::colvec& sigmasq, const arma::colvec& phi, const arma::colvec& zeta, const arma::mat& ga, const arma::colvec& kasq, const arma::mat& grid_ne_idx, const arma::mat& grid_ne_dist, const arma::colvec& probs, const bool& verbose, const int& nreport, const bool sam);
RcppExport SEXP _nnmp_predSkewGausNNMP_covars_lapar(SEXP obsSEXP, SEXP nneSEXP, SEXP ref_XXSEXP, SEXP nonref_XXSEXP, SEXP ref_par_lablesSEXP, SEXP nonref_par_labelsSEXP, SEXP DDSEXP, SEXP bbSEXP, SEXP laSEXP, SEXP sigmasqSEXP, SEXP phiSEXP, SEXP zetaSEXP, SEXP gaSEXP, SEXP kasqSEXP, SEXP grid_ne_idxSEXP, SEXP grid_ne_distSEXP, SEXP probsSEXP, SEXP verboseSEXP, SEXP nreportSEXP, SEXP samSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< const int& >::type nne(nneSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ref_XX(ref_XXSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type nonref_XX(nonref_XXSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type ref_par_lables(ref_par_lablesSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type nonref_par_labels(nonref_par_labelsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type DD(DDSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type bb(bbSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type la(laSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type sigmasq(sigmasqSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ga(gaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type kasq(kasqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type grid_ne_idx(grid_ne_idxSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type grid_ne_dist(grid_ne_distSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type probs(probsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int& >::type nreport(nreportSEXP);
    Rcpp::traits::input_parameter< const bool >::type sam(samSEXP);
    rcpp_result_gen = Rcpp::wrap(predSkewGausNNMP_covars_lapar(obs, nne, ref_XX, nonref_XX, ref_par_lables, nonref_par_labels, DD, bb, la, sigmasq, phi, zeta, ga, kasq, grid_ne_idx, grid_ne_dist, probs, verbose, nreport, sam));
    return rcpp_result_gen;
END_RCPP
}
// predCopNNMP_simple
List predCopNNMP_simple(const arma::colvec& obs, const int& cop_family, const arma::cube& cop_param, const int& marg_family, const arma::mat marg_param, const int& nne, const arma::mat& DD, const arma::colvec zeta, const arma::mat& ga, const arma::colvec kasq, const arma::mat grid_ne_idx, const arma::mat grid_ne_dist, const arma::colvec& probs, const bool& verbose, const int& nreport, const bool sam);
RcppExport SEXP _nnmp_predCopNNMP_simple(SEXP obsSEXP, SEXP cop_familySEXP, SEXP cop_paramSEXP, SEXP marg_familySEXP, SEXP marg_paramSEXP, SEXP nneSEXP, SEXP DDSEXP, SEXP zetaSEXP, SEXP gaSEXP, SEXP kasqSEXP, SEXP grid_ne_idxSEXP, SEXP grid_ne_distSEXP, SEXP probsSEXP, SEXP verboseSEXP, SEXP nreportSEXP, SEXP samSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< const int& >::type cop_family(cop_familySEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type cop_param(cop_paramSEXP);
    Rcpp::traits::input_parameter< const int& >::type marg_family(marg_familySEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type marg_param(marg_paramSEXP);
    Rcpp::traits::input_parameter< const int& >::type nne(nneSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type DD(DDSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ga(gaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type kasq(kasqSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type grid_ne_idx(grid_ne_idxSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type grid_ne_dist(grid_ne_distSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type probs(probsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int& >::type nreport(nreportSEXP);
    Rcpp::traits::input_parameter< const bool >::type sam(samSEXP);
    rcpp_result_gen = Rcpp::wrap(predCopNNMP_simple(obs, cop_family, cop_param, marg_family, marg_param, nne, DD, zeta, ga, kasq, grid_ne_idx, grid_ne_dist, probs, verbose, nreport, sam));
    return rcpp_result_gen;
END_RCPP
}
// predCopNNMP_covars
List predCopNNMP_covars(const arma::colvec& obs, const int& cop_family, const arma::cube& cop_param, const int& marg_family, const arma::mat ref_marg_param1, const arma::mat ref_marg_param2, const arma::mat nonref_marg_param1, const arma::mat nonref_marg_param2, const int& nne, const arma::mat& DD, const arma::colvec zeta, const arma::mat& ga, const arma::colvec kasq, const arma::mat grid_ne_idx, const arma::mat grid_ne_dist, const arma::colvec& probs, const bool& verbose, const int& nreport, const bool sam);
RcppExport SEXP _nnmp_predCopNNMP_covars(SEXP obsSEXP, SEXP cop_familySEXP, SEXP cop_paramSEXP, SEXP marg_familySEXP, SEXP ref_marg_param1SEXP, SEXP ref_marg_param2SEXP, SEXP nonref_marg_param1SEXP, SEXP nonref_marg_param2SEXP, SEXP nneSEXP, SEXP DDSEXP, SEXP zetaSEXP, SEXP gaSEXP, SEXP kasqSEXP, SEXP grid_ne_idxSEXP, SEXP grid_ne_distSEXP, SEXP probsSEXP, SEXP verboseSEXP, SEXP nreportSEXP, SEXP samSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< const int& >::type cop_family(cop_familySEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type cop_param(cop_paramSEXP);
    Rcpp::traits::input_parameter< const int& >::type marg_family(marg_familySEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type ref_marg_param1(ref_marg_param1SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type ref_marg_param2(ref_marg_param2SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type nonref_marg_param1(nonref_marg_param1SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type nonref_marg_param2(nonref_marg_param2SEXP);
    Rcpp::traits::input_parameter< const int& >::type nne(nneSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type DD(DDSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ga(gaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type kasq(kasqSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type grid_ne_idx(grid_ne_idxSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type grid_ne_dist(grid_ne_distSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type probs(probsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int& >::type nreport(nreportSEXP);
    Rcpp::traits::input_parameter< const bool >::type sam(samSEXP);
    rcpp_result_gen = Rcpp::wrap(predCopNNMP_covars(obs, cop_family, cop_param, marg_family, ref_marg_param1, ref_marg_param2, nonref_marg_param1, nonref_marg_param2, nne, DD, zeta, ga, kasq, grid_ne_idx, grid_ne_dist, probs, verbose, nreport, sam));
    return rcpp_result_gen;
END_RCPP
}
// logitGausWeight1_cpp
List logitGausWeight1_cpp(const int& nne, const arma::mat& dist, const double& zeta, const arma::colvec& mu, const arma::colvec& ka);
RcppExport SEXP _nnmp_logitGausWeight1_cpp(SEXP nneSEXP, SEXP distSEXP, SEXP zetaSEXP, SEXP muSEXP, SEXP kaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type nne(nneSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type dist(distSEXP);
    Rcpp::traits::input_parameter< const double& >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type ka(kaSEXP);
    rcpp_result_gen = Rcpp::wrap(logitGausWeight1_cpp(nne, dist, zeta, mu, ka));
    return rcpp_result_gen;
END_RCPP
}
// logitGausWeight2_cpp
List logitGausWeight2_cpp(const int& nne, const arma::mat& dist, const double& zeta, const arma::colvec& mu, const arma::colvec& ka);
RcppExport SEXP _nnmp_logitGausWeight2_cpp(SEXP nneSEXP, SEXP distSEXP, SEXP zetaSEXP, SEXP muSEXP, SEXP kaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type nne(nneSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type dist(distSEXP);
    Rcpp::traits::input_parameter< const double& >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type ka(kaSEXP);
    rcpp_result_gen = Rcpp::wrap(logitGausWeight2_cpp(nne, dist, zeta, mu, ka));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_nnmp_dBiCopMar_cpp", (DL_FUNC) &_nnmp_dBiCopMar_cpp, 8},
    {"_nnmp_dBiCopMar_cpp2", (DL_FUNC) &_nnmp_dBiCopMar_cpp2, 8},
    {"_nnmp_updateCopLabel", (DL_FUNC) &_nnmp_updateCopLabel, 14},
    {"_nnmp_updateCeCopAux", (DL_FUNC) &_nnmp_updateCeCopAux, 8},
    {"_nnmp_cdfCopNNMP_simple", (DL_FUNC) &_nnmp_cdfCopNNMP_simple, 11},
    {"_nnmp_updateGnnmpSPE", (DL_FUNC) &_nnmp_updateGnnmpSPE, 8},
    {"_nnmp_updateGnnmpLabel", (DL_FUNC) &_nnmp_updateGnnmpLabel, 10},
    {"_nnmp_pplcGausNNMP", (DL_FUNC) &_nnmp_pplcGausNNMP, 6},
    {"_nnmp_dicGausNNMP", (DL_FUNC) &_nnmp_dicGausNNMP, 5},
    {"_nnmp_updateSNnnmpLabel", (DL_FUNC) &_nnmp_updateSNnnmpLabel, 10},
    {"_nnmp_updateSNnnmpLabel_covars_lapar", (DL_FUNC) &_nnmp_updateSNnnmpLabel_covars_lapar, 13},
    {"_nnmp_refNe", (DL_FUNC) &_nnmp_refNe, 3},
    {"_nnmp_gridNe", (DL_FUNC) &_nnmp_gridNe, 4},
    {"_nnmp_predGausNNMP", (DL_FUNC) &_nnmp_predGausNNMP, 17},
    {"_nnmp_predSkewGausNNMP_simple", (DL_FUNC) &_nnmp_predSkewGausNNMP_simple, 15},
    {"_nnmp_predSkewGausNNMP_covars_lapar", (DL_FUNC) &_nnmp_predSkewGausNNMP_covars_lapar, 20},
    {"_nnmp_predCopNNMP_simple", (DL_FUNC) &_nnmp_predCopNNMP_simple, 16},
    {"_nnmp_predCopNNMP_covars", (DL_FUNC) &_nnmp_predCopNNMP_covars, 19},
    {"_nnmp_logitGausWeight1_cpp", (DL_FUNC) &_nnmp_logitGausWeight1_cpp, 5},
    {"_nnmp_logitGausWeight2_cpp", (DL_FUNC) &_nnmp_logitGausWeight2_cpp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_nnmp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
