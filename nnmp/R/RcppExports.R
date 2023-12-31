# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

dBiCopMar_cpp <- function(xx, yy, cop_family, cop_param, marg_family, mar_param, ce = FALSE, logd = FALSE) {
    .Call(`_nnmp_dBiCopMar_cpp`, xx, yy, cop_family, cop_param, marg_family, mar_param, ce, logd)
}

dBiCopMar_cpp2 <- function(xx, yy, cop_family, cop_param, marg_family, mar_param, ce = FALSE, logd = FALSE) {
    .Call(`_nnmp_dBiCopMar_cpp2`, xx, yy, cop_family, cop_param, marg_family, mar_param, ce, logd)
}

updateCopLabel <- function(yy, yy_ne, nne, weight, rho, cop_family, marg_family, yy_param1, yy_param2, yy_ne_param1, yy_ne_param2, logit_mu, logit_ka, cutoff) {
    .Call(`_nnmp_updateCopLabel`, yy, yy_ne, nne, weight, rho, cop_family, marg_family, yy_param1, yy_param2, yy_ne_param1, yy_ne_param2, logit_mu, logit_ka, cutoff)
}

updateCeCopAux <- function(yy, oo, rho, ne_label, cop_family, mar_family, yy_param1, yy_param2) {
    .Call(`_nnmp_updateCeCopAux`, yy, oo, rho, ne_label, cop_family, mar_family, yy_param1, yy_param2)
}

cdfCopNNMP_simple <- function(grid_ne_obs, grid_ne_dist, cop_family, cop_param, marg_family, marg_param, zeta, ga, kasq, DD, probs) {
    .Call(`_nnmp_cdfCopNNMP_simple`, grid_ne_obs, grid_ne_dist, cop_family, cop_param, marg_family, marg_param, zeta, ga, kasq, DD, probs)
}

updateGnnmpSPE <- function(yy, XX, bb, sre, sre_ne_label, rho, tausq, sigmasq) {
    .Call(`_nnmp_updateGnnmpSPE`, yy, XX, bb, sre, sre_ne_label, rho, tausq, sigmasq)
}

updateGnnmpLabel <- function(sre, sre_ne, nne, weight, mu, sigmasq, rho, logit_mu, logit_ka, cutoff) {
    .Call(`_nnmp_updateGnnmpLabel`, sre, sre_ne, nne, weight, mu, sigmasq, rho, logit_mu, logit_ka, cutoff)
}

pplcGausNNMP <- function(yy, XX, bb, zz, tausq, weight = 1.0) {
    .Call(`_nnmp_pplcGausNNMP`, yy, XX, bb, zz, tausq, weight)
}

dicGausNNMP <- function(yy, XX, bb, zz, tausq) {
    .Call(`_nnmp_dicGausNNMP`, yy, XX, bb, zz, tausq)
}

updateSNnnmpLabel <- function(yy, yy_ne, nne, weight, la, sigmasq, rho, logit_mu, logit_ka, cutoff) {
    .Call(`_nnmp_updateSNnnmpLabel`, yy, yy_ne, nne, weight, la, sigmasq, rho, logit_mu, logit_ka, cutoff)
}

updateSNnnmpLabel_covars_lapar <- function(yy, yy_ne, nne, weight, nu, nu_ne, la, la_ne, sigmasq, rho, logit_mu, logit_ka, cutoff) {
    .Call(`_nnmp_updateSNnnmpLabel_covars_lapar`, yy, yy_ne, nne, weight, nu, nu_ne, la, la_ne, sigmasq, rho, logit_mu, logit_ka, cutoff)
}

refNe <- function(nne, coords, obs) {
    .Call(`_nnmp_refNe`, nne, coords, obs)
}

gridNe <- function(nne, coords, grid, obs) {
    .Call(`_nnmp_gridNe`, nne, coords, grid, obs)
}

predGausNNMP <- function(nne, XX, DD, bb, zz, sigmasq, tausq, phi, zeta, ga, kasq, grid_ne_idx, grid_ne_dist, probs, verbose, nreport, sam = FALSE) {
    .Call(`_nnmp_predGausNNMP`, nne, XX, DD, bb, zz, sigmasq, tausq, phi, zeta, ga, kasq, grid_ne_idx, grid_ne_dist, probs, verbose, nreport, sam)
}

predSkewGausNNMP_simple <- function(obs, nne, DD, la, sigmasq, phi, zeta, ga, kasq, grid_ne_idx, grid_ne_dist, probs, verbose, nreport, sam = FALSE) {
    .Call(`_nnmp_predSkewGausNNMP_simple`, obs, nne, DD, la, sigmasq, phi, zeta, ga, kasq, grid_ne_idx, grid_ne_dist, probs, verbose, nreport, sam)
}

predSkewGausNNMP_covars_lapar <- function(obs, nne, ref_XX, nonref_XX, ref_par_lables, nonref_par_labels, DD, bb, la, sigmasq, phi, zeta, ga, kasq, grid_ne_idx, grid_ne_dist, probs, verbose, nreport, sam = FALSE) {
    .Call(`_nnmp_predSkewGausNNMP_covars_lapar`, obs, nne, ref_XX, nonref_XX, ref_par_lables, nonref_par_labels, DD, bb, la, sigmasq, phi, zeta, ga, kasq, grid_ne_idx, grid_ne_dist, probs, verbose, nreport, sam)
}

predCopNNMP_simple <- function(obs, cop_family, cop_param, marg_family, marg_param, nne, DD, zeta, ga, kasq, grid_ne_idx, grid_ne_dist, probs, verbose, nreport, sam = FALSE) {
    .Call(`_nnmp_predCopNNMP_simple`, obs, cop_family, cop_param, marg_family, marg_param, nne, DD, zeta, ga, kasq, grid_ne_idx, grid_ne_dist, probs, verbose, nreport, sam)
}

predCopNNMP_covars <- function(obs, cop_family, cop_param, marg_family, ref_marg_param1, ref_marg_param2, nonref_marg_param1, nonref_marg_param2, nne, DD, zeta, ga, kasq, grid_ne_idx, grid_ne_dist, probs, verbose, nreport, sam = FALSE) {
    .Call(`_nnmp_predCopNNMP_covars`, obs, cop_family, cop_param, marg_family, ref_marg_param1, ref_marg_param2, nonref_marg_param1, nonref_marg_param2, nne, DD, zeta, ga, kasq, grid_ne_idx, grid_ne_dist, probs, verbose, nreport, sam)
}

logitGausWeight1_cpp <- function(nne, dist, zeta, mu, ka) {
    .Call(`_nnmp_logitGausWeight1_cpp`, nne, dist, zeta, mu, ka)
}

logitGausWeight2_cpp <- function(nne, dist, zeta, mu, ka) {
    .Call(`_nnmp_logitGausWeight2_cpp`, nne, dist, zeta, mu, ka)
}

