// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// sim_stBD
Rcpp::List sim_stBD(SEXP sbr, SEXP sdr, SEXP numbsim, Rcpp::NumericVector n_tips, Rcpp::NumericVector gsa_stop_mult);
RcppExport SEXP _treeducken_sim_stBD(SEXP sbrSEXP, SEXP sdrSEXP, SEXP numbsimSEXP, SEXP n_tipsSEXP, SEXP gsa_stop_multSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sbr(sbrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sdr(sdrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type numbsim(numbsimSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type n_tips(n_tipsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type gsa_stop_mult(gsa_stop_multSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_stBD(sbr, sdr, numbsim, n_tips, gsa_stop_mult));
    return rcpp_result_gen;
END_RCPP
}
// sim_stBD_t
Rcpp::List sim_stBD_t(SEXP sbr, SEXP sdr, SEXP numbsim, SEXP t);
RcppExport SEXP _treeducken_sim_stBD_t(SEXP sbrSEXP, SEXP sdrSEXP, SEXP numbsimSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sbr(sbrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sdr(sdrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type numbsim(numbsimSEXP);
    Rcpp::traits::input_parameter< SEXP >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_stBD_t(sbr, sdr, numbsim, t));
    return rcpp_result_gen;
END_RCPP
}
// sim_ltBD
Rcpp::List sim_ltBD(Rcpp::List species_tree, SEXP gbr, SEXP gdr, SEXP lgtr, SEXP num_loci, Rcpp::String transfer_type);
RcppExport SEXP _treeducken_sim_ltBD(SEXP species_treeSEXP, SEXP gbrSEXP, SEXP gdrSEXP, SEXP lgtrSEXP, SEXP num_lociSEXP, SEXP transfer_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type species_tree(species_treeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type gbr(gbrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type gdr(gdrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type lgtr(lgtrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type num_loci(num_lociSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type transfer_type(transfer_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_ltBD(species_tree, gbr, gdr, lgtr, num_loci, transfer_type));
    return rcpp_result_gen;
END_RCPP
}
// sim_cophyBD_ana
Rcpp::List sim_cophyBD_ana(SEXP hbr, SEXP hdr, SEXP sbr, SEXP sdr, SEXP s_disp_r, SEXP s_extp_r, SEXP host_exp_rate, SEXP cosp_rate, SEXP time_to_sim, SEXP numbsim, Rcpp::NumericVector host_limit, Rcpp::CharacterVector hs_mode);
RcppExport SEXP _treeducken_sim_cophyBD_ana(SEXP hbrSEXP, SEXP hdrSEXP, SEXP sbrSEXP, SEXP sdrSEXP, SEXP s_disp_rSEXP, SEXP s_extp_rSEXP, SEXP host_exp_rateSEXP, SEXP cosp_rateSEXP, SEXP time_to_simSEXP, SEXP numbsimSEXP, SEXP host_limitSEXP, SEXP hs_modeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type hbr(hbrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type hdr(hdrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sbr(sbrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sdr(sdrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type s_disp_r(s_disp_rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type s_extp_r(s_extp_rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type host_exp_rate(host_exp_rateSEXP);
    Rcpp::traits::input_parameter< SEXP >::type cosp_rate(cosp_rateSEXP);
    Rcpp::traits::input_parameter< SEXP >::type time_to_sim(time_to_simSEXP);
    Rcpp::traits::input_parameter< SEXP >::type numbsim(numbsimSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type host_limit(host_limitSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type hs_mode(hs_modeSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_cophyBD_ana(hbr, hdr, sbr, sdr, s_disp_r, s_extp_r, host_exp_rate, cosp_rate, time_to_sim, numbsim, host_limit, hs_mode));
    return rcpp_result_gen;
END_RCPP
}
// sim_cophyBD
Rcpp::List sim_cophyBD(SEXP hbr, SEXP hdr, SEXP sbr, SEXP sdr, SEXP host_exp_rate, SEXP cosp_rate, SEXP time_to_sim, SEXP numbsim, Rcpp::NumericVector host_limit, Rcpp::CharacterVector hs_mode, Rcpp::LogicalVector mutualism);
RcppExport SEXP _treeducken_sim_cophyBD(SEXP hbrSEXP, SEXP hdrSEXP, SEXP sbrSEXP, SEXP sdrSEXP, SEXP host_exp_rateSEXP, SEXP cosp_rateSEXP, SEXP time_to_simSEXP, SEXP numbsimSEXP, SEXP host_limitSEXP, SEXP hs_modeSEXP, SEXP mutualismSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type hbr(hbrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type hdr(hdrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sbr(sbrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sdr(sdrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type host_exp_rate(host_exp_rateSEXP);
    Rcpp::traits::input_parameter< SEXP >::type cosp_rate(cosp_rateSEXP);
    Rcpp::traits::input_parameter< SEXP >::type time_to_sim(time_to_simSEXP);
    Rcpp::traits::input_parameter< SEXP >::type numbsim(numbsimSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type host_limit(host_limitSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type hs_mode(hs_modeSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type mutualism(mutualismSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_cophyBD(hbr, hdr, sbr, sdr, host_exp_rate, cosp_rate, time_to_sim, numbsim, host_limit, hs_mode, mutualism));
    return rcpp_result_gen;
END_RCPP
}
// sim_msc
Rcpp::List sim_msc(SEXP species_tree, SEXP ne, SEXP num_sampled_individuals, SEXP num_genes, Rcpp::LogicalVector rescale, Rcpp::NumericVector mutation_rate, Rcpp::NumericVector generation_time);
RcppExport SEXP _treeducken_sim_msc(SEXP species_treeSEXP, SEXP neSEXP, SEXP num_sampled_individualsSEXP, SEXP num_genesSEXP, SEXP rescaleSEXP, SEXP mutation_rateSEXP, SEXP generation_timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type species_tree(species_treeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ne(neSEXP);
    Rcpp::traits::input_parameter< SEXP >::type num_sampled_individuals(num_sampled_individualsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type num_genes(num_genesSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type rescale(rescaleSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mutation_rate(mutation_rateSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type generation_time(generation_timeSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_msc(species_tree, ne, num_sampled_individuals, num_genes, rescale, mutation_rate, generation_time));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_treeducken_sim_stBD", (DL_FUNC) &_treeducken_sim_stBD, 5},
    {"_treeducken_sim_stBD_t", (DL_FUNC) &_treeducken_sim_stBD_t, 4},
    {"_treeducken_sim_ltBD", (DL_FUNC) &_treeducken_sim_ltBD, 6},
    {"_treeducken_sim_cophyBD_ana", (DL_FUNC) &_treeducken_sim_cophyBD_ana, 12},
    {"_treeducken_sim_cophyBD", (DL_FUNC) &_treeducken_sim_cophyBD, 11},
    {"_treeducken_sim_msc", (DL_FUNC) &_treeducken_sim_msc, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_treeducken(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
