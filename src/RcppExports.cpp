// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// treeducken
int treeducken(std::string params_file);
RcppExport SEXP _treeducken_treeducken(SEXP params_fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type params_file(params_fileSEXP);
    rcpp_result_gen = Rcpp::wrap(treeducken(params_file));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_treeducken_treeducken", (DL_FUNC) &_treeducken_treeducken, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_treeducken(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
