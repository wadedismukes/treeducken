// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
#include "Engine.h"

// [[Rcpp::export]]
int treeducken(std::string params_file) {
    run_treeducken(params_file);
    return 0;
}

// [[Rcpp::export]]
Rcpp::List sim_sptree_bdp(SEXP sbr_, SEXP sdr_, SEXP numbsim_, SEXP n_tips_){
    return bdsim_species_tree(sbr_, sdr_, numbsim_, n_tips_);
}