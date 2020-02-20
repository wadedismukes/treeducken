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
#include <string.h>

// [[Rcpp::export]]
int treeducken(std::string params_file) {
    run_treeducken(params_file);
    return 0;
}

// [[Rcpp::export]]
Rcpp::List sim_sptree_bdp(SEXP sbr_, SEXP sdr_, SEXP numbsim_, SEXP n_tips_){
    double sbr = as<double>(sbr_);
    double sdr = as<double>(sdr_);
    unsigned numbsim = as<int>(numbsim_);
    unsigned n_tips = as<int>(n_tips_);
    if(sbr <= 0.0)
        stop("'sbr' must be bigger than 0.0.");
    if(sbr < sdr)
        stop("'sbr' must be greater than 'sdr'");
    if(numbsim < 1)
        stop("'numbsim' must be more than 1");
    if(sdr < 0.0)
        stop("'sdr' must be 0.0 or greater.");
    if(n_tips < 1)
        stop("'n_tips' must be greater than 1.");
    return bdsim_species_tree(sbr, sdr, numbsim, n_tips);
}

// [[Rcpp::export]]
Rcpp::List sim_cophylo_bdp(SEXP hbr_,
                           SEXP hdr_,
                           SEXP sbr_,
                           SEXP sdr_,
                           SEXP host_exp_rate_,
                           SEXP cosp_rate_,
                           SEXP timeToSimTo_,
                           SEXP numbsim_){
    double hbr = as<double>(hbr_);
    double hdr = as<double>(hdr_);
    double sbr = as<double>(sbr_);
    double sdr = as<double>(sdr_);
    double cosp_rate = as<double>(cosp_rate_);
    double host_exp_rate = as<double>(host_exp_rate_);
    double timeToSimTo = as<double>(timeToSimTo_);
    int numbsim = as<int>(numbsim_);

    if(hbr < 0.0){
         stop("'hbr' must be positive or 0.0.");
    }
    if(hbr < hdr){
        stop("'hbr' must be greater than 'hdr'.");
    }
    if(hdr < 0.0){
        stop("'hdr' must be a positive value or 0.0.");
    }
    if(host_exp_rate < 0.0){
        stop("'host_exp_rate' must be a positive value or 0.0.");
    }
    if(numbsim < 1){
        stop("'numbsim' must be larger than 1.");
    }
    if(cosp_rate < 0.0)
        stop("'cosp_rate' must be a positive value or 0.0.");
    if(timeToSimTo < 0.0)
        stop("'timeToSimTo' must be a positive value or 0.0.");

    return sim_host_symb_treepair(hbr,
                                  hdr,
                                  sbr,
                                  sdr,
                                  host_exp_rate,
                                  cosp_rate,
                                  timeToSimTo,
                                  numbsim);
}

