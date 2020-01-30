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
Rcpp::List sim_locustree_bdp(SEXP species_tree_,
                          SEXP gbr_,
                          SEXP gdr_,
                          SEXP lgtr_,
                          SEXP numLoci_){
    Rcpp::List species_tree = as<Rcpp::List>(species_tree_);
    double gbr = as<double>(gbr_);
    double gdr = as<double>(gdr_);
    double lgtr = as<double>(lgtr_);
    int numLoci = as<int>(numLoci_);
    if(strcmp(species_tree.attr("class"), "phylo") != 0){
         stop("species_tree must be an object of class phylo'.");
    }
    if(gbr <= 0.0){
         stop("'gbr' must be positive or 0.0.");
    }
    if(gbr < gdr){
        stop("'gbr' must be greater than 'gdr'.");
    }
    if(gdr <= 0.0){
        stop("'gdr' must be a positive value or 0.0.");
    }
    if(lgtr <= 0.0){
        stop("'lgtr' must be a positive value or 0.0.");
    }
    if(numLoci < 1){
        stop("'numLoci' must be larger than 1.");
    }
    SpeciesTree* specTree = new SpeciesTree(species_tree);

    return sim_locus_tree(specTree, gbr, gdr, lgtr, numLoci);
}