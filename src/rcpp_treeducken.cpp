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
//' Simulates species tree using constant rate birth-death process
//'
//' @details Forward simulates to a number of tips. This function does so using
//'     the general algorithm of Hartmann et al. 2010.
//' @param sbr_ species birth rate (i.e. speciation rate)
//' @param sdr_ species death rate (i.e. extinction rate)
//' @param numbsim_ number of species trees to simulate
//' @param n_tips_ number of tips to simulate to
//' @return List of objects of the tree class (as implemented in APE)
//' @references
//' K. Hartmann, D. Wong, T. Stadler. Sampling trees from evolutionary models.
//'     Syst. Biol., 59(4): 465-476, 2010.
//'
//' T. Stadler. Simulating trees on a fixed number of extant species.
//'     Syst. Biol., 60: 676-684, 2011.
//' @examples
//' mu <- 0.5 # death rate
//' lambda <- 2.0 # birth rate
//' numb_replicates <- 10
//' numb_extant_tips <- 4
//' # simulate trees under the GSA so first simulates a tree with
//' # numb_extant_tips * 100 tips counting each time we have a tree with 10 tips
//' # then randomly picks one of those trees
//'
//' sim_sptree_bdp(sbr_ = lambda,
//'                 sdr_ = mu,
//'                 numbsim_ = numb_replicates,
//'                 n_tips_ = numb_extant_tips)
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
//' Simulates species tree using constant rate birth-death process to a time
//'
//' @details Forward simulates a tree until a provided time is reached.
//' @param sbr_ species birth rate (i.e. speciation rate)
//' @param sdr_ species death rate (i.e. extinction rate)
//' @param numbsim_ number of species trees to simulate
//' @param t_ time to simulate to
//' @return List of objects of the tree class (as implemented in APE)
//' @references
//' K. Hartmann, D. Wong, T. Stadler. Sampling trees from evolutionary models.
//'     Syst. Biol., 59(4): 465-476, 2010.
//'
//' T. Stadler. Simulating trees on a fixed number of extant species.
//'     Syst. Biol., 60: 676-684, 2011.
//' @examples
//' mu <- 0.5 # death rate
//' lambda <- 2.0 # birth rate
//' numb_replicates <- 10
//' time <- 4
//'
//' sim_sptree_bdp(sbr_ = lambda,
//'                 sdr_ = mu,
//'                 numbsim_ = numb_replicates,
//'                 t_ = time)
// [[Rcpp::export]]
Rcpp::List sim_sptree_bdp_time(SEXP sbr_, SEXP sdr_, SEXP numbsim_, SEXP t_){
    double sbr = as<double>(sbr_);
    double sdr = as<double>(sdr_);
    unsigned numbsim = as<int>(numbsim_);
    double t = as<double>(t_);
    if(sbr <= 0.0)
        stop("'sbr' must be bigger than 0.0.");
    if(sbr < sdr)
        stop("'sbr' must be greater than 'sdr'");
    if(numbsim < 1)
        stop("'numbsim' must be more than 1");
    if(sdr < 0.0)
        stop("'sdr' must be 0.0 or greater.");
    if(t <= 0.0)
        stop("'t' must be greater than 0.");
    return sim_bdsimple_species_tree(sbr, sdr, numbsim, t);
}
//' Simulates locus tree using constant rate birth-death-transfer process
//'
//' @details Given a species tree simulates a locus or gene family tree along
//'     the species tree.
//' @param species_tree species tree to simulate along
//' @param gbr gene birth rate
//' @param gdr gene death rate
//' @param lgtr gene trasnfer rate
//' @param num_loci number of locus trees to simulate
//' @return List of objects of the tree class (as implemented in APE)
//' @references
//' Rasmussen MD, Kellis M. Unified modeling of gene duplication, loss, and
//'     coalescence using a locus tree. Genome Res. 2012;22(4):755–765.
//'     doi:10.1101/gr.123901.111
//' @examples
//' # first simulate a species tree
//' mu <- 0.5 # death rate
//' lambda <- 2.0 # birth rate
//' numb_replicates <- 10
//' numb_extant_tips <- 4
//' # simulate trees under the GSA so first simulates a tree with
//' # numb_extant_tips * 100 tips counting each time we have a tree with 10 tips
//' # then randomly picks one of those trees
//'
//' sp_tree <- sim_sptree_bdp(sbr_ = lambda,
//'                 sdr = mu,
//'                 numbsim = numb_replicates,
//'                 n_tips = numb_extant_tips)
//'
//' gene_br <- 1.0
//' gene_dr <- 0.2
//' transfer_rate <- 0.2
//' sim_locustree_bdp(species_tree = sp_tree,
//'                   gbr = gene_br,
//'                   gdr = gene_dr,
//'                   lgtr = transfer_rate,
//'                   num_loci = 10)
// [[Rcpp::export]]
Rcpp::List sim_locustree_bdp(SEXP species_tree,
                             SEXP gbr,
                             SEXP gdr,
                             SEXP lgtr,
                             SEXP num_loci){
    Rcpp::List species_tree_ = as<Rcpp::List>(species_tree);
    double gbr_ = as<double>(gbr);
    double gdr_ = as<double>(gdr);
    double lgtr_ = as<double>(lgtr);
    unsigned numLoci = as<int>(num_loci);
    if(gbr_ < 0.0)
        stop("'gbr' must be a positive number or 0.0");
    if(gbr_ < gdr_){
        if(gbr_ == gdr_)
            stop("'gbr' can only be equal to 'gdr' if both are 0.0");
        stop("'gbr' must be greater than 'gdr'");
    }
    if(numLoci < 1)
        stop("'num_loci' must be greater than 1");
    if(lgtr_ < 0.0)
        stop("'lgtr' must be a positive number or 0.0");
    if(gdr_ < 0.0)
        stop("'gdr' must be greater than or equal to 0.0");
    if(strcmp(species_tree_.attr("class"), "phylo") != 0)
        stop("species_tree must be an object of class phylo'.");
    SpeciesTree* specTree = new SpeciesTree(species_tree_);
    return sim_locus_tree(specTree, gbr_, gdr_, lgtr_, numLoci);
}
//' Simulates a cophylogenetic system using a paired birth-death process
//'
//' @details Simulates a cophylogenetic system using birth-death processes. The
//'     host tree is simulated following a constant rate birth-death process
//'     with an additional parameter - the cospeciation rate. This rate works as
//'     the speciation rate with the additional effect that if cospeciation
//'     occurs the symbiont tree also speciates. The symbiont tree is related to
//'     the host tree via an association matrix that describes which lineages
//'     are associated with which. The symbiont tree has an independent
//'     birth-death process with the addition of a host shift speciation rate
//'     that allows for the addition of more associated hosts upon symbiont
//'     speciation.
//' @param species_tree_ species tree to simulate along
//' @param hbr_ host tree birth rate
//' @param hdr_ host tree death rate
//' @param sbr_ symbiont tree birth rate
//' @param sdr_ symbiont tree death rate
//' @param host_exp_rate_ host shift speciation rate
//' @param cosp_rate_ cospeciation rate
//' @param timeToSimTo_ time units to simulate until
//' @param numbsim_ number of replicates
//' @return A list containing the `host_tree`, the `symbiont_tree`, the
//'     association matrix at present, and the history of events that have
//'     occurred.
//' @examples
//' # first simulate a species tree
//' host_mu <- 0.5 # death rate
//' host_lambda <- 2.0 # birth rate
//' numb_replicates <- 10
//' time <- 2.9
//' symb_mu <- 0.2
//' symb_lambda <- 0.4
//' host_shift_rate <- 0.1
//' cosp_rate <- 2.0
//'
//' cophylo_pair <- sim_cophylo_bdp(hbr_ = host_lambda,
//'                            hdr_ = host_mu
//'                            cosp_rate_ = cosp_rate,
//'                            host_exp_rate_ = host_shift_rate,
//'                            sdr_ = symb_mu,
//'                            sbr_ = symb_lambda,
//'                            numbsim_ = numb_replicates,
//'                            timeToSimTo_ = time)
//'
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

// [[Rcpp::export]]
Rcpp::List sim_locustree_genetree_mlc(SEXP species_tree,
                             SEXP gbr,
                             SEXP gdr,
                             SEXP lgtr,
                             SEXP num_loci,
                             SEXP num_sampled_individuals,
                             SEXP popsize,
                             SEXP num_genes_per_locus){
    Rcpp::List species_tree_ = as<Rcpp::List>(species_tree);
    double gbr_ = as<double>(gbr);
    double gdr_ = as<double>(gdr);
    double lgtr_ = as<double>(lgtr);
    unsigned numLoci = as<int>(num_loci);

    int num_sampled_individuals_ = as<int>(num_sampled_individuals);
    int popsize_ = as<int>(popsize);
    int num_genes_ = as<int>(num_genes_per_locus);

    if(gbr_ < 0.0)
        stop("'gbr' must be a positive number or 0.0");
    if(gbr_ < gdr_)
        stop("'gbr' must be greater than 'gdr'");
    if(numLoci < 1)
        stop("'num_loci' must be greater than or equal to 1");
    if(lgtr_ < 0.0)
        stop("'lgtr' must be a positive number or 0.0");
    if(gdr_ < 0.0)
        stop("'gdr' must be greater than or equal to 0.0");
    if(popsize_ < 1)
        stop("'popsize' must be greater than or equal to 1");
    if(num_genes_ < 1)
        stop("'num_genes_per_locus' must be greater than or equal to 1");
    if(num_sampled_individuals_ > popsize_)
        stop("'num_sampled_individuals' must be less than or equal to 'popsize'");
    if(num_sampled_individuals_ < 1)
        stop("'num_sampled_individuals' must be greater than or equal to 1");
    if(strcmp(species_tree_.attr("class"), "phylo") != 0)
        stop("species_tree must be an object of class phylo'.");
    SpeciesTree* specTree = new SpeciesTree(species_tree_);
    return sim_locus_tree_gene_tree(specTree,
                          gbr_,
                          gdr_,
                          lgtr_,
                          numLoci,
                          popsize_,
                          num_sampled_individuals_,
                          num_genes_);
}

