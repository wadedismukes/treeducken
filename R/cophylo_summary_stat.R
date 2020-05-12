#' @describeIn cophylo_summary_stat Calculates the summary statistics for one index of the list of cophylogenetic objects
#'
#' @return A vector consisting of (in order) cospeciations, host speciations, host extinctions, symbiont speciations, symbiont extinctions, parafit statistic, and parafit p-value
cophylo_summary_stat_by_indx <- function(cophylo_obj, cophylo_obj_indx){
    print(cophylo_obj_indx)
    events <- cophylo_obj[[cophylo_obj_indx]]$event_history
    event_types <- levels(events$Event.Type)
    num_each_event <- tabulate(events$Event.Type)
    names(num_each_event) <- event_types
    cospeciations <- num_each_event["C"]
    host_speciations <- num_each_event["HG"]
    host_extinctions <- num_each_event["HL"]
    symbiont_speciations <- num_each_event["SG"]
    symbiont_extinctions <- num_each_event["SL"]
    host_expansions <- NULL # not sure why I don't have a specific tag for this...
    if(length(cophylo_obj[[cophylo_obj_indx]]$host_tree$tip.label) < 3 ||
       length(cophylo_obj[[cophylo_obj_indx]]$symb_tree$tip.label) < 3)
    {
        parafits <- NA
        parafit_test <- NA
    }
    else
    {
        parafits <- treeducken::parafit_stat(cophylo_obj[[cophylo_obj_indx]]$host_tree,
                                             cophylo_obj[[cophylo_obj_indx]]$symb_tree,
                                             cophylo_obj[[cophylo_obj_indx]]$association_mat)
        parafit_test <- treeducken::parafit_test(cophylo_obj[[cophylo_obj_indx]]$host_tree,
                                                 cophylo_obj[[cophylo_obj_indx]]$symb_tree,
                                                 cophylo_obj[[cophylo_obj_indx]]$association_mat,
                                                 parafits)
    }
    c(cospeciations,
      host_speciations,
      host_extinctions,
      symbiont_speciations,
      symbiont_extinctions,
     #host_expansions,
      parafits,
      parafit_test)
}
#' Calculates summary statistics for cophylogenetyic objects
#'
#' @description For cophylogenetic objects produced in treeducken via `sim_cophylo_bdp`, calculates the numbers of different events of interest. In addition, calculates and tests the ParaFit test.
#'
#' @param cophylo_obj The cophylogenetic object produced via `sim_cophylo_bdp`
#' @param cophylo_obj_indx The index with `cophylo_obj` for `cophylo_summary_stat_by_indx`
#'
#' @return A dataframe containing statistics relevant to cophylogenetic analysis
#' @examples
#'
cophylo_summary_stat <- function(cophylo_obj){
    num_cophylo_obj <- length(cophylo_obj)
    stat_df <- data.frame(matrix(0, nrow = num_cophylo_obj, ncol = 7))
    for(i in 1:num_cophylo_obj){
        stat_df[i,] <- treeducken::cophylo_summary_stat_by_indx(cophylo_obj, i)
    }
    colnames(stat_df) <- c("Cospeciations",
                           "Host.Speciations",
                           "Host.Extinctions",
                           "Symbiont.Speciations",
                           "Symbiont.Extinctions",
                        #   "Host Expansions",
                            "Parafit.Stat",
                            "Parafit P-value")
    #stat_df[is.na(stat_df)] <- 0

    stat_df
}
#' Calculate the ParafitGlobal statistic on 2 trees and their association matrix
#'
#' @description Calculate the ParafitGlobal Statistic to be used in the hypothesis test described in Legendre et al. (2002).
#' The null hypothesis of this test being that the evolution of the two trees together with their associations at the present
#' have been independent.
#' @param host_tr The host tree of class "phy"
#' @param symb_tr The symbiont tree of class "phy
#' @param assoc_mat Association matrix between the extant tips of `host_tr` and `symb_tr`
#' @details
#' `parafit_stat` drops any non-extant tips from the tree. Then the phylogenetic distance matrix is obtained for both host and symbiont tree.
#' Next the principal coordinates are found for the host and symbiont distance matrices before these PCoA vectors are used in the following
#' matrix multiplication following Legendre et al. (2002): D = H t(A) A. The trace is then found of this to get our ParaFitGlobal Statistic.
#'
#' The test function `parafit_test` performs a row-wise permutation of the association matrix as described in Legendre et al. 2002. This is
#' performed a number of times set by the user (default is 999) and a p-value is output.
#'
#' The value from this is input into the test function. Note that this gives only the raw statistic unlike `ape::parafit`. That is the
#' only reason it is implemented here in treeducken (similar to `treeducken::cherries`).
#' @examples
#' tr_pair <- sim_cophylo_bdp(hbr_=0.1,
#'                           hdr_=0.05,
#'                           sdr_=0.1,
#'                           host_exp_rate_=0.4,
#'                           sbr_=0.05,
#'                           cosp_rate_=1.0,
#'                           numbsim_=100,
#'                           timeToSimTo_=1)
#' # maybe we are interested in only cophylogenetic object 42
#' ht <- tr_pair[[42]]$host_tree
#' st <- tr_pair[[42]]$symb_tree
#' A <- tr_pair[[42]]$association_mat
#' pfs <- parafit_stat(host_tr = ht, symb_tr = st, assoc_mat = A)
#'
#' parafit_test(ht, st, A, pfs, reps = 99)
#' @seealso parafit_test
#' @references
#' Legendre, P., Y. Desdevises and E. Bazin. 2002. A statistical test for host-parasite coevolution. Systematic Biology, 51(2), 217–234.
parafit_stat <- function(host_tr, symb_tr, assoc_mat){
    host_tree <- geiger::drop.extinct(host_tr, tol= 0.001)
    symb_tree <- geiger::drop.extinct(symb_tr, tol = 0.001)
    H <- ape::cophenetic.phylo(host_tree)
    S <- ape::cophenetic.phylo(symb_tree)
    H_pcoas <- ape::pcoa(H)
    S_pcoas <- ape::pcoa(S)
# remember the weirdness with pca$values corresponding to cospeciations
    D <- t(H_pcoas$vectors) %*% t(assoc_mat) %*% S_pcoas$vectors
    sum(diag(D)^2)
}
#' @describeIn parafit_stat Perform ParaFit Hypothesis Test
#' @param D the statistic calculated using `parafit_stat`
#' @param reps Number of permutations to perform on the association matrix for the hypothesis test
#' @return A p-value for the hypothesis test described above
parafit_test <- function(host_tr, symb_tr, assoc_mat, D, reps = 999){
    null_dist <- vector(length = reps)
    for(i in 1:reps){
        shuffled_A <- t(apply(assoc_mat, 1, sample))
        null_dist[i] <- parafit_stat(host_tr, symb_tr, shuffled_A)
    }
    append(null_dist, D)
    length(null_dist[null_dist >= D]) / reps
}


almost_parafit_stat <- function(host_tr, symb_tr, assoc_mat){
  host_tree <- geiger::drop.extinct(host_tr, tol= 0.001)
  symb_tree <- geiger::drop.extinct(symb_tr, tol = 0.001)
  H <- ape::cophenetic.phylo(host_tree)
  S <- ape::cophenetic.phylo(symb_tree)
  H_eigen <- eigen(H)
  S_eigen <- eigen(S)
  # remember the weirdness with pca$values corresponding to cospeciations
  D <- t(H_eigen$vectors) %*% t(assoc_mat) %*% S_eigen$vectors
  sum(diag(D)^2)
}
