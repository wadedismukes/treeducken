#' Plot host and symbiont pair with current associations
#'
#' This function plots a host and symbiont tree given the object returned by
#' `sim_cophylo_bdp`.
#'
#' @param tr_pair_obj a tree pair object returned by `sim_cophylo_bdp`
#' @return a plot of the host and symbiont tree with extant interactions
#' @examples
#'
#' tr_pair <- sim_cophylo_bdp(hbr_=0.1,
#'                            hdr_=0.05,
#'                            sdr_=0.1,
#'                            host_exp_rate_=0.4,
#'                            sbr_=0.05,
#'                            cosp_rate_=1.0,
#'                            numbsim_=10,
#'                            timeToSimTo_=4)
#' plot_paired_trees_with_associations(tr_pair[[1]])
plot_paired_trees_with_associations <- function(tr_pair_obj){
    host_tree <- geiger::drop.extinct(tr_pair_obj$host_tree, tol= 0.001)
    symb_tree <- geiger::drop.extinct(tr_pair_obj$symb_tree, tol = 0.001)
    assoc_mat <- tr_pair_obj$association_mat
    rownames(assoc_mat) <- symb_tree$tip.label
    colnames(assoc_mat) <- host_tree$tip.label

    m <- apply(assoc_mat, 2, function(x) which(x == 1))
    tip_label_tab <- matrix(ncol=2)
    for(i in 1:length(m)){
        for(j in 1:length(m[i][[1]])){
            r <- c(names(m[i]), names(m[i][[1]][j]))
            tip_label_tab <- rbind(tip_label_tab, r)
        }
    }
    tip_label_table <- tip_label_tab[-c(1),]

    plot(obj <- phytools::cophylo(host_tree,
                        symb_tree,
                        assoc = tip_label_table))

}
