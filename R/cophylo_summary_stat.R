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


parafit_test <- function(host_tr, symb_tr, assoc_mat, D, reps = 999){
    null_dist <- vector(length = reps)
    for(i in 1:reps){
        shuffled_A <- t(apply(assoc_mat, 1, sample))
        null_dist[i] <- parafit_stat(host_tr, symb_tr, shuffled_A)
    }
    append(null_dist, D)
    length(null_dist[null_dist >= D]) / reps
}
