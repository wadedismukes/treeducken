# locus_tree_fun <- function(tree) {
sim_multilocus_coal <- function(locus_tree,
                                effective_pop_size,
                                generation_time = 1e-6,
                                num_genes) {
    # first rescale the locus tree to be not in time but in coalescent units
    # tree is in time units so to get to generations we want

    ipp <- 1

    locus_tree_time_length <- max(ape::node.depth.edgelength(locus_tree))
    v <- effective_pop_size * generation_time
    locus_tree <- geiger::rescale(locus_tree,
                                model = "depth",
                                locus_tree_time_length * v)

    if(!(any(grep("D[A-Z]", locus_tree$node.label)))) {
        return(treeducken::sim_multispecies_coal(
                                                locus_tree,
                                                v,
                                                num_sampled_individuals = ipp,
                                                num_genes = num_genes))
    }

    # split tree into subtrees at each duplication
    loc_trees_subtree_all <- ape::subtrees(locus_tree)
    ul_trees <- unlist(loc_trees_subtree_all, recursive = FALSE)
    node_labels <- ul_trees[seq(from = 6, to = length(ul_trees), by = 7)]
    indices_of_dup <- which(node_labels$node.label != "")
    locus_trees_by_dup <- rev(loc_trees_subtree_all[c(1, indices_of_dup)])

    # add an artificial root for treeducken
    for(i in seq_len(length(locus_trees_by_dup) - 1)) {
        locus_trees_by_dup[[i]]$root.edge <- 0.0
    }

    # get the tip labels of each locus created by duplications
    num_locus <- length(locus_trees_by_dup)
    ul_dup_loc_trees <- unlist(locus_trees_by_dup, recursive = FALSE)
    tip_labels_dup_trees <- ul_dup_loc_trees[seq(from = 4, to = length(ul_dup_loc_trees), by = 8)]
    mlc_df <- list(length = length(locus_trees_by_dup))
    # loop through the locus trees simulate using MSC then prune the loci
    # that have that locus as a subtree to one tip (to simulate that that coalescent
    # event has occurred and won't affect the sims)
    for(i in seq_len(length(locus_trees_by_dup))) {
        if(is.null(locus_trees_by_dup[[i]]$root.edge)) {
            locus_trees_by_dup[[i]]$root.edge <- 0.0
        }

        mlc_df[[i]]  <- treeducken::sim_multispecies_coal(
                                                locus_trees_by_dup[[i]],
                                                v,
                                                num_sampled_individuals = ipp,
                                                num_genes = num_genes)
        # which tips to remove
        tips_to_remove <- locus_trees_by_dup[[i]]$tip.label[2:ape::Ntip(locus_trees_by_dup[[i]])]
        # figure out which loci have those tips
        subtrees_w_tips_to_drop <- which(
                                        unlist(
                                            lapply(
                                                mapply(function(x,y) x %in% y,
                                                       tip_labels_dup_trees,
                                                       tips_to_remove)
                                                , any)))
        # loci to prune
        trees_to_prune <- locus_trees_by_dup[subtrees_w_tips_to_drop]
        # prune the loci making sure to keep the fake root
        pruned_trees <- lapply(trees_to_prune,
                               function(x) ape::drop.tip(x,
                                                         tip = tips_to_remove,
                                                         rooted = TRUE))
        # update the list of locus trees
        locus_trees_by_dup[subtrees_w_tips_to_drop] <- pruned_trees
    }
    mlc_df
}
