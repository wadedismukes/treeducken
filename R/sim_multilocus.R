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

    # split tree into subtrees at each duplication
    locus_trees_by_dup <- seperate_into_loci(locus_tree)

    mlc_df <- list(length = length(locus_trees_by_dup))
    # loop through the locus trees simulate using MSC then prune the loci
    # that have that locus as a subtree to one tip (to simulate that
    # that coalescent
    # event has occurred and won't affect the sims)
    for(i in seq_len(length(locus_trees_by_dup))) {
        # add artificial root for treeducken
        if(is.null(locus_trees_by_dup[[i]]$root.edge)) {
            locus_trees_by_dup[[i]]$root.edge <- 0.0
        }

        mlc_df[[i]]  <- treeducken::sim_multispecies_coal(
                                                locus_trees_by_dup[[i]],
                                                v,
                                                num_sampled_individuals = ipp,
                                                num_genes = num_genes)
        # collapse those loci on all
        locus_trees_by_dup <- collapse_locus_subtree(locus_trees_by_dup,
                                                    locus_trees_by_dup[[i]])
    }
    mlc_df
}


seperate_into_loci <- function(locus_tree) {
    if(!(any(grep("D[A-Z]", locus_tree$node.label)))) {
        return(locus_tree)
    }

    # split tree into subtrees at each duplication
    loc_trees_subtree_all <- ape::subtrees(locus_tree)
    ul_trees <- unlist(loc_trees_subtree_all, recursive = FALSE)
    node_labels <- ul_trees[which(names(ul_trees) == "node.label")]
    indices_of_dup <- which(node_labels$node.label != "")
    rev(loc_trees_subtree_all[c(1, indices_of_dup)])

}
# SUBTREES
# given locus tree with duplications points marked with D + whatever
# give daughter trees + mother tree

# LOOP THROUGH SUBTREES


# COALESCE


# COLLAPSE THAT SUBTREE ON TREES CONTAINING THAT SUBTREE
collapse_locus_subtree <- function(list_of_subtrees,
                                   locus_to_collapse) {
    tip_labels_subtrees <- get_tip_labels_tree_list(list_of_subtrees)
    tips_to_remove <- locus_to_collapse$tip.label[-1]
        # figure out which loci have those tips

    tree_indices_to_drop <- lapply(tip_labels_subtrees,
                                   function(x) x %in% tips_to_remove)

    self_check <-  which(
                        unlist(
                            lapply(tree_indices_to_drop, all))
                    , useNames = FALSE)
    if(length(self_check) != 0)
        tree_indices_to_drop <- tree_indices_to_drop[-self_check]

    tree_indices_to_drop <- which(
                                unlist(
                                    lapply(tree_indices_to_drop, any
                            , useName = FALSE)))
    trees_to_prune <- list_of_subtrees[tree_indices_to_drop]
    pruned_trees <- suppressWarnings(lapply(trees_to_prune,
                                function(x) ape::drop.tip(x,
                                                         tip = tips_to_remove,
                                                         rooted = TRUE)))
    list_of_subtrees[tree_indices_to_drop] <- pruned_trees
    list_of_subtrees
}

get_tip_labels <- function(tree) {
    tree$tip.label
}

get_tip_labels_tree_list <- function(multi_tree) {
    ul_multi_tree <- unlist(multi_tree, recursive = FALSE)
    ul_multi_tree[which(names(ul_multi_tree) == "tip.label")]
}
