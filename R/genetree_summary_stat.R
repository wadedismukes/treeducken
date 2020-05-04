

genetree_summary_stat <- function(locus_tree_gene_tree_obj, locus_tree_indx){
     genetrees <- locus_tree_gene_tree_obj[[locus_tree_indx]]$gene.trees
     locus_tree <- locus_tree_gene_tree_obj[[locus_tree_indx]]$locus.tree
     num_genetrees <- length(genetrees)
     colless <- vector(length = num_genetrees)
     gamma_locus <- rep(ape::gammaStat(locus_tree), times = length(genetrees))
     gamma <- vector(length = num_genetrees)
     sackin <- vector(length = num_genetrees)
     tmrca <- vector(length = num_genetrees)
     for(i in 1:num_genetrees){
         colless[i] <- apTreeshape::colless(genetrees[[i]])
         gamma[i] <- ape::gammaStat(genetrees[[i]])
         sackin[i] <- apTreeshape::sackin(genetrees[[i]])
         tmrca[i] <- max(phytools::nodeHeights(genetrees[[i]]))
     }
     data.frame(colless, sackin, tmrca, gamma_locus, gamma)
}


