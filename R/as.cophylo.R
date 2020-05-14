#'
#'
#'
#'
as.cophylo <- function(hostTree, symbTree, assocMat, eventHistory = NULL){
    if(!(identical(class(hostTree), "phylo")))
        stop("`hostTree` input is not of class `phylo`")
    if(!(identical(class(symbTree), "phylo")))
        stop("`symbTree` input is not of class `phylo`")
    if(!(identical(class(assocMat), "matrix")))
        stop("`assocMat` input is not of class `matrix`")
    if(!is.null(eventHistory)){
        if(!(identical(class(eventHistory), "matrix")))
            stop("`eventHistory` input is not of class `data.frame`")
    }
    pruned_host_tree <- geiger::drop.extinct(hostTree, tol = 0.0001)
    pruned_symb_tree <- geiger::drop.extinct(symbTree, tol = 0.0001)
    nExtHostTips <- length(pruned_host_tree$tip.label)
    nExtSymbTips <- length(pruned_symb_tree$tip.label)
    if(nExtHostTips != ncol(assocMat))
        stop("number of extant tips in 'hostTree' does not match number cols in 'assocMat'")
    if(nExtSymbTips != nrow(assocMat))
        stop("number of extant tips in 'symbTree' does not match number rows in 'assocMat'")
    cophylo <- list("host_tree" = hostTree,
                    "symb_tree" = symbTree,
                    "association_mat" = assocMat,
                    "event_history" = eventHistory)
    class(cophylo) <- "cophylo"
    cophylo
}


as.cophylo <- function(x){
    if (identical(class(x), "cophylo")) return(x)
    UseMethod("as.cophylo")
}