#' Converts an object into an object of type cophylo
#'
#' @description Functions for converting either a list of three components (host_tree, symb_tree, and association_mat) into an object of class cophylo
#' Otherwise turns arguments into the cophylo object if inputting a hostTree of type `phylo`, a symbiont tree of type `phylo`, and a matrix of type eventHistory.
#'
#' @param hostTree An object of type `phylo`
#' @param symbTree An object of type `phylo`
#' @param assocMat A matrix with rows being extant symbionts and columns being extant hosts
#' @param eventHistory An optional data frame of four columns: Symbiont Index, Host Index, Event Type (see details), and Event Time
#' @return An object of type cophylo
#'
#' @details The association matrix must be with rows equal to the number of extant symbionts and columns equal to the number of extant hosts.
#' Non-zero values in this matrix indicate associations (typically this will be a matrix of just zeros and ones).
#'
#' The eventHistory parameter has four columns: Symbiont Index, Host Index, Event Type (see details), and Event Time.
#' The indexing of the first two columns should follow the indexing of the `phylo` objects `hostTree` and `symbTree`.
#' The types of events are as follows:
#' * HG - a host speciation event
#' * HL - a host extinction event
#' * C - a cospeciation event
#' * SG - a symbiont speciation event
#' * SL - a symbiont extinction event
#' * AG - an association gain between symbiont x and host y
#' * AL - an association loss between symbiont x and host y
#'
#' @seealso is.cophylo
#'
#' @examples
#' # maybe have the gopher dataset here and the fig one and convert them into the
#' # data structure. would be a nice proof of concept
#'
#TODO: add gopher dataset, and other classic cophylo datasets
as.cophylo <- function(hostTree, symbTree, assocMat, eventHistory = NULL){
    if(!(identical(class(hostTree), "phylo")))
        stop("`hostTree` input is not of class `phylo`")
    if(!(identical(class(symbTree), "phylo")))
        stop("`symbTree` input is not of class `phylo`")
    if(!(identical(class(assocMat), "matrix")))
        stop("`assocMat` input is not of class `matrix`")
    if(!is.null(eventHistory)){
        if(!(identical(class(eventHistory), "data.frame")))
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