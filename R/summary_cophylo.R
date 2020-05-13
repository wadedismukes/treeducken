## This file was modified from 'summary.phylo.R' (2019-01-30)  in the R-package
## `ape' version 5.3 with copyright 2003-2019 Emmanuel Paradis, 2006 Ben Bolker,
## and Klaus Schliep 2016 under a GNU Public Lisence 2.0.
##
## I modified this code here for my own uses.
## Distributed under GNU public license (>=2)
## Copyright 2020 Wade Dismukes
##
## Print Summary of a Cophylogeny set and "multiCophylo" operators
host_tree <- function(phy) UseMethod("host_tree")

host_tree.cophylo <- function(phy) ape::print.phylo(phy)


host_tree.multiCophylo <- function(phy){
    sapply(unclass(phy), host_tree.cophylo)
}

symb_tree <- function(phy) UseMethod("symb_tree")

symb_tree.cophylo <- function(phy) ape::print.phylo(phy)


symb_tree.multiCophylo <- function(phy){
    sapply(unclass(phy), symb_tree.cophylo)
}

association_mat <- function(cophy) UseMethod("association_mat")

association_mat.cophylo <- function(cophy) dim(cophy$association_mat)

association_mat.multiCophylo <- function(cophy) {
    sapply(unclass(cophy), association_mat.cophylo)
}

event_history <- function(cophy) UseMethod("event_history")

event_history.cophylo <- function(cophy) summary(cophy$event_history)

event_history.multiCophylo <- function(cophy){
    sapply(unclass(cophy), event_history.cophylo)
}


summary.cophylo <- function(object, ...){
    cat("\nSet of host and symbiont tree:", deparse(substitute(object)), "\n\n")
    # below is a modified version of ape's summary.phylo function
    ht_tips <- length(object$host_tree$tip.label)
    ht_nodes <- object$host_tree$Nnode
    cat("\nHost tree: ")
    cat("  Number of tips:", ht_tips, "\n")
    cat("  Number of nodes:", ht_nodes, "\n")
    if (is.null(object$host_tree$edge.length))
        cat("  No branch lengths.\n")
    else {
        cat("  Branch lengths:\n")
        print(summary(object$host_tree$edge.length)[-4])
    }
    if (is.null(object$host_tree$root.edge))
        cat("  No root edge.\n")
    else
        cat("  Root edge:", object$host_tree$root.edge, "\n")
    if (ht_tips <= 10) {
        cat("  Tip labels:", object$host_tree$tip.label[1], "\n")
        cat(paste("             ", object$host_tree$tip.label[-1]), sep = "\n")
    }
    else {
        cat("  First ten tip labels:", object$host_tree$tip.label[1], "\n")
        cat(paste("                       ", object$host_tree$tip.label[2:10]), sep = "\n")
    }
    if (is.null(object$node.label))
        cat("  No node labels.\n")
    else {
        if (ht_nodes <= 10) {
            cat("  Node labels:", object$host_tree$node.label[1], "\n")
            cat(paste("              ", object$host_tree$node.label[-1]), sep = "\n")
        }
        else {
            cat("  First ten node labels:", object$host_tree$node.label[1], "\n")
            cat(paste("                        ", object$host_tree$node.label[2:10]), sep = "\n")

        }
    }
    # symbiont tree summary
    st_tips <- length(object$symb_tree$tip.label)
    st_nodes <- object$symb_tree$Nnode
    cat("\n\nSymb tree: ")
    cat("  Number of tips:", st_tips, "\n")
    cat("  Number of nodes:", st_nodes, "\n")
    if (is.null(object$symb_tree$edge.length))
        cat("  No branch lengths.\n")
    else {
        cat("  Branch lengths:\n")
        print(summary(object$symb_tree$edge.length)[-4])
    }
    if (is.null(object$symb_tree$root.edge))
        cat("  No root edge.\n")
    else
        cat("  Root edge:", object$symb_tree$root.edge, "\n")
    if (st_tips <= 10) {
        cat("  Tip labels:", object$symb_tree$tip.label[1], "\n")
        cat(paste("             ", object$symb_tree$tip.label[-1]), sep = "\n")
    }
    else {
        cat("  First ten tip labels:", object$symb_tree$tip.label[1], "\n")
        cat(paste("                       ", object$symb_tree$tip.label[2:10]), sep = "\n")
    }
    if (is.null(object$node.label))
        cat("  No node labels.\n")
    else {
        if (st_nodes <= 10) {
            cat("  Node labels:", object$symb_tree$node.label[1], "\n")
            cat(paste("              ", object$symb_tree$node.label[-1]), sep = "\n")
        }
        else {
            cat("  First ten node labels:", object$symb_tree$node.label[1], "\n")
            cat(paste("                        ", object$symb_tree$node.label[2:10]), sep = "\n")

        }
    }


    cat("\n\nAssociation Matrix")
    cat("\n There are ", nrow(object$association_mat), " rows (i.e. extant symbionts.")
    cat("\n There are ", ncol(object$association_mat), " cols (i.e. extant hosts.")

    if(is.null(object$event_history)){
        cat("\n No event history.\n")
    }
    else{
        cat("Event history summary: \n\n", summary(object$event_history))
    }
}

print.cophylo <- function(x, printlen = 24, ...){
    cat("\n Host Tree:\n\n")
    ape::print.phylo(x$host_tree)
    cat("\n")
    cat("\n Symbiont Tree:\n\n")
    ape::print.phylo(x$symb_tree)
    cat("\n Association Matrix: \n\n")
    cat("\n There are ", nrow(x$association_mat), " rows (i.e. extant symbionts).")
    cat("\n There are ", ncol(x$association_mat), " cols (i.e. extant hosts).")
}


print.multiCophylo <- function(c, details = FALSE, ...){
    num_cophys <- length(x)
    cat(num_cophys, "cophylogenetic", ifelse(num_cophys > 1, "sets.\n", "set.\n"))
    if(details){
        for(i in 1:num_cophys){
            cat("Cophylogenetic set", i, ":\n")
            cat("\tSymbiont tree has ", length(x[[i]]$symb_tree$tip.label), " tips.\n")
            cat("\tHost tree has ", length(x[[i]]$host_tree$tip.label), " tips.\n")
        }
    }
}

`$.multiCophylo` <- function(x, name) x[[name]]

"[.multiCophylo" <- function(x, i)
{
    oc <- oldClass(x)
    class(x) <- NULL
    structure(x[i], host_tree = attr(x, "host_tree"),
              class = oc)
}

str.multiCophylo <- function(object, ...)
{
    class(object) <- NULL
    cat('Class "multiCophylo"\n')
    str(object, ...)
}


.c_cophylo_single <- function(cophy) structure(list(cophy), class = "multiCophylo")

c.cophylo <- function(..., recursive = TRUE)
{
    obj <- list(...)
    classes <- lapply(obj, class)
    iscophylo <- sapply(classes, function(x) "cophylo" %in% x)
    if (all(iscophylo)) {
        class(obj) <- "multiCophylo"
        return(obj)
    }
    if (!recursive) return(obj)
    ismulti <- sapply(classes, function(x) "multiCophylo" %in% x)
    if (all(iscophylo | ismulti)) {
        for (i in which(iscophylo)) obj[[i]] <- .c_cophylo_single(obj[[i]])
        ## added by Klaus:
        for (i in which(ismulti)) obj[[i]] <- .uncompressSetTipLabel(obj[[i]])
        obj <- .makeMultiCophyloFromObj(obj)
    } else {
        warning('some objects not of class "cophylo" or "multiCophylo": argument recursive=TRUE ignored')
    }
    obj
}

# modified from function by Klaus Schliep
.makeMultiCophyloFromObj <- function(obj)
{
    n <- length(obj)
    N <- lengths(obj, FALSE)
    cs <- c(0, cumsum(N))
    x <- vector("list", cs[length(cs)])
    for (i in 1:n) {
        a <- cs[i] + 1L
        b <- cs[i + 1L]
        x[a:b] <- obj[[i]]
    }
    class(x) <- "multiCophylo"
    x
}

c.multiCophylo <- function(..., recursive = TRUE)
{
    obj <- list(...)
    if (!recursive) return(obj)
    classes <- lapply(obj, class)
    iscophylo <- sapply(classes, function(x) "cophylo" %in% x)
    ismulti <- sapply(classes, function(x) "multiCophylo" %in% x)
    if (!all(iscophylo | ismulti)) {
        warning('some objects not of class "cophylo" or "multiPhylo": argument recursive=TRUE ignored')
        return(obj)
    }
    for (i in which(iscophylo)) obj[[i]] <- .c_cophylo_single(obj[[i]])
    ## added by Klaus
    for (i in which(ismulti)) obj[[i]] <- .uncompressSetTipLabel(obj[[i]])
    .makeMultiCophyloFromObj(obj)
}

.uncompressSetTipLabel <- function(x)
{
    HostLab <- attr(x$host_tree, "TipLabel")
    SymbLab <- attr(x$symb_tree, "TipLabel")
    if (is.null(HostLab)) return(x)
    if (is.null(SymbLab)) return(x)
    class(x) <- NULL
    for (i in 1:length(x)){
        x[[i]]$host_tree$tip.label <- HostLab
        x[[i]]$symb_tree$tip.label <- SymbLab
    }
    class(x) <- "multiCophylo"
    attr(x$host_Tree, "TipLabel") <- NULL
    attr(x$symb_tree, "TipLabel") <- NULL
    x
}

`[<-.multiCophylo` <- function(x, ..., value)
{
    ## recycling is allowed so no need to check: length(value) != length(..1)

    ## check that all elements in 'value' inherit class "phylo"
    test <- unlist(lapply(value, function(xx) !inherits(xx, "cophylo")))
    if (any(test))
        stop("at least one element in 'value' is not of class 'cophylo'.")

    oc <- oldClass(x)
    class(x) <- NULL

    if (is.null(attr(x, "host_tree"))) {
        x[..1] <- value
        class(x) <- oc
        return(x)
    }

    if (is.null(attr(x, "symb_tree"))) {
        x[..1] <- value
        class(x) <- oc
        return(x)
    }

    x[..1] <- 0L # in case x needs to be elongated
    class(x) <- oc
    j <- 1L
    for (i in ..1) {
        ## x is of class "multiPhylo", so this uses the operator below:
        x[[i]] <- value[[j]]
        j <- j + 1L
    }
    x
}

"[[.multiCophylo" <- function(x, i)
{
    class(x) <- NULL
    cophy <- x[[i]]
    if (!is.null(attr(x, "host_tree")))
        cophy$host_tree <- attr(x, "host_tree")
    if (!is.null(attr(x, "symb_tree")))
        cophy$symb_tree <- attr(x, "symb_tree")
    if (!is.null(attr(x, "association_mat")))
        cophy$association_mat <- attr(x, "association_mat")
    cophy
}

`[[<-.multiCophylo` <- function(x, ..., value)
{
    if (!inherits(value, "cophylo"))
        stop('trying to assign an object not of class "cophylo" into an object of class "multiCophylo".')

    oc <- oldClass(x)
    class(x) <- NULL

    ht <- attr(x, "host_tree")

    if (!is.null(ht)) {
        n <- length(ht$tip.label)
        if (n != length(value$host_tree$tip.label))
            stop("tree with different number of tips than those in the list (which all have the same labels; maybe you want to uncompress them)")

        o <- match(value$host_tree$tip.label, ht$tip.label)
        if (any(is.na(o)))
            stop("tree tip labels do not match with those in the list; maybe you want to uncompress them.")
        value$host_tree$tip.label <- NULL
        ie <- match(o, value$host_tree$edge[, 2])
        value$host_tree$edge[ie, 2] <- 1:n
    }

    st <- attr(x, "host_tree")

    if (!is.null(st)) {
        n <- length(st$tip.label)
        if (n != length(value$symb_tree$tip.label))
            stop("Host tree with different number of tips than those in the list (which all have the same labels; maybe you want to uncompress them)")

        o <- match(value$symb_tree$tip.label, st$tip.label)
        if (any(is.na(o)))
            stop("tree tip labels do not match with those in the list; maybe you want to uncompress them.")
        value$symb_tree$tip.label <- NULL
        ie <- match(o, value$symb_tree$edge[, 2])
        value$symb_tree$edge[ie, 2] <- 1:n
    }
    x[[..1]] <- value
    class(x) <- oc
    x
}

`$<-.multiCophylo` <- function(x, ..., value)
{
    x[[..1]] <- value
    x
}