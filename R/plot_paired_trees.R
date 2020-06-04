#' Plot host and symbiont pair with current associations
#'
#' This function plots a host and symbiont tree given the object returned by
#' `sim_cophylo_bdp`.
#'
#' @param x a tree pair object returned by `sim_cophylo_bdp`
#' @param ... other plotting parameters
#' @return a plot of the host and symbiont tree with extant interactions
#' @examples
#'
#' tr_pair <- sim_cophylo_bdp(hbr = 0.1,
#'
#'                            hdr = 0.05,
#'                            sdr = 0.1,
#'                            host_exp_rate = 0.4,
#'                            sbr = 0.05,
#'                            cosp_rate =1.0,
#'                            numbsim = 10,
#'                            time_to_sim = 2)
#' plot.cophylo(tr_pair[[1]])
plot.cophylo <- function(x, ...){
    if (!inherits(x,"cophylo"))
        stop("x should be an object of class 'cophylo'.")
    # plotting parameters
    plot.new()
    if (hasArg(mar))
        mar <- list()$mar
    else
        mar <- c(0,0,0,0)
    if (hasArg(xlim))
        xlim <- list()$xlim
    else
        xlim <- c(-0.5,0.5)
    if (hasArg(scale.bar))
        scale.bar <- list()$scale.bar
    else
        scale.bar <- rep(0,2)
    if (hasArg(ylim))
        ylim <- list()$ylim
    else
        ylim <- if (any(scale.bar > 0)) c(-0.1,1) else c(0,1)
    if (hasArg(link.type))
        link.type <- list()$link.type
    else
        link.type <- "straight"
    if (hasArg(link.lwd))
        link.lwd <- list()$link.lwd
    else
        link.lwd <- 1
    if (hasArg(link.col))
        link.col <- list()$link.col
    else
        link.col <- "red"
    if (hasArg(link.lty))
        link.lty <- list()$link.lty
    else
        link.lty <- "dashed"
    if (hasArg(edge.col))
        edge.col <- list()$edge.col
    else
        edge.col <- list(left = rep("black",nrow(x$host_tree$edge)),
                         right = rep("black",nrow(x$symb_tree$edge)))
    obj <- list(...)
    if (is.null(obj$part))
        obj$part <- 0.4
    par(mar = mar)
    plot.window(xlim = xlim, ylim = ylim)
    leftArgs <- rightArgs <- obj
    leftArgs$edge.col <- edge.col$left
    rightArgs$edge.col <- edge.col$right
    if (!is.null(obj$fsize)) {
        if (length(obj$fsize) > 1) {
            leftArgs$fsize <- obj$fsize[1]
            rightArgs$fsize <- obj$fsize[2]
            sb.fsize <- ifelse(length(obj$fsize) > 2, obj$fsize[3], 1)
        }
        else
            sb.fsize <- 1
    }
    else
        sb.fsize <- 1
    # end of plotting parameters setup
    # now want to plot the trees
    #host_tr_plot <- treeducken::phylogram(x$host_tree)
    host_tr_plot <- do.call("phylogram", c(list(tree = x$host_tree), leftArgs))
    left <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
    symb_tr_plot <- do.call("phylogram", c(list(tree = x$symb_tree, direction = "leftwards"), leftArgs))
    #symb_tr_plot <- treeducken::phylogram(x$symb_tree, direction = "leftwards")
    right <- get("last_plot.phylo",envir = ape::.PlotPhyloEnv)
    # trees plotted
    # want to make links for associations
    if (!is.null(x$association_mat))
        make_links(x,
                   c(host_tr_plot, symb_tr_plot),
                   link.type,
                   link.lwd,
                   link.col,
                   link.lty)
    else
        cat("No associations provided.\n")
    if (any(scale.bar > 0))
        add_scalebar(x, scale.bar, sb.fsize)
    assign("last_plot.cophylo", list(left = left,right = right), envir = ape::.PlotPhyloEnv)

}

add_scalebar <- function(obj, scale.bar, fsize) {
    # host scale
    if (scale.bar[1] > 0) {
        s1 <- (0.4 - max(fsize*strwidth(obj$host_tree$tip.label)))/max(ape::node.depth.edgelength(obj$host_tree))
        lines(c(-0.5,-0.5 + scale.bar[1]*s1), rep(-0.05,2))
        lines(rep(-0.5,2), c(-0.05,-0.06))
        lines(rep(-0.5 + scale.bar[1]*s1,2),c(-0.05,-0.06))
        text(mean(c(-0.5,-0.5 + scale.bar[1]*s1)),rep(-0.05,2), scale.bar[1], pos = 1)
    }
    # symbiont scale
    if (scale.bar[2] > 0) {
        s2 <- (0.4 - max(fsize*strwidth(obj$symb_tree$tip.label)))/max(ape::node.depth.edgelength(obj$symb_tree))
        lines(c(0.5 - scale.bar[2]*s2,0.5),rep(-0.05,2))
        lines(rep(0.5 - scale.bar[2]*s2,2),c(-0.05,-0.06))
        lines(rep(0.5,2), c(-0.05,-0.06))
        text(mean(c(0.5 - scale.bar[2]*s2,0.5)), rep(-0.05,2), scale.bar[2], pos = 1)
    }
}

# internal function to plot
#  modified from Liam Revell phytools package
#  under GPL v. 2
phylogram <- function(tree,
                      part = 1,
                      direction = "rightwards",
                      fsize = 1,
                      ftype = "i",
                      lwd = 1,
                      ...){
    if (hasArg(pts))
        pts <- list(...)$pts
    else
        pts <- TRUE
    if (hasArg(edge.col))
        edge.col <- list(...)$edge.col
    else
        edge.col <- rep("black",nrow(tree$edge))
    if (hasArg(tip.lwd))
        tip.lwd <- list(...)$tip.lwd
    else
        tip.lwd <- 1
    if (hasArg(tip.lty))
        tip.lty <- list(...)$tip.lty
    else
        tip.lty <- "dotted"
    if (hasArg(tip.len))
        tip.len <- list(...)$tip.len
    else
        tip.len <- 0.1
    if (pts == TRUE && tip.len == 0)
        tip.len <- 0.1
    d <- ifelse(direction == "rightwards", 1, -1)
    ## sub "_" for " "
    tree$tip.label <- gsub("_"," ", tree$tip.label)
    ## check if edge lenths
    if (ftype == "o")
        fsize <- 0
    num_tips <- length(tree$tip.label) # n
    sh <- fsize*strwidth(tree$tip.label) # sh
    node_depths <- ape::node.depth.edgelength(tree)
    H <- cbind(node_depths[tree$edge[,1]], node_depths[tree$edge[,2]])

    tip_heights <- sapply(1:num_tips, function(i,x,e) x[which(e == i)],
                          x = H[,2],
                          e = tree$edge[,2]) + tip.len*max(H)
    tree$edge.length <- tree$edge.length/max(tip_heights/(part - sh))


    ## reorder cladewise to assign tip positions
    cw_tree <- ape::reorder.phylo(tree,"cladewise")
    y <- vector(length = num_tips + cw_tree$Nnode)
    y[cw_tree$edge[cw_tree$edge[,2] <= num_tips,2]] <- 0:(num_tips - 1)/(num_tips - 1)
    ## reorder pruningwise for post-order traversal
    po_tree <- ape::reorder.phylo(tree,"postorder")

    internal_edges <- unique(po_tree$edge[,1])
    ## compute vertical position of each edge
    for (i in 1:length(internal_edges)) {
        y_hold <- y[po_tree$edge[which(po_tree$edge[,1] == internal_edges[i]),2]]
        y[internal_edges[i]] <- mean(range(y_hold))
    }
    ## compute start & end points of each edge
    X <- ape::node.depth.edgelength(cw_tree) - 0.5
    edge_mat <- cbind(X[cw_tree$edge[,1]], X[cw_tree$edge[,2]])
    ## plot horizontal edges
    for (i in 1:nrow(edge_mat)) {
        lines(d*edge_mat[i,],
              rep(y[cw_tree$edge[i,2]], 2),
              lwd = lwd,
              lend = 2,
              col = edge.col[i])
    }

    ## plot vertical relationships
    for (i in 1:(num_tips + tree$Nnode)) {
        edge_id <- which(cw_tree$edge[,1] == i)
        p <- if (i %in% cw_tree$edge[,2]) which(cw_tree$edge[,2] == i) else NULL

        if (!is.null(p)) {
            x_coords <- c(edge_mat[edge_id,1],edge_mat[p,2])
            y_coords <- sort(c(y[cw_tree$edge[edge_id,2]], y[cw_tree$edge[p,2]]))
        }
        else {
            x_coords <- c(edge_mat[edge_id,1],edge_mat[edge_id[1],1])
            y_coords <- sort(c(y[cw_tree$edge[edge_id,2]],
                               mean(y[cw_tree$edge[edge_id,2]])))
        }
        segments(x0 = d*x_coords[1:(length(x_coords) - 1)],
                 y0 = y_coords[1:(length(y_coords) - 1)],
                 x1 = d*x_coords[2:length(x_coords)],
                 y1 = y_coords[2:length(y_coords)],
                 lwd = lwd,
                 lend = 2,
                 col = edge.col[edge_id])
    }
    h <- part - 0.5 - tip.len*(max(edge_mat) - min(edge_mat)) - fsize*strwidth(tree$tip.label)
    ## plot links to tips
    for (i in 1:num_tips) {
        lines(d*c(edge_mat[which(cw_tree$edge[,2] == i),2],
                  h[i] + tip.len*(max(edge_mat) - min(edge_mat))),
              rep(y[i],2),
              lwd = tip.lwd,
              lty = tip.lty)
        if (pts)
            points(d*edge_mat[which(cw_tree$edge[,2] == i),2], y[i], pch = 16, cex = pts*0.7*sqrt(lwd))
    }
    ## plot tip labels
    font <- which(c("off","reg","b","i","bi") == ftype) - 1
    if (font > 0) {
        for (i in 1:num_tips) {
            make_textbox(d*(h[i] + fsize*strwidth(tree$tip.label[i]) + tip.len*(max(X) - min(X))),
                         y[i],
                         tree$tip.label[i],
                         pos = ifelse(d < 0, 4, 2),
                         offset = 0.1,
                         cex = fsize,
                         font = font)
        }
    }
    phylo_plot <- list(type = "phylogram",
                       use.edge.length = TRUE,
                       node.pos = 1,
                       show.tip.label = ifelse(ftype != "off", TRUE, FALSE),
                       show.node.label = FALSE,
                       font = ftype,
                       cex = fsize,
                       adj = 0,
                       srt = 0,
                       no.margin = FALSE,
                       label.offset = 0,
                       x.lim = par()$usr[1:2],
                       y.lim = par()$usr[3:4],
                       direction = direction,
                       tip.color = "black",
                       Ntip = length(cw_tree$tip.label),
                       Nnode = cw_tree$Nnode,
                       edge = cw_tree$edge,
                       xx = d*sapply(1:(length(cw_tree$tip.label) + cw_tree$Nnode),
                                     function(x,y,z) y[match(x,z)], y = X, z = cw_tree$edge),
                       yy = y)
    assign("last_plot.phylo", phylo_plot, envir = ape::.PlotPhyloEnv)
    ## return rightmost or leftmost edge of tip labels
    invisible(d*max(h + fsize*strwidth(tree$tip.label) + tip.len*(max(X) - min(X))))
}

make_textbox <- function(x, y, label, pos, offset, cex, font){
    rect(x,y - 0.5*strheight(label,cex = cex,font = font),
         x + if (pos == 4) strwidth(label,cex = cex,font = font) else -strwidth(label,cex = cex,font = font),
         y + 0.5*strheight(label,cex = cex,font = font),border = FALSE,
         col = if (par()$bg %in% c("white","transparent")) "white" else par()$bg)
    text(x = x, y = y, label = label, pos = pos, offset = offset, cex = cex, font = font)
}
## plot links between tip taxa according to assoc
## written by Liam J. Revell 2015, 2016, 2019

make_links <- function(obj,
                       x,
                       link.type = "curved",
                       link.lwd = 1,
                       link.col = "red",
                       link.lty = "dashed") {
    host_tree <- geiger::drop.extinct(obj$host_tree)
    symb_tree <- geiger::drop.extinct(obj$symb_tree)
    rownames(obj$association_mat) <- symb_tree$tip.label
    colnames(obj$association_mat) <- host_tree$tip.label
    postorder_ht <- ape::reorder.phylo(obj$host_tree, order = "postorder")
    postorder_st <- ape::reorder.phylo(obj$symb_tree, order = "postorder")
    print(postorder_ht$edge)
    print(postorder_st$edge)

    associations <- which(obj$association_mat == 1, arr.ind = TRUE)
    print(associations)
    # plotting parameters for links between trees
    if (length(link.lwd) == 1)
        link.lwd <- rep(link.lwd, nrow(associations))
    if (length(link.col) == 1)
        link.col <- rep(link.col, nrow(associations))
    if (length(link.lty) == 1)
        link.lty <- rep(link.lty, nrow(associations))
    # end plotting parameteres for links between trees

    # loop through the rows of associations
    # this is a matrix with col 1 giving row of obj$association_mat
    # and col 2 giving col of obj$associatino_mat
    # aka col1 = index of symb_tree$tip.label (so pruned of extinct tips)
    # and col2 = index of host_tree$tip.label (so pruned of extinct tips)
    for (i in 1:nrow(associations)) {
        symbs <- which(postorder_st$edge == associations[i,1])

        hosts <- which(postorder_ht$edge == associations[i,2])
        print(hosts)
        print(symbs)
        y <- c((i - 1)/(length(obj$host_tree$tip.label) - 1),
               (i - 1)/(length(obj$symb_tree$tip.label) - 1))

        if (link.type == "straight")
            lines(x, y, lty = link.lty[i], lwd = link.lwd[i], col = link.col[i])
        else if (link.type == "curved")
            drawCurve(x, y, lty = link.lty[i], lwd = link.lwd[i], col = link.col[i])
    }
}
## function to draw sigmoidal links
## modified from phytools which is
## modified from https://stackoverflow.com/questions/32046889/connecting-two-points-with-curved-lines-s-ish-curve-in-r

drawCurve <- function(x, y, scale=0.01, ...){
    x1 <- x[1]
    x2 <- x[2]
    y1 <- y[1]
    y2 <- y[2]
    curve(plogis(x,
                 scale = scale,
                 location = (x1 + x2)/2)*(y2 - y1) + y1,
          x1,
          x2,
          add = TRUE,
          ...)
}

# this function was taken from Liam Revell's phytools
# package copied under GNU public license 2
#' @describeIn plot.cophylo
#' @param x object of class multiCophylo
#' @export
plot.multiCophylo <- function(x,...){
    par(ask = TRUE)
    for (i in 1:length(x)) {
        plot.cophylo(x[[i]],...)
    }
}
