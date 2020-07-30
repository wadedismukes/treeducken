## cophy_plot.R (2014-04-07)

##   Plots two phylogenetic trees face to
##   face with the links between the tips

## Copyright 2008-2010 Damien de Vienne

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

plot.cophy <-
    function(cophy_obj,
             use_edge_length = TRUE, space = 1,
             length_line = 1, gap = 0, type = "phylogram",
             col = par("fg"), lwd = par("lwd"), lty = par("lty"),
             show_tip_label = TRUE, font = 3, fsize = 1.0,
            ...) {

    if (!inherits(cophy_obj, "cophy"))
        stop("cophy_obj should be an object of class 'cophy'.")


    host <- cophy_obj$host_tree
    symb <- cophy_obj$symb_tree
    host_tree_pruned <- geiger::drop.extinct(host)
    symb_tree_pruned <- geiger::drop.extinct(symb)
    rownames(cophy_obj$association_mat) <- symb_tree_pruned$tip.label
    colnames(cophy_obj$association_mat) <- host_tree_pruned$tip.label
    assoc <- which(cophy_obj$association_mat == 1, arr.ind = TRUE)
    assoc <- cbind(host_tree_pruned$tip.label[assoc[, 2]],
                   symb_tree_pruned$tip.label[assoc[, 1]])

    draw_cophy(host,
               symb,
               assoc = assoc,
               use_edge_length = use_edge_length,
               space = space,
               length_line = length_line,
               gap = gap,
               type = type,
               return = FALSE,
               col = col,
               lwd = lwd,
               lty = lty,
               show_tip_label = show_tip_label,
               font = font,
               fsize = fsize)
}


# x is the host
# y is the symb
draw_cophy <-
    function(x,
             y,
             assoc = assoc,
             use_edge_length = use_edge_length,
             space = space, 
             length_line = length_line,
             gap = gap,
             type = type,
             return = return,
             col = col,
             lwd=lwd,
             lty=lty,
             show_tip_label = show_tip_label,
             font = font,
             fsize = fsize,
             ...) {
    res <- list()

###choice of the minimum space between the trees
    left <- min(nchar(x$tip.label, type = "width")) * fsize
    right <- min(nchar(y$tip.label, type = "width")) * fsize
    space_min <- left + right + gap * 2
    if ((space <= 0) || (space < space_min)) space <- space_min
    n_tip_x <- ape::Ntip(x)
    n_tip_y <- ape::Ntip(y)
    res$n_tip_x <- n_tip_x
    res$n_tip_y <- n_tip_y
    # a is coordinates of host
    a <- ape::plotPhyloCoor(x,
                            use_edge_length = use_edge_length,
                            type = type)
    res$a <- a
    # b is coordinates of symbiont
    b <- ape::plotPhyloCoor(y,
                            use_edge_length = use_edge_length,
                            direction = "leftwards",
                            type = type)
###for the two trees to have the extreme leaves at the same ordinate.
    a[, 2] <- a[, 2] - min(a[, 2])
    b[, 2] <- b[, 2] - min(b[, 2])
    res$b <- b
    b2 <- b
    b2[, 1] <- b[seq_len(nrow(b)), 1] * (max(a[, 1]) / max(b[, 1])) +
        space + max(a[, 1])
    b2[, 2] <- b[seq_len(nrow(b)), 2] * (max(a[, 2]) / max(b[, 2]))
    res$b2 <- b2
    c <- matrix(ncol = 2, nrow = nrow(a) + nrow(b))
    c[seq_len(nrow(a)), ] <- a[seq_len(nrow(a)), ]
    c[nrow(a) + seq_len(nrow(b)), 1] <- b2[, 1]
    c[nrow(a) + seq_len(nrow(b)), 2] <- b2[, 2]
    res$c <- c
    plot(c, xlim = NULL, ylim = NULL, log = "", main = NULL,
        sub = NULL, xlab = NULL, ylab = NULL, ann = FALSE, axes = FALSE,
        frame.plot = FALSE)
 ###segments for cladograms
   if (type == "cladogram") {
        for (i in 1:(nrow(a) - 1)) segments(a[x$edge[i, 1], 1],
            a[x$edge[i, 1], 2], a[x$edge[i, 2], 1], a[x$edge[i,
                2], 2], col = "red")
        for (i in 1:(nrow(b) - 1))
            segments(b2[y$edge[i, 1], 1], b2[y$edge[i, 1], 2],
                     b2[y$edge[i, 2], 1], b2[y$edge[i, 2], 2])
    }
###segments for phylograms
    if (type == "phylogram") {
        for (i in (n_tip_x + 1):nrow(a)) {
            l <- length(x$edge[x$edge[, 1] == i, ][, 1])
            for (j in 1:l) {
                segments(a[x$edge[x$edge[, 1] == i, ][1, 1],
                  1], a[x$edge[x$edge[, 1] == i, 2], 2][1], a[x$edge[x$edge[,
                  1] == i, ][1, 1], 1], a[x$edge[x$edge[, 1] ==
                  i, 2], 2][j])
                segments(a[x$edge[x$edge[, 1] == i, ][1, 1], 1],
                         a[x$edge[x$edge[, 1] == i, 2], 2][j],
                         a[x$edge[x$edge[, 1] == i, 2], 1][j],
                         a[x$edge[x$edge[, 1] == i, 2], 2][j])
            }
        }
        for (i in (n_tip_y + 1):nrow(b)) {
            l <- length(y$edge[y$edge[, 1] == i, ][, 1])
            for (j in 1:l) {
                segments(b2[y$edge[y$edge[, 1] == i, ][1, 1], 1],
                         b2[y$edge[y$edge[, 1] == i, 2], 2][1],
                         b2[y$edge[y$edge[, 1] == i, ][1, 1], 1],
                         b2[y$edge[y$edge[, 1] == i, 2], 2][j])
                segments(b2[y$edge[y$edge[, 1] == i, ][1, 1], 1],
                         b2[y$edge[y$edge[, 1] == i, 2], 2][j],
                         b2[y$edge[y$edge[, 1] == i, 2], 1][j],
                         b2[y$edge[y$edge[, 1] == i, 2], 2][j])
            }
        }
    }
    if (show_tip_label) {
        # add host tips
        make_textbox(a[1:n_tip_x, 1], a[1:n_tip_x, 2],
                    label = x$tip.label,
                    pos = 4,
                    offset = 0.1,
                    cex = fsize,
                    font = font)
        # symb tips
        make_textbox(b2[1:n_tip_y, 1], b2[1:n_tip_y, 2],
                    label = y$tip.label,
                    pos = 2,
                    offset = 0.1,
                    cex = fsize,
                    font = font)


    }
###
## Plot links between associated taxa.
## Takes into account the size of the character strings of the taxa names.


    lsa <- 1:n_tip_x
    lsb <- 1:n_tip_y
    decx <- array(nrow(assoc))
    decy <- array(nrow(assoc))


    #colors
    if (length(col) == 1) colors <- c(rep(col, nrow(assoc)))
    else if (length(col) >= nrow(assoc)) colors <- col
    else  colors <- c(rep(col, as.integer(nrow(assoc) / length(col)) + 1))

    #lwd
    if (length(lwd) == 1) lwidths <- c(rep(lwd, nrow(assoc)))
    else if (length(lwd) >= nrow(assoc)) lwidths <- lwd
    else  lwidths <- c(rep(lwd, as.integer(nrow(assoc) / length(lwd)) + 1))

    #lty
    if (length(lty) == 1) ltype <- c(rep(lty, nrow(assoc)))
    else if (length(lty) >= nrow(assoc)) ltype <- lty
    else  ltype <- c(rep(lty, as.integer(nrow(assoc) / length(lty)) + 1))


    for (i in seq_len(nrow(assoc))) {
        if (show_tip_label) {
            decx[i] <- strwidth(x$tip.label[lsa[x$tip.label == assoc[i, 1]]]) * fsize + 1 * fsize
            decy[i] <- strwidth(y$tip.label[lsb[y$tip.label == assoc[i, 2]]]) * fsize + 1 * fsize
        } else {
            decx[i] <- decy[i] <- 0
        }
        segments(a[lsa[x$tip.label == assoc[i, 1]], 1] + decx[i],
                 a[lsa[x$tip.label == assoc[i, 1]], 2],
                 b2[lsb[y$tip.label == assoc[i, 2]], 1] - decy[i],
                 b2[lsb[y$tip.label == assoc[i, 2]], 2],
                 col = colors[i], lwd = lwidths[i], lty = ltype[i])
    }
    if (return == TRUE)  res
}



#' Internal tree plot function
#' @description internal function to make textbox for tip labels modified from Liam Revell phytools package under GPL v. 2
#'
#' @param x x coordinates
#' @param y y coordinates
#' @param label Labels as vector of strings
#' @param pos Position in plot environment
#' @param offset How offset from tips
#' @param cex a numerical vector giving the amount by which plotting characters and symbols should be scaled relative to the default
#' @param font font choice
#'
#'
#' @references
#' Revell, L.J. (2012), phytools: an R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution, 3: 217-223. doi:10.1111/j.2041-210X.2011.00169.x
#' @keywords Internal
make_textbox <- function(x, y, label, pos, offset, cex, font) {
    rect(x, y - 0.5 * strheight(label, cex = cex, font = font),
         x + if (pos == 4) strwidth(label, cex = cex, font = font) 
             else -strwidth(label, cex = cex, font = font),
         y + 0.5 * strheight(label, cex = cex, font = font), border = FALSE,
         col = if (par()$bg %in% c("white", "transparent")) "white"
               else par()$bg)
    text(x = x,
         y = y,
         label = label,
         pos = pos,
         offset = offset,
         cex = cex,
         font = font)
}


## function to draw sigmoidal links
## modified from phytools which is
## modified from https://stackoverflow.com/questions/32046889/connecting-two-points-with-curved-lines-s-ish-curve-in-r
## plot links between tip taxa according to assoc
## written by Liam J. Revell 2015, 2016, 2019
#' Curve draw function
#' @description internal function to draw curved links between tips modified from Liam Revell phytools package under GPL v. 2
#'
#' @param x x positions on graph
#' @param y y positions on graph
#' @param scale Scale of the logistic (which is where the curve comes from)
#' @param ... Other plotting parameters
#'
#' @references
#' Revell, L.J. (2012), phytools: an R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution, 3: 217-223. doi:10.1111/j.2041-210X.2011.00169.x
#' @keywords Internal
draw_curve <- function(x, y, scale=0.01, ...){
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
#' @describeIn plot.cophy Plots multiple cophy plots
#' @param x object of class multiCophy
#' @export
plot.multiCophy <- function(x,...){
    par(ask = TRUE)
    for (i in 1:length(x)) {
        plot.cophy(x[[i]],...)
    }
}
