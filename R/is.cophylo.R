#' Test for the cophylogenetic set object
#' @description Tests if an object is of class `cophylo`
#'
#' @param cophy an object to test to see if it is of class `cophylo`
#' @details Checks that an object is of class `cophylo`. For multicophy checks that
#' the class is `multiCophylo` and that each element is of class `cophylo`.
#' @return A logical vector
#' @seealso `as.cophylo`
#' @examples
#' h_lambda <- 1.0
#' h_mu <- 0.3
#' c_lambda <- 0.0
#' s_lambda <- 1.0
#' s_mu <- 0.3
#' s_her <- 0.0
#' host_symb_sets <- sim_cophylo_bdp(hbr = h_lambda,
#'                                   hdr = h_mu,
#'                                   sbr = s_lambda,
#'                                   cosp_rate = c_lambda,
#'                                   sdr = s_mu,
#'                                   host_exp_rate = s_her,
#'                                   timeToSimTo = 2.0,
#'                                   numbsim = 1)
#' is.cophylo(host_symb_sets[[1]])
#' is.multiCophylo(host_symb_sets)
is.cophylo <- function(cophy){
    inherits(cophy, "cophylo")
}
#' @describeIn is.cophylo Tests for `multiCophylo` composed of `cophylo` objects
#' @param multiCophylo an object to test for multiCophy
is.multiCophylo <- function(multiCophy){
    t <- inherits(multiCophy, "multiCophylo")
    if(t){
        tt <- sapply(unclass(cophy), inherits, what = "cophylo")
    }
    all(c(t, tt))
}