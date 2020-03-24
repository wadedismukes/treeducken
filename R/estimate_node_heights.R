#' Calculate expected leaves of a species tree
#'
#' @details Calculates the expected number of leaves for a birth-death
#'    simulation given a speciation and extinction rate and a time.
#'
#' @param lambda speciation rate
#' @param mu extinction rate
#' @param t time to simulate until
#' @return The expected number of leaves
#' @references
#' Mooers, A., Gascuel, O., Stadler, T., Li, H., & Steel, M. (2012).
#'     Branch lengths on birth–death trees and the expected loss of phylogenetic
#'      diversity. Systematic biology, 61(2), 195-203.
#' @examples
#' spec_rate <- 1.0
#' ext_rate <- 0.5
#' time <- 2
#' calculate_expected_leaves_sptree(spec_rate, ext_rate, time)
calculate_expected_leaves_sptree <- function(lambda,
                                      mu,
                                      t){
    2 * exp((lambda - mu)*t)
}
#' Calculate expected leaves of a locus tree
#'
#' @details Calculates the expected number of leaves for a birth-death
#'    simulation given a gene birth and death rate, a time, and the number of
#'    leaves on the species tree that the locus tree is to be simulated upon.
#'
#' @param t time to simulate until (the length of the species tree)
#' @param dup_rate gene birth rate
#' @param loss_rate gene death rate
#' @param num_species number of leaves on the species tree
#' @return The expected number of leaves
#' @references
#' Mallo, D., de Oliveira Martins, L., & Posada, D. (2016). SimPhy: phylogenomic
#'      simulation of gene, locus, and species trees. Systematic biology, 65(2),
#'       334-344.

#' @examples
#' gene_birth_rate <- 1.0
#' gene_death_rate <- 0.5
#' time <- 2
#' num_species <- 10
#' calculate_expected_leaves_locustree(time,
#'                                     gene_birth_rate,
#'                                     gene_death_rate,
#'                                     num_species)
calculate_expected_leaves_locustree <- function(t,
                                                dup_rate,
                                                loss_rate,
                                                num_species){
    f_numer <- (loss_rate * (1 - exp(-(dup_rate - loss_rate) * t)))
    f_denom <- dup_rate - loss_rate * (exp(-(dup_rate - loss_rate) * t))
    f <- f_numer / f_denom
    (num_species * exp(t * (dup_rate - loss_rate))) / (1 - f^(num_species - 2))
}
#' Calculate expected time to branching point of a species tree
#'
#' @details Calculates the expected time to branching point of a species tree
#'     for a birth-death simulation given a speciation and extinction rate and a
#'      a number of leaves, and a branching point. By default this branching
#'      point is 1 which corresponds to the root, and thus this function will
#'      return the depth of the tree.
#'
#'
#'
#' @param lambda speciation rate
#' @param mu extinction rate
#' @param t time to simulate until
#' @return The expected number of leaves
#' @references
#' Gernhard, T. (2008). The conditioned reconstructed process. Journal of
#'     theoretical biology, 253(4), 769-778.
#' @examples
#' spec_rate <- 1.0
#' ext_rate <- 0.5
#' n <- 10
#' estimate_node_heights(spec_rate, ext_rate, n)
#'
#' estimate_node_heigths(spec_rate, ext_rate, n, k = 5)
estimate_node_heights <- function(lambda,
                                  mu,
                                  k = 1,
                                  n){
    if(lambda <= 0.0)
        stop("'lambda' is lower than 0. it must be greater than 0.")
    if(mu > lambda)
        stop("'mu' is greater than 'lambda', please correc this.")
    if(mu == 0.0){
        s <- 0
        for(i in k+1:n){
            s <- s + 1/(lambda*i)
        }
        speciation_event_time <- s
    }
    else if(mu == lambda){
        speciation_event_time <- (n-k)/(lambda*k)
    }
    else{
        rho <- mu / lambda
        first_part <- ((k + 1) / lambda) * exp(lchoose(n, k + 1)) * (-1)^k
        out_sum <- 0
        for(i in 0:(n-k-1)){
            s1 <- exp(lchoose(n - k - 1, i))
            s2 <- 1 / ((k + i + 1) * rho)
            s3 <- exp((k+i) * log(1/rho -1))    # (1/rho - 1)^(k+i)

            in_sum <- 0
            for(j in 1:k+i){
                ins1 <- exp(lchoose(k+i, j))
                ins2 <- (-1^j)/j
                ins3 <- 1 - (1/(1-rho))^j
                print(ins1)
                print(ins2)
                print(ins3)
                in_sum <- in_sum + ins1*ins2*ins3
            }
            s4 <- log(1/(1-rho)) - in_sum
            print(s4)
            out_sum <- out_sum + s1*s2*s3*s4
        }
        speciation_event_time <- first_part * out_sum
    }
    speciation_event_time
}


