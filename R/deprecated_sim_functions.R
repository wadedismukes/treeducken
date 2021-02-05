#' @export
#' @rdname sim_st_bdp
sim_sptree_bdp<- function(sbr,
                          sdr,
                          numbsim,
                          n_tips,
                          gsa_stop_mult = 10) {
    warning("please use sim_st_bdp instead of sim_sptree_bdp", call. = FALSE)
    sim_st_bdp(sbr, sdr, numbsim, n_tips, gsa_stop_mult)
}
#' @export
#' @rdname sim_st_bdp_t
sim_sptree_bdp_time <- function(sbr,
                            sdr,
                            numbsim,
                            t) {
    warning("please use sim_st_bdp_t instead of sim_sptree_bdp_time", call. = FALSE)
    sim_st_bdp_t(sbr, sdr, numbsim, t)
}

#' @export
#' @rdname sim_lt_bdp
sim_locustree_bdp <- function(species_tree,
                             gbr,
                             gdr,
                             lgtr,
                             num_loci,
                             transfer_type = "random") {
    warning("please use sim_lt_bdp instead of sim_locustree_bdp", call. = FALSE)
    sim_lt_bdp(species_tree, gbr, gdr, lgtr, num_loci, transfer_type)
}
#' @export
#' @rdname sim_cbdp_ana
sim_cophylo_bdp_ana <- function(hbr,
                                hdr,
                                sbr,
                                sdr,
                                symb_dispersal_rate,
                                symb_extirpation_rate,
                                host_exp_rate,
                                cosp_rate,
                                time_to_sim,
                                numbsim,
                                host_limit = 0) {
    warning("please use sim_cbdp_ana instead of sim_cophylo_bdp_ana", call. = FALSE)
    sim_cbdp_ana(hbr,
                 hdr,
                 sbr,
                 sdr,
                 symb_dispersal_rate,
                 symb_extirpation_rate,
                 host_exp_rate,
                 cosp_rate,
                 time_to_sim,
                 numbsim,
                 host_limit)
}
#' @export
#' @rdname sim_cbdp
sim_cophylo_bdp <- function(hbr,
                            hdr,
                            sbr,
                            sdr,
                            host_exp_rate,
                            cosp_rate,
                            time_to_sim,
                            numbsim,
                            host_limit = 0) {
    warning("please use sim_cbdp instead of sim_cophylo_bdp", call. = FALSE)
    sim_cbdp(hbr,
                 hdr,
                 sbr,
                 sdr,
                 host_exp_rate,
                 cosp_rate,
                 time_to_sim,
                 numbsim,
                 host_limit)
}
#' @export
#' @rdname sim_msc
sim_multispecies_coal <- function(species_tree,
                                  ne,
                                  num_sampled_individuals,
                                  num_genes,
                                  rescale = TRUE,
                                  mutation_rate = 1) {

    warning("please use sim_msc() instead of sim_multispecies_coal()", call.=FALSE)
    sim_msc(species_tree,
            ne,
            num_sampled_individuals,
            num_genes,
            rescale,
            mutation_rate)
}