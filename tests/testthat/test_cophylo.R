


get_length_host_tree <- function(cophy){
    max(apTreeshape:::nodes_depths_ordonnes(cophy$host_tree)) + cophy$host_tree$root.edge
}

get_all_host_tree_lengths <- function(multiCoph){
    unique(sapply(multiCoph, get_length_host_tree))
}

# test that tree has correct length (simple)
test_that("sim_cophylo_bdp produces the right length trees", {
    expect_equal(get_all_host_tree_lengths(sim_cophylo_bdp(hbr = 0.5,
                                        hdr = 0.3,
                                        sbr = 1.0,
                                        sdr = 0.3,
                                        host_exp_rate = 0.0,
                                        cosp_rate = 0.5,
                                        timeToSimTo = 2.0,
                                        numbsim = 10)), 2.0)
    expect_equal(get_all_host_tree_lengths(sim_cophylo_bdp(hbr = 0.5,
                                        hdr = 0.3,
                                        sbr = 1.0,
                                        sdr = 0.3,
                                        host_exp_rate = 0.0,
                                        cosp_rate = 0.5,
                                        timeToSimTo = 1.5,
                                        numbsim = 10)), 1.5)
    expect_equal(get_all_host_tree_lengths(sim_cophylo_bdp(hbr = 0.5,
                                        hdr = 0.3,
                                        sbr = 1.0,
                                        sdr = 0.3,
                                        host_exp_rate = 0.0,
                                        cosp_rate = 0.5,
                                        timeToSimTo = 1.0,
                                        numbsim = 10)), 1.0)
})

# test that tree has correct extant tips
test_that("sim_cophylo_bdp produces the right number of trees", {
    expect_equal(length(sim_cophylo_bdp(hbr = 0.5,
                                        hdr = 0.3,
                                        sbr = 1.0,
                                        sdr = 0.3,
                                        host_exp_rate = 0.0,
                                        cosp_rate = 0.5,
                                        timeToSimTo = 1.5,
                                        numbsim = 10)), 10)
    expect_equal(length(sim_cophylo_bdp(hbr = 0.5,
                                        hdr = 0.3,
                                        sbr = 1.0,
                                        sdr = 0.3,
                                        host_exp_rate = 0.0,
                                        cosp_rate = 0.5,
                                        timeToSimTo = 1.5,
                                        numbsim = 5)), 5)
    expect_equal(length(sim_cophylo_bdp(hbr = 0.5,
                                        hdr = 0.3,
                                        sbr = 1.0,
                                        sdr = 0.3,
                                        host_exp_rate = 0.0,
                                        cosp_rate = 0.5,
                                        timeToSimTo = 1.5,
                                        numbsim = 20)), 20)
})

is_host_and_symbiont_the_same <- function(t, n){
    pair <- sim_cophylo_bdp(hbr = 0.0,
                            hdr = 0.0,
                            sbr = 0.0,
                            sdr = 0.0,
                            host_exp_rate = 0.0,
                            cosp_rate = 1.0,
                            timeToSimTo = t,
                            numbsim = n)
    pair_true <- vector(length = n)
    for(i in 1:n){
        pair_true <- is.equal.phylo(pair[[i]]$host_tree, pair[[i]]$symb_tree)
    }
    all(pair_true == TRUE)
}

test_that("sim_cophylo_bdp produces identical symbiont and host trees when only cospeciation is present", {
    expect_true(is_host_and_symbiont_the_same(t = 1.5, n = 10))
})