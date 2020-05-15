


get_length_host_tree <- function(cophy){
    max(phytools::nodeHeights(cophy$host_tree)) + cophy$host_tree$root.edge
}

get_all_host_tree_lengths <- function(multiCoph){
    mean(sapply(multiCoph, get_length_host_tree))
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
