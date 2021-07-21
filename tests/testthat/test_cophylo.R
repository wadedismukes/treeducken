get_number_hosts <- function(cophylo_obj, hl) {
    assoc_mats <- association_mat.multiCophy(cophylo_obj)
    x <- lapply(assoc_mats, colSums)
    y <- lapply(x, function(z) {any(z > hl)})
    any(unlist(y))
}


test_that("host_limit is set correctly", {
    expect_equal(get_number_hosts(sim_cophyBD(hbr = 0.5,
                                                  hdr = 0.3,
                                                  sbr = 1.0,
                                                  sdr = 0.15,
                                                  host_exp_rate = 0.15,
                                                  cosp_rate = 0.5,
                                                  time_to_sim = 2.0,
                                                  numbsim = 10,
                                                  host_limit = 2), 2), FALSE)
    expect_equal(get_number_hosts(sim_cophyBD(hbr = 0.5,
                                                  hdr = 0.3,
                                                  sbr = 1.0,
                                                  sdr = 0.15,
                                                  host_exp_rate = 0.15,
                                                  cosp_rate = 0.5,
                                                  time_to_sim = 2.0,
                                                  numbsim = 10,
                                                  host_limit = 3), 3), FALSE)
    expect_equal(get_number_hosts(sim_cophyBD(hbr = 0.5,
                                                  hdr = 0.3,
                                                  sbr = 1.0,
                                                  sdr = 0.15,
                                                  host_exp_rate = 0.15,
                                                  cosp_rate = 0.5,
                                                  time_to_sim = 2.0,
                                                  numbsim = 10,
                                                  host_limit = 4), 4), FALSE)
})

test_that("host_limit is set correctly for anagenetic", {
    expect_equal(get_number_hosts(sim_cophyBD_ana(hbr = 0.5,
                                                  hdr = 0.3,
                                                  sbr = 1.0,
                                                  sdr = 0.15,
                                                  host_exp_rate = 0.15,
                                                  s_disp_r = 0.2,
                                                  s_extp_r = 0.05,
                                                  cosp_rate = 0.5,
                                                  time_to_sim = 2.0,
                                                  numbsim = 10,
                                                  host_limit = 2), 2), FALSE)
    expect_equal(get_number_hosts(sim_cophyBD_ana(hbr = 0.5,
                                                  hdr = 0.3,
                                                  sbr = 1.0,
                                                  sdr = 0.15,
                                                  s_disp_r = 0.2,
                                                  s_extp_r = 0.05,
                                                  host_exp_rate = 0.15,
                                                  cosp_rate = 0.5,
                                                  time_to_sim = 2.0,
                                                  numbsim = 10,
                                                  host_limit = 3), 3), FALSE)
    expect_equal(get_number_hosts(sim_cophyBD_ana(hbr = 0.5,
                                                  hdr = 0.3,
                                                  sbr = 1.0,
                                                  sdr = 0.15,
                                                  s_disp_r = 0.2,
                                                  s_extp_r = 0.05,
                                                  host_exp_rate = 0.15,
                                                  cosp_rate = 0.5,
                                                  time_to_sim = 2.0,
                                                  numbsim = 10,
                                                  host_limit = 4), 4), FALSE)
})


get_length_host_tree <- function(cophy) {
    max(ape::node.depth.edgelength(cophy$host_tree)) + cophy$host_tree$root.edge
}

get_all_host_tree_lengths <- function(multiCoph) {
    min(sapply(multiCoph, get_length_host_tree))
}

# test that tree has correct length (simple)
test_that("sim_cophyBD produces the right length trees", {
    expect_equal(get_all_host_tree_lengths(sim_cophyBD(hbr = 0.5,
                                        hdr = 0.3,
                                        sbr = 1.0,
                                        sdr = 0.3,
                                        host_exp_rate = 0.0,
                                        cosp_rate = 0.5,
                                        time_to_sim = 2.0,
                                        numbsim = 10)), 2.0)
    expect_equal(get_all_host_tree_lengths(sim_cophyBD(hbr = 0.5,
                                        hdr = 0.3,
                                        sbr = 1.0,
                                        sdr = 0.3,
                                        host_exp_rate = 0.0,
                                        cosp_rate = 0.5,
                                        time_to_sim = 1.5,
                                        numbsim = 10)), 1.5)
    expect_equal(get_all_host_tree_lengths(sim_cophyBD(hbr = 0.5,
                                        hdr = 0.3,
                                        sbr = 1.0,
                                        sdr = 0.3,
                                        host_exp_rate = 0.0,
                                        cosp_rate = 0.5,
                                        time_to_sim = 1.0,
                                        numbsim = 10)), 1.0)
})

# test that tree has correct extant tips
test_that("sim_cophyBD produces the right number of trees", {
    expect_equal(length(sim_cophyBD(hbr = 0.5,
                                        hdr = 0.3,
                                        sbr = 1.0,
                                        sdr = 0.3,
                                        host_exp_rate = 0.0,
                                        cosp_rate = 0.5,
                                        time_to_sim = 1.5,
                                        numbsim = 10)), 10)
    expect_equal(length(sim_cophyBD(hbr = 0.5,
                                        hdr = 0.3,
                                        sbr = 1.0,
                                        sdr = 0.3,
                                        host_exp_rate = 0.0,
                                        cosp_rate = 0.5,
                                        time_to_sim = 1.5,
                                        numbsim = 5)), 5)
    expect_equal(length(sim_cophyBD(hbr = 0.5,
                                        hdr = 0.3,
                                        sbr = 1.0,
                                        sdr = 0.3,
                                        host_exp_rate = 0.0,
                                        cosp_rate = 0.5,
                                        time_to_sim = 1.5,
                                        numbsim = 20)), 20)
})

is_host_and_symbiont_the_same <- function(t, n) {
    pair <- sim_cophyBD(hbr = 0.0,
                            hdr = 0.0,
                            sbr = 0.0,
                            sdr = 0.0,
                            host_exp_rate = 0.0,
                            cosp_rate = 1.0,
                            time_to_sim = t,
                            numbsim = n)
    pair_true <- vector(length = n)
    for(i in 1:n){
        pair_true[i] <- ape::all.equal.phylo(pair[[i]]$host_tree,
                                     pair[[i]]$symb_tree,
                                     use.tip.label = FALSE)
    }
    all(pair_true == TRUE)
}

test_that("sim_cophyBD produces identical symbiont and host trees when only cospeciation is present", {
    expect_true(is_host_and_symbiont_the_same(t = 1.5, n = 10))
})

any_host_rows_zero <- function(t, n) {
    pair <- sim_cophyBD(hbr = 1.0,
                        hdr = 0.75,
                        sbr = 1.0,
                        sdr = 0.75,
                        host_exp_rate = 1.0,
                        cosp_rate = 1.0,
                        time_to_sim = t,
                        numbsim = n,
                        mutualism = TRUE)
    mats <- association_mat(pair)
    y <- lapply(mats, FUN = function(x) {
        all(rowSums(x) != 0)
    })
    all(unlist(y))
}

test_that("sim_cophyBD mutualism produces no association matrices with empty rows", {
    expect_true(any_host_rows_zero(t = 1.0, n = 10))
})
