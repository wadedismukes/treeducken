library(treeducken)


# test species tree output is a list of trees with correct length
test_that("sim_sptree_bdp produces the right number of trees", {
    expect_equal(length(sim_sptree_bdp(1.0, 0.5, 10, 10)), 10)
    expect_equal(length(sim_sptree_bdp(1.0, 0.5, 5, 10)), 5)
    expect_equal(length(sim_sptree_bdp(1.0, 0.5, 20, 10)), 20)
})

get_number_extant_tips <- function(tr){
    tip_vec <- vector(length = length(tr))
    for(i in 1:length(tr)){
        pruned_tr <- geiger::drop.extinct(tr[[i]], tol = 0.001)
        tip_vec[i] <- length(pruned_tr$tip.label)
    }
    mean(tip_vec)
}

# test that tree has correct extant tips (gsa)
test_that("sim_sptree_bdp produces the right number of extant tips", {
    expect_equal(get_number_extant_tips(sim_sptree_bdp(1.0, 0.5, 10, 10)), 10)
      expect_equal(get_number_extant_tips(sim_sptree_bdp(1.0, 0.5, 10, 5)), 5)
    expect_equal(get_number_extant_tips(sim_sptree_bdp(1.0, 0.5, 10, 20)), 20)
})

# test that species tree produces tree within correct distribution (gsa)
get_treesim_treedepth_dist <- function(sbr, sdr, nt, reps){
  trees <- TreeSim::sim.bd.taxa(lambda = sbr, mu = sdr, n = nt, numbsim = reps)
  max(phytools::nodeHeights(trees))
}


# test that tree has correct extant tips
test_that("sim_sptree_bdp_time produces the right number of trees", {
    expect_equal(length(sim_sptree_bdp_time(1.0, 0.5, 10, 2.0)), 10)
    expect_equal(length(sim_sptree_bdp_time(1.0, 0.5, 5, 2.0)), 5)
    expect_equal(length(sim_sptree_bdp_time(1.0, 0.5, 20, 2.0)), 20)
})


get_length_tree <- function(tr){
    tree_depth <- vector(length = length(tr))
    for(i in 1:length(tr)){
        tree_depth[i] <- max(phytools::nodeHeights(tr[[i]])) + tr[[i]]$root.edge
    }
    mean(tree_depth)
}
# test that tree has correct length (simple)
test_that("sim_sptree_bdp_time produces the right length trees", {
    expect_equal(get_length_tree(sim_sptree_bdp_time(1.0, 0.5, 10, 1.0)), 1.0)
    expect_equal(get_length_tree(sim_sptree_bdp_time(1.0, 0.5, 10, 2.0)), 2.0)
    expect_equal(get_length_tree(sim_sptree_bdp_time(1.0, 0.5, 10, 5.0)), 5.0)
})

# test that species tree produces tree within correct distribution (simple)
# test_that("sim_sptree_bdp_time produces trees under the right distribution"){
#
# }
