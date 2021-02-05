library(treeducken)


# test species tree output is a list of trees with correct length
test_that("sim_st_bdp produces the right number of trees", {
    expect_equal(length(sim_st_bdp(1.0, 0.5, 10, 10)), 10)
    expect_equal(length(sim_st_bdp(1.0, 0.5, 5, 10)), 5)
    expect_equal(length(sim_st_bdp(1.0, 0.5, 20, 10)), 20)
})

get_number_extant_tips <- function(tr) {
    tip_vec <- vector(length = length(tr))
    for(i in 1:length(tr)) {
        pruned_tr <- treeducken::drop_extinct(tr[[i]], tol = 0.0001)
        tip_vec[i] <- length(pruned_tr$tip.label)
    }
    min(tip_vec)
}

# test that tree has correct extant tips (gsa)
test_that("sim_st_bdp produces the right number of extant tips", {
    expect_equal(get_number_extant_tips(sim_st_bdp(1.0, 0.5, 100, n_tips = 10)), 10)
    expect_equal(get_number_extant_tips(sim_st_bdp(1.0, 0.5, 100, n_tips = 5)), 5)
    expect_equal(get_number_extant_tips(sim_st_bdp(1.0, 0.5, 100, n_tips = 20)), 20)
})


# test that tree has correct extant tips
test_that("sim_st_bdp_t produces the right number of trees", {
    expect_equal(length(sim_st_bdp_t(1.0, 0.5, 10, 2.0)), 10)
    expect_equal(length(sim_st_bdp_t(1.0, 0.5, 5, 2.0)), 5)
    expect_equal(length(sim_st_bdp_t(1.0, 0.5, 20, 2.0)), 20)
})


get_length_tree <- function(tr){
    max(ape::node.depth.edgelength(tr)) + tr$root.edge
}

get_all_tree_lengths <- function(multiTree){
    min(sapply(multiTree, get_length_tree))
}

# test that tree has correct length (simple)
test_that("sim_st_bdp_t produces the right length trees", {
    expect_equal(get_all_tree_lengths(sim_st_bdp_t(1.0, 0.5, 1000, 1.0)), 1.0)
    expect_equal(get_all_tree_lengths(sim_st_bdp_t(1.0, 0.5, 1000, 2.0)), 2.0)
    expect_equal(get_all_tree_lengths(sim_st_bdp_t(1.0, 0.5, 1000, 3.0)), 3.0)
})
