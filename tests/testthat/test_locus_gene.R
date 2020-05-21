# test species tree output is a list of trees with correct length
sim_test_spt_loct <- function(gbr, gdr, lgtr, numLoci, species_tree_len = 5.0){
    tr <- sim_sptree_bdp_time(0.1, 0.05, 1, species_tree_len)
    sim_locustree_bdp(tr[[1]], gbr = 0.0, gdr = 0.0, lgtr = 0.0, num_loci = numLoci)
}

sim_test_spt_loct_equality <- function(gene_birth = 0.0, gene_death = 0.0, transfers = 0.0, numLoci){
    tr <- sim_sptree_bdp_time(0.1, 0.05, 1, 5.0)
    loctr <- sim_locustree_bdp(tr[[1]], gbr = gene_birth, gdr = gene_death, lgtr = transfers, num_loci = numLoci)
    ape::all.equal.phylo(tr[[1]], loctr)
}

test_that("sim_locustree_bdp produces the right number of trees", {
    expect_equal(length(sim_test_spt_loct(1.0, 0.5, 0.5, 10)), 10)
    expect_equal(length(sim_test_spt_loct(1.0, 0.5, 0.5, 5)), 5)
    expect_equal(length(sim_test_spt_loct(1.0, 0.5, 0.5, 50)), 50)
})



# test that input species tree is the same as locus tree if every param is set
# 0.0
test_that("sim_locustree_bdp returns tree concordant with species tree when locus tree parameters are 0.0",{
    expect_true(sim_test_spt_loct_equality(gbr = 0.0, gdr = 0.0, lgtr = 0.0, numLoci = 20))
})

# test that length is greater than or equal to input species tree
test_that("sim_locustree_bdp returns a tree greather than or equal in length to species tree", {
    expect_gte(min(sim_test_spt_loct(1.0,
                                        0.5,
                                        0.5,
                                        10,
                                        3.0)), 3.0)
    expect_gte(min(sim_test_spt_loct(1.0,
                                     0.5,
                                     0.5,
                                     10,
                                     5.0)), 5.0)
    expect_gte(min(sim_test_spt_loct(1.0,
                                     0.5,
                                     0.5,
                                     10,
                                     4.0)), 4.0)
})