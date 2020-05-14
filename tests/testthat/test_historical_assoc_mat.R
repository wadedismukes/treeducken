library(treeducken)

hist_assoc_tester <- function(reps = 1, time){
    host_mu <- 0.0 # death rate
    host_lambda <- 0.0 # birth rate
    numb_replicates <- reps
    time <- time
    symb_mu <- 0.0
    symb_lambda <- 0.0
    host_shift_rate <- 0.0
    cosp_rate <- 2.0

    cophylo_pair <- sim_cophylo_bdp(hbr = host_lambda,
                                    hdr = host_mu,
                                    cosp_rate = cosp_rate,
                                    host_exp_rate = host_shift_rate,
                                    sdr = symb_mu,
                                    sbr = symb_lambda,
                                    numbsim = numb_replicates,
                                    timeToSimTo = time)
    build_historical_association_matrix(t = time/2, cophylo_pair[[1]])
}

test_that("build historical association matrix gives I if only cospeciation is on", {
    x <- hist_assoc_tester(time = 1.5)
    expect_equal(det(x), 1)
})


