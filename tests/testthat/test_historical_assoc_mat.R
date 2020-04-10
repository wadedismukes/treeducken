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

    cophylo_pair <- sim_cophylo_bdp(hbr_ = host_lambda,
                                    hdr_ = host_mu,
                                    cosp_rate_ = cosp_rate,
                                    host_exp_rate_ = host_shift_rate,
                                    sdr_ = symb_mu,
                                    sbr_ = symb_lambda,
                                    numbsim_ = numb_replicates,
                                    timeToSimTo_ = time)
    build_historical_association_matrix(t = time/2, cophylo_pair[[1]])
}

test_that("build historical association matrix gives I if only cospeciation is on", {
    x <- hist_assoc_tester(time = 2)
    expect_equal(det(x), det(diag(nrow = nrow(x), ncol = ncol(x))))
})


