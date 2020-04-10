host_mu <- 0.5 # death rate
host_lambda <- 1.0 # birth rate
numb_replicates <- 10
time <- 1.0
symb_mu <- 0.2
symb_lambda <- 0.4
host_shift_rate <- 0.1
cosp_rate <- 1.0

cophylo_pair <- sim_cophylo_bdp(hbr_ = host_lambda,
                          hdr_ = host_mu,
                          cosp_rate_ = cosp_rate,
                          host_exp_rate_ = host_shift_rate,
                          sdr_ = symb_mu,
                          sbr_ = symb_lambda,
                          numbsim_ = numb_replicates,
                          timeToSimTo_ = time)
time <- 0.5
for(i in 1:10){
    assoc_mat_at_t <- build_historical_association_matrix(t = time, tr_pair_obj = cophylo_pair[[i]])
    print(i)
}
