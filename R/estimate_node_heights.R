
calculate_expected_leaves <- function(lambda,
                                      mu,
                                      t){
    2 * exp((lambda - mu)*t)
}
estimate_node_heights <- function(lambda,
                                  mu,
                                  k = 1,
                                  n){
    if(lambda <= 0.0)
        stop("'lambda' is lower than 0. it must be greater than 0.")
    if(mu > lambda)
        stop("'mu' is greater than 'lambda', please correc this.")
    if(mu == 0.0){
        s <- 0
        for(i in k+1:n){
            s <- s + 1/(lambda*i)
        }
        speciation_event_time <- s
    }
    else if(mu == lambda){
        speciation_event_time <- (n-k)/(lambda*k)
    }
    else{
        rho <- mu / lambda
        first_part <- ((k + 1) / lambda) * exp(lchoose(n, k + 1)) * (-1)^k
        out_sum <- 0
        for(i in 0:(n-k-1)){
            s1 <- exp(lchoose(n - k - 1, i))
            s2 <- 1 / ((k + i + 1) * rho)
            s3 <- exp((k+i) * log(1/rho -1))    # (1/rho - 1)^(k+i)

            in_sum <- 0
            for(j in 1:k+i){
                ins1 <- exp(lchoose(k+i, j))
                ins2 <- (-1^j)/j
                ins3 <- 1 - (1/(1-rho))^j
                print(ins1)
                print(ins2)
                print(ins3)
                in_sum <- in_sum + ins1*ins2*ins3
            }
            s4 <- log(1/(1-rho)) - in_sum
            print(s4)
            out_sum <- out_sum + s1*s2*s3*s4
        }
        speciation_event_time <- first_part * out_sum
    }
    speciation_event_time
}


#estimate_node_heights(1.0, 0.1, 1, 4)
