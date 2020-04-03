#' Reconstruct historical association matrix
#'
#' @details Given a time and a tree pair object produced by the `sim_cophylo_bdp`
#'     object will produce the association matrix at that time point for the
#'     tree object
#'
#' @param time The time of interest
#' @param tr_pair_obj The tree pair object from `sim_cophylo_bdp`
#' @return Matrix of the associations at given time
#' @examples
#' host_mu <- 0.5 # death rate
#' host_lambda <- 2.0 # birth rate
#' numb_replicates <- 10
#' time <- 2.9
#' symb_mu <- 0.2
#' symb_lambda <- 0.4
#' host_shift_rate <- 0.1
#' cosp_rate <- 2.0
#'
#' cophylo_pair <- sim_cophylo_bdp(hbr_ = host_lambda,
#'                            hdr_ = host_mu,
#'                            cosp_rate_ = cosp_rate,
#'                            host_exp_rate_ = host_shift_rate,
#'                            sdr_ = symb_mu,
#'                            sbr_ = symb_lambda,
#'                            numbsim_ = numb_replicates,
#'                            timeToSimTo_ = time)
#' t <- 1.4
#' assoc_mat_at_t <- build_historical_association_matrix(time = t, tr_pair_obj = cophylo_pair[[1]])
#'
build_historical_association_matrix <- function(time, tr_pair_obj){

    host_tr <- tr_pair_obj$host_tree
    symb_tr <- tr_pair_obj$symb_tree
    assoc_mat <- tr_pair_obj$association_mat
    events <- tr_pair_obj$event_history

    event_times <- unique(events$Event.Time)
    i_time <- 1.0
    prev_mat <- matrix(1, nrow = 1, ncol = 1)



    time_slice_count <- 1
    t <- 0.0
    while(t < i_time){
        curr_events <- subset(events, Event.Time == event_times[time_slice_count])
        t <- event_times[time_slice_count]

        hosts_in_slice <- unique(curr_events$Host.Index)
        symbs_in_slice <- unique(curr_events$Symb.Index)
        curr_mat <- matrix(0, nrow = length(symbs_in_slice),
                           ncol = length(hosts_in_slice))
        curr_mat[1:nrow(prev_mat), 1:ncol(prev_mat)]


        curr_mat
        time_slice_count <- time_slice_count + 1
    }
    curr_mat

}

















