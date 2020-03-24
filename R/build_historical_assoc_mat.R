


build_historical_association_matrix <- function(time, tr_pair_obj){

    host_tr <- tr_pair[[2]]$host_tree
    symb_tr <- tr_pair[[2]]$symb_tree
    assoc_mat <- tr_pair[[2]]$association_mat
    events <- tr_pair[[2]]$event_history

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

