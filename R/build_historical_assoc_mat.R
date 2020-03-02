


build_historical_association_matrix <- function(time, tr_pair_obj){
    assoc_df <- tr_pair[[1]]$event_history
    event_times <- unique(tr_pair[[1]]$event_history$Event.Time)
    i_time <- 1.0
    t <- 0.0

    assoc_mat <- matrix(1)
    rownames(assoc_mat) <- assoc_df$Symb.Index[1]
    colnames(assoc_mat) <- assoc_df$Host.Index[1]

    time_slice_count <- 1
    while(t < i_time){
        time_slice_count <- time_slice_count + 1
        t <- event_times[time_slice_count]
        crrnt_evnts <- subset(assoc_df, Event.Time <= event_times[time_slice_count])
    }
}