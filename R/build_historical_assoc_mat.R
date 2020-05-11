#' Reconstruct historical association matrix
#'
#' @details Given a time and a tree pair object produced by the `sim_cophylo_bdp`
#'     object will produce the association matrix at that time point for the
#'     tree object
#'
#' @param t The time of interest
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
build_historical_association_matrix <- function(t, tr_pair_obj){

    times <- unique(tr_pair_obj$event_history$Event.Time)
    if(t > max(times)){
        stop("The association matrix you are looking for is the present day one.")
    }
    else{
        events <- tr_pair_obj$event_history
        curr_indx <- min(which(t <= times)) - 1
        init_mat <- matrix(1, nrow = 1, ncol = 1)
        colnames(init_mat) <- as.character(length(tr_pair_obj$host_tree$tip.label) + 1)
        rownames(init_mat) <- as.character(length(tr_pair_obj$symb_tree$tip.label) + 1)
        prev_mat <- init_mat
        for(i in 2:curr_indx){
            curr_events <- subset(events, Event.Time == times[i])
            # what is first row of curr_events
            curr_events
            main_event <- curr_events[1,]$Event.Type
           # print(main_event)
            if(main_event == "HG"){
                # make new matrix with dimensions of prev_mat same matrix
                # remove column indicated by curr_events[1,]$Host.Index
                # make 2 new column vectors name them correctly based on
                # curr_events[2,]$Host.Index whichever rows correspond to 4 in matrx are 1
                # curr_events[3,]$Host.Index
                curr_mat <- matrix(0, nrow = nrow(prev_mat), ncol = ncol(prev_mat) + 1)
                rownames(curr_mat) <- as.character(rownames(prev_mat))
                prev_col_names <- colnames(prev_mat)
                prev_mat <- as.matrix(prev_mat[,-which(colnames(prev_mat) == as.character(curr_events[1,]$Host.Index))])
                if(length(prev_mat) > 0)
                    curr_mat[1:nrow(prev_mat), 1:ncol(prev_mat)] <- prev_mat

                prev_col_names_cut <- prev_col_names[which(prev_col_names != as.character(curr_events[1,]$Host.Index))]
                # if(length(prev_col_names_cut) > 0)
                #     colnames(prev_mat) <- prev_col_names_cut
                col_names <- c(prev_col_names_cut, curr_events[2:3,]$Host.Index)
                colnames(curr_mat) <- as.character(col_names)
                for(j in 1:nrow(curr_events)){
                    type <- curr_events[j,]$Event.Type
                    if(type == "AG"){
                        curr_mat[as.character(curr_events[j,]$Symbiont.Index),
                                 as.character(curr_events[j,]$Host.Index)] <- 1
                    }
                }
            }
            else if(main_event == "HL"){
                curr_mat <- matrix(0, nrow = nrow(prev_mat), ncol = ncol(prev_mat) - 1)
                rownames(curr_mat) <- rownames(prev_mat)
                prev_col_names <- colnames(prev_mat)
                prev_mat <- as.matrix(prev_mat[,-which(colnames(prev_mat) == as.character(curr_events[1,]$Host.Index))])
                prev_col_names_cut <- prev_col_names[which(prev_col_names != as.character(curr_events[1,]$Host.Index))]
                colnames(prev_mat) <- prev_col_names_cut
                curr_mat <- prev_mat
                colnames(curr_mat) <- as.character(prev_col_names_cut)
            }
            else if(main_event == "SG"){
                curr_mat <- matrix(0, nrow = nrow(prev_mat) + 1, ncol = ncol(prev_mat))
                colnames(curr_mat) <- as.character(colnames(prev_mat))
                prev_row_names <- rownames(prev_mat)
                prev_mat <- as.matrix(prev_mat[-which(rownames(prev_mat) == as.character(curr_events[1,]$Symbiont.Index)),], )
                if(length(prev_mat) > 0)
                    curr_mat[1:nrow(prev_mat), 1:ncol(prev_mat)] <- prev_mat

                prev_row_names_cut <- prev_row_names[which(prev_row_names != as.character(curr_events[1,]$Symbiont.Index))]
                # if(length(prev_row_names_cut) > 0)
                #     rownames(prev_mat) <- prev_row_names_cut
                row_names <- c(prev_row_names_cut, curr_events[2:3,]$Symbiont.Index)
                rownames(curr_mat) <- as.character(row_names)
                for(j in 2:nrow(curr_events)){
                    type <- curr_events[j,]$Event.Type
                    if(type == "AG"){
                        curr_mat[as.character(curr_events[j,]$Symbiont.Index),
                                 as.character(curr_events[j,]$Host.Index)] <- 1
                    }
                }
            }
            else if(main_event == "SL"){
                curr_mat <- matrix(0, nrow = nrow(prev_mat) - 1, ncol = ncol(prev_mat))
                colnames(curr_mat) <- colnames(prev_mat)
                prev_row_names <- rownames(prev_mat)
                prev_mat <- as.matrix(prev_mat[-which(rownames(prev_mat) == as.character(curr_events[1,]$Symbiont.Index)),])
                prev_row_names_cut <- prev_row_names[which(prev_row_names != as.character(curr_events[1,]$Symbiont.Index))]
                rownames(prev_mat) <- rownames(prev_row_names_cut)
                curr_mat <- prev_mat
                rownames(curr_mat) <- as.character(prev_row_names_cut)
            }
            else{
                curr_mat <- matrix(0, nrow = nrow(prev_mat) + 1, ncol = ncol(prev_mat) + 1)
                curr_mat[nrow(prev_mat), ncol(prev_mat)] <- 1
                curr_mat[nrow(prev_mat) + 1, ncol(prev_mat) + 1] <- 1
                prev_col_names <- colnames(prev_mat)
                prev_row_names <- rownames(prev_mat)

                prev_mat <- as.matrix(prev_mat[-which(rownames(prev_mat) == as.character(curr_events[1,]$Symbiont.Index)),
                                               -which(colnames(prev_mat) == as.character(curr_events[1,]$Host.Index))])


                prev_col_names_cut <- prev_col_names[which(prev_col_names != as.character(curr_events[1,]$Host.Index))]
                prev_row_names_cut <- prev_row_names[which(prev_row_names != as.character(curr_events[1,]$Symbiont.Index))]

                if(length(prev_mat) > 0)
                    curr_mat[1:nrow(prev_mat), 1:ncol(prev_mat)] <- prev_mat

                # if(length(prev_col_names_cut) > 0)
                #     colnames(prev_mat) <- prev_col_names_cut
                # if(length(prev_row_names_cut) > 0)
                #     rownames(prev_mat) <- prev_row_names_cut
                col_names <- c(prev_col_names_cut, curr_events[2:3,]$Host.Index)
                row_names <- c(prev_row_names_cut, curr_events[2:3,]$Symbiont.Index)
                rownames(curr_mat) <- as.character(row_names)
                colnames(curr_mat) <- as.character(col_names)
                for(j in 2:nrow(curr_events)){
                    type <- curr_events[j,]$Event.Type
                    if(type == "AG"){
                        curr_mat[as.character(curr_events[j,]$Symbiont.Index),
                                 as.character(curr_events[j,]$Host.Index)] <- 1
                    }
                }
            }
            # if C delete both, add 2 new columns and rows
            # else if HG delete col, add 2 new columns
            # else if SG delete col, add 2 new columns
            # else if SL delete col
            # else if HL delete col
            prev_mat <- curr_mat

         #   print(prev_mat)
        }
    }
    prev_mat
}

















