### functions to fit ssr and make predictions

#' Assemble a list of control parameters for the ssr function
#' 
#' @param lag the number of lags to include
#' @param dist_fn is a function to compute distances between vectors
#' @param dist_fn_args is a named list of arguments for dist_fn
#' 
#' @return the (at this point, unvalidated) list of control parameters
ssr_control <- function(lag=1,
        dist_fn=dist,
        dist_fn_args=list(method="euclidean")) {
    control <- list()

    control$lag <- lag

    control$dist_fn <- dist_fn
    control$dist_fn_args <- dist_fn_args

    return(control)
}

#' Estimate the parameters for ssr.
#' This is not yet implemented since we haven't exactly settled on
#' what the parameters are or how we're going to estimate them.
#' 
#' @param train_data a vector of consecutive observations
#' @param control a list of parameters controlling how the fit is done.
#'     See the documentation for ssr_control.
#' 
#' @return an object representing an estimated ssr model
ssr <- function(train_data
        control=ssr_control()) {
    ## do some estimation process in here and return the results?
    lag_hat <- 1
    theta_hat <- 1
    
    return(list(control=control), lag = lag_hat, theta = theta_hat)
}

#' Make predictions from an estimated ssr model.
#' 
#' @param ssr_fit is an object representing a fitted ssr model
#' @param train_data is a vector of training data points, assumed to be
#'     temporally contiguous
#' @param predict_data is a vector of data points to use in prediction
#' @param prediction steps is a vector specifying the number of steps ahead to
#'     perform prediction
#' @param k the maximum number of non-zero weights.
#' 
#' @return a list with two components:
#'     weights: a list with one component for each element of prediction_steps.
#'         component p is a T_train by T_predict matrix, where
#'         T_train = length(train_data) and T_predict = nrow(predict_data).
#'     train_data: a copy of the train_data argument.  
ssr_predict <- function(ssr_fit,
        train_data,
        predict_data,
        prediction_steps,
        k=length(train_data)) {
    ## three steps to prediction:
    ##  1) compute distances between all points of the form
    ##     train_data[i:(i + ssr_fit$lag)] and predict_data[j:(j + ssr_fit$lag)]
    ##  2) translate distances to weights
    ##  3) trace forward to get observations corresponding to each weight.
    
    ## step 1 -- compute distances between lagged observation vectors from
    ## train_data and predict_data.  Entry (i, j) is distance between
    ## train_data[i - lag, ... , i] and predict_data[j - lag, ... , j]
    dists <- compute_pairwise_lagged_obs_distances(train_data, predict_data,
            ssr_fit$lag, ssr_fit$control$dist_fn, ssr_fit$control$dist_fn_args)

    ## steps 2 and 3 -- form weights matrix for each prediction_step.
    ## Entry (i, j) is weight of ith training case for jth prediction case.
    ## Three things to think about:
    ##  1) for lags ssr_fit$lag and prediction_step, only nonzero weights at
    ##     indices i s.t. 
    ##     (a) i > ssr_fit$lag -- because we lag the training data, the first
    ##         few observations can't be used
    ##     (b) i <= length(train_data) - prediction_step -- because we look
    ##         forward to make predictions, the last few observations can't be
    ##         used
    ##  2) enforce the k argument -- at most k non-zero weights
    weights <- lapply(prediction_steps, function(prediction_step) {
        ## prediction_step = prediction lag
        temp <- sapply(seq_len(ncol(dists)), function(j) {
            ## j = prediction case index in lagged observations
            ## step 2 -- compute weights
            
            ## get inds such that we can look ahead the required number
            ## prediction_step of time points:
            ## ssr_fit$lag < i <= length(train_data) - prediction_step
            inds <- seq(from=ssr_fit$lag + 1,
                    to=length(train_data) - prediction_step)
            
            ## enforce at most k non-zero weights -- choose the smallest distances <=> largest weights
            inds <- inds[get_inds_smallest_k(dists[inds - ssr_fit$lag, j],
                    min(k, length(inds)))]
            
            ## weights default to 0, log = -Inf
            log_w_j <- rep(-Inf, length(train_data))

            ## fill in selected values -- push weights forward prediction_step
            ## time units.  I'm following the formulas in Perretti et al. here,
            ## but I think there may be cancellation -- we should work it out.
            ## w_ij = exp(-theta * d_ij / a_j), where
            ## a_j = (1/n) sum_i d_ij is the average distance between
            ## prediction case j and all training points.
            dists_j <- dists[inds - ssr_fit$lag, j]
            log_w_j[inds + prediction_step] <- (-1 * ssr_fit$theta * dists_j /
                            mean(dists_j))
            
            ## re-normalize and return weights for prediction case j
            norm_const <- logspace_sum_matrix_rows(matrix(log_w_j, nrow = 1))
            return(exp(log_w_j - norm_const))
        })

        ## return weights for all prediction cases prediction_step steps ahead
        ## padded by NA columns for prediction cases where there weren't enough
        ## preceeding cases to formed lagged obs. vectors
        return(cbind(matrix(NA, nrow=length(train_data), ncol=ssr_fit$lag),
                temp))
    })
    names(weights) <- paste0("lag_", prediction_steps)

    return(list(weights=weights,
        centers=train_data,
        prediction_steps=prediction_steps))
}

#' This is a wrapper for the ssr_predict function to perform prediction from
#' the dengue data sets for the competition.  Computes KDE predictions for
#' the number of cases in the weeks indexed by t + prediction_step, where
#'   - t is specified by last_obs_season and last_obs_week
#'   - prediction_step varies over the values in prediction_steps
#' Currently the function assumes that data have already been smoothed on the
#' log scale and a data frame variable named smooth_log_cases is available.
#' Log cases are used to do the SSR and get weights for each time point.
#' KDE estimates are then computed both on the log scale using smooth_log_cases
#' and also on the original data scale using exp(smooth_log_cases) as the
#' centers.
#' 
#' @param last_obs_season is the season for the last observed week before we
#'     want to make a prediction.  Has the form "2008/2009"
#' @param last_obs_week is the last observed week of the season specified by
#'     last_obs_season before we want to make a prediction.  An integer
#'     between 1 and 52.
#' @param tr_lag is the number of lags to use in forming the state space
#'     representation via lagged observations
#' @param prediction_steps is a vector giving the number of steps ahead
#'     to do prediction from the week before first_predict_week
#' @param data is the data set to use
#' @param prediction_types is a character vector indicating the prediction types
#'     to perform: may contain one or more of "pt" and "density"
#' 
#' @return list of data frames with predictions.  one list component per entry
#'     in prediction_types argument.  Density estimates are from KDE along a
#'     grid of values for total_counts for each week that we predicted at.
ssr_predict_dengue_stepsahead_one_week <- function(last_obs_season,
    last_obs_week,
    theta,
    tr_lag,
    prediction_steps,
    data,
    prediction_types = c("pt", "density")) {
    ## get indices for prediction data
    ## in order to predict at time t + prediction_step,
    ##  - data must start at time t - tr_lag so that we can get tr_lag + 1
    ##    values in forming the lagged observation vector
    ##  - we need the number of data points equal to tr_lag + 1
    ## first_predict_ind is the first index in the data data frame
    ## that needs to be included in the prediction data
    last_obs_ind <- which(data$season == last_obs_season &
            data$season_week == last_obs_week)
    first_predict_data_ind <- last_obs_ind - tr_lag
    predict_data_inds <- seq(from=first_predict_data_ind, length=tr_lag + 1)
    
    ## for training inds, eliminate anything in the given season(s) where
    ## we are performing prediction.
    predict_seasons <- unique(data$season[last_obs_ind + prediction_steps])
    train_data_inds <- seq_len(nrow(data))[!(data$season %in% predict_seasons)]
    
    ## train_inds may include non-adjacent time intervals,
    ## for example if we are predicting for a season in the middle of the data
    ## perform prediction separately for these non-adjacent intervals, since
    ## ssr_predict currently assumes the provided train_data are all temporally
    ## adjacent in forming lagged observation vectors.
    
    ## get points in the train_data_inds vector with non-adjacent values
    ## THIS SHOULD BE UPDATED TO CHECK season AND season_week VARIABLES, NOT
    ## TRAIN_DATA_INDS.  If there are multiple non-adjacent time points that have
    ## been dropped (e.g., you're doing 2 levels of cross validation),
    ## the current code fails.
    split_pts <- which((train_data_inds[- length(train_data_inds)] + 1) !=
            train_data_inds[-1])
    split_pts <- c(0, split_pts, length(train_data_inds))
    
    ## form list with one component for each group of adjacent train data inds
    train_data_inds_by_segment <- lapply(seq_len(length(split_pts) - 1),
        function(i) {
            return(train_data_inds[seq(from = split_pts[i] + 1,
                        to = split_pts[i + 1])])
        }
    )
    
    ## do prediction for each segment of contiguous training inds
    ssr_pred_by_segment <- lapply(train_data_inds_by_segment,
        function(restricted_train_data_inds) {
            return(ssr_predict(
                    ssr_fit = list(control = ssr_control(),
                        theta = theta,
                        lag = tr_lag),
                    train_data = data$smooth_log_cases[restricted_train_data_inds],
                    predict_data = data$smooth_log_cases[predict_data_inds],
                    prediction_steps = prediction_steps))
        }
    )
    
    ## combine prediction results from different segments of training inds
    ssr_pred_combined <- list(
        ## weights are list by prediction step of combined weights matrices
        weights = lapply(seq_along(prediction_steps),
            function(psi) {
                rbind.fill.matrix(lapply(ssr_pred_by_segment,
                        function(li) {
                            li$weights[[psi]]
                        }))
            }),
        ## centers are vector of combined center vectors
        centers = unlist(lapply(ssr_pred_by_segment,
                function(li) {
                    li$centers
                }))
    )
    
    ## a function to get point predictions for one week
    get_pt_predictions_one_week <- function(prediction_step_i) {
        # get weighted mean
        pt_est_log_scale <- weighted.mean(ssr_pred_combined$centers,
            ssr_pred_combined$weights[[prediction_step_i]][, 2])
        pt_est_orig_scale <- weighted.mean(exp(ssr_pred_combined$centers),
            ssr_pred_combined$weights[[prediction_step_i]][, 2])
        
        # time at which prediction was performed
        season_and_week <- get_prediction_season_week(last_obs_ind,
            prediction_steps[prediction_step_i], data)
        
        return(data.frame(est_total_cases_log_scale = pt_est_log_scale,
                est_total_cases_orig_scale = pt_est_orig_scale,
                est_total_cases_orig_scale_from_log_scale = exp(pt_est_log_scale),
                season = season_and_week$season,
                season_week = season_and_week$week,
                prediction_step = prediction_steps[prediction_step_i]
            ))
    }
    
    ## a function to get kde predictions for one week
    get_kde_predictions_one_week <- function(prediction_step_i) {
        ## do weighted kde
        kde_est_log_scale <- density(ssr_pred_combined$centers,
            weights = ssr_pred_combined$weights[[prediction_step_i]][, 2],
            bw = "SJ")
        kde_est_orig_scale <- density(exp(ssr_pred_combined$centers),
            weights = ssr_pred_combined$weights[[prediction_step_i]][, 2],
            bw = "SJ")
        
        # time at which prediction was performed
        season_and_week <- get_prediction_season_week(last_obs_ind,
            prediction_steps[prediction_step_i], data)
        
        return(data.frame(log_total_cases = kde_est_log_scale$x,
                total_cases = kde_est_orig_scale$x,
                est_density_log_scale = kde_est_log_scale$y,
                est_density_orig_scale = kde_est_orig_scale$y,
                est_density_orig_scale_from_log_scale = kde_est_log_scale$y / kde_est_orig_scale$x,
                season = season_and_week$season,
                season_week = season_and_week$week,
                prediction_step = prediction_steps[prediction_step_i]
            ))
    }
    
    ## get predictions for all specified prediction_steps
    result <- list()
    
    if("pt" %in% prediction_types) {
        result$pt_preds <- rbind.fill(
            lapply(seq_along(prediction_steps), get_pt_predictions_one_week)
        )
        result$pt_preds$prediction_step <-
            factor(result$pt_preds$prediction_step)
    }
    
    if("density" %in% prediction_types) {
        result$density_preds <- rbind.fill(
            lapply(seq_along(prediction_steps), get_kde_predictions_one_week)
        )
        result$density_preds$prediction_step <-
            factor(result$density_preds$prediction_step)
    }
    
    return(result)
}


#' Get week and season when we're predicting based on the index of the last obs
#' and the number of steps forward we're predicting
#' 
#' @param last_obs_ind index in the data data frame of the last observation
#'     before we start predicting
#' @param prediction_step the number of steps forward we're predicting
#' @param data the data frame, with columns named season_week and season
#' 
#' @return a list with two components indicating the time of prediction:
#'     1) week is an integer from 1 to 52 and
#'     2) season is a string of the form "2008/2009"
get_prediction_season_week <- function(last_obs_ind, prediction_step, data) {
    wk <- (data$season_week[last_obs_ind] + prediction_step) %% 52
    if(wk == 0) {
        wk <- 52
    }
    seasons_advanced <- (data$season_week[last_obs_ind] + prediction_step
            - wk) / 52
    start_season_last_obs <- as.integer(substr(
            as.character(data$season[last_obs_ind]),
            start = 1,
            stop = 4))
    season <- start_season_last_obs + seasons_advanced
    season <- paste0(season, "/", season + 1)
    
    return(list(week = wk, season = season))
}

#' Get the indices of the smallest k elements of v.  This code currently assumes
#' that k >= length(v)
#' 
#' @param v a vector
#' @param k number of indices to return
#' 
#' @return a vector of length k containing the indices of the k smallest
#'   elements of v, in ascending order.
get_inds_smallest_k <- function(v, k) {
    return(order(v, decreasing=FALSE)[seq_len(k)])
}

#' Convert an observation vector to a matrix of lagged observation vectors
#'
#' @param orig_data a vector to convert
#' @param lag the number of lags
get_lagged_obs_matrix <- function(orig_data, lag) {
    ## orig_data is a vector, lag is an integer >= 0
    n <- length(orig_data)

    lag <- as.integer(lag)
    if(is.na(lag) || lag < 0 || lag > n - 1) {
        stop("invalid lag")
    }

    ## there's probably a much cleaner way to do this
    if(lag == n - 1) {
        dim(orig_data) <- c(1, n)
        return(orig_data)
    } else {
        return(sapply(seq(from=0, to=lag),
            function(l) {
                orig_data[seq(from=1 + l, to=n - lag + l)]
            }))
    }
}

#' Compute pairwise distances between the lagged observation vectors 
#' formed from v1 and v2.
#' 
#' @param v1 first vector
#' @param v2 second vector
#' @param lag number of lags to use in SSR
#' @param dist_fn a function to compute distances between vectors.  It should
#'     accept an argument x, a matrix with rows the vectors to compute the
#'     distance between
#' @param dist_fn_args a named list of arguments to dist_fn
#' 
#' @return a matrix, where entry (i, j) is the distance between
#'     v1[i, ..., i + lag] and v2[j, ..., j + lag]
compute_pairwise_lagged_obs_distances <- function(v1,
        v2,
        lag,
        dist_fn,
        dist_fn_args) {
    ## this is a first pass -- wastes memory and processing power,
    ## but easy to understand and change

    ## form matrices with lagged observations
    m1 <- get_lagged_obs_matrix(v1, lag)
    m2 <- get_lagged_obs_matrix(v2, lag)

    ## compute result
    return(outer(seq_len(nrow(m1)), seq_len(nrow(m2)),
        FUN=Vectorize(function(i, j) {
            dist_fn_args$x <- rbind(m1[i, ], m2[j, ])
            do.call(dist_fn, dist_fn_args)
        })))
}
