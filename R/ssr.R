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
#' @param train_data
#' @param predict_data
#' @param control
#' 
#' @return an object representing an estimated ssr model
ssr <- function(train_data, predict_data, control=ssr_control()) {

    ## do some estimation process in here and return the results?

    return(list(control=control))
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
ssr_predict <- function(ssr_fit, train_data, predict_data, prediction_steps, k=length(train_data)) {
    ## three steps to prediction:
    ##  1) compute distances between all points of the form train_data[i:(i + ssr_fit$lag)] and predict_data[j:(j + ssr_fit$lag)]
    ##  2) translate distances to weights
    ##  3) trace forward to get observations corresponding to each weight.
    
    ## step 1 -- compute distances between lagged observation vectors from train_data and predict_data.
    ## Entry (i, j) is distance between train_data[i - lag, ... , i] and predict_data[j - lag, ... , j]
    dists <- compute_pairwise_lagged_obs_distances(train_data, predict_data, ssr_fit$lag, ssr_fit$control$dist_fn, ssr_fit$control$dist_fn_args)

    ## steps 2 and 3 -- form weights matrix for each prediction_step.
    ## Entry (i, j) is weight of ith training case for jth prediction case.
    ## Three things to think about:
    ##  1) for lags ssr_fit$lag and prediction_step, only nonzero weights at indices i s.t. 
    ##     (a) i > ssr_fit$lag -- because we lag the training data, the first few observations can't be used
    ##     (b) i <= length(train_data) - prediction_step -- because we look ahead to make predictions, the last few observations can't be used
    ##  2) enforce the k argument -- at most k non-zero weights
    weights <- lapply(prediction_steps, function(prediction_step) { # prediction_step = prediction lag
        temp <- sapply(seq_len(ncol(dists)), function(j) { # j = prediction case index in lagged observations
            ## step 2 -- compute weights
            
            ## get inds such that we can look ahead the required number prediction_step of time points:
            ## ssr_fit$lag < i <= length(train_data) - prediction_step
            inds <- seq(from=ssr_fit$lag + 1, to=length(train_data) - prediction_step)
            
            ## enforce at most k non-zero weights -- choose the smallest distances <=> largest weights
            inds <- inds[get_inds_smallest_k(dists[inds - ssr_fit$lag, j], min(k, length(inds)))]
            
            ## weights default to 0
            w_j <- rep(0, length(train_data))

            ## fill in selected values -- push weights forward prediction_step time units
            ## I'm following the formulas in Perretti et al. here, but I think there may be cancellation -- we should work it out.
            ## w_ij = exp(-theta * d_ij / a_j), where a_j = (1/n) sum_i d_ij is the average distance between prediction case j and all training points.
            dists_j <- dists[inds - ssr_fit$lag, j]
            w_j[inds + prediction_step] <- exp(-1 * ssr_fit$theta * dists_j / mean(dists_j))
            
            ## re-normalize and return weights for prediction case j
            return(w_j / sum(w_j))
        })

        ## return weights for all prediction cases prediction_step steps ahead
        ## padded by NA columns for prediction cases where there weren't enough preceeding cases to formed lagged obs. vectors
        return(cbind(matrix(NA, nrow=length(train_data), ncol=ssr_fit$lag), temp))
    })
    names(weights) <- paste0("lag_", prediction_steps)

    return(list(weights=weights, centers=train_data, prediction_steps=prediction_steps))
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

    return(sapply(seq(from=0, to=lag),
        function(l) {
            orig_data[seq(from=1 + l, to=n - lag + l)]
        }))
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
    ## this is a first pass -- wastes memory and processing power, but easy to understand and change

    ## form matrices with lagged observations
    m1 <- get_lagged_obs_matrix(v1, lag)
    m2 <- get_lagged_obs_matrix(v2, lag)

    ## m1 and m2 are matrices with the same number of columns, dist_fn is a function to compute distances, dist_fn_args are arguments to dist_fn
    ## entry (i, j) of the result has the distance between row i of train and row j of predict.
    return(outer(seq_len(nrow(m1)), seq_len(nrow(m2)),
        FUN=Vectorize(function(i, j) {
            dist_fn_args$x <- rbind(m1[i, ], m2[j, ])
            do.call(dist_fn, dist_fn_args)
        })))
}
