library(ssr)

context("ssr prediction functions")

test_that("compute_pairwise_lagged_obs_distances works", {
    ## slow manually computed result that I know is right, using Euclidean distance
    testa <- cbind(1:3, 2:4, 3:5)
    testb <- cbind(2:8, 3:9, 4:10)
    testres <- matrix(NA, nrow=nrow(testa), ncol=nrow(testb))
    for(i in seq_len(nrow(testa))) {
        for(j in seq_len(nrow(testb))) {
            testres[i, j] <- sqrt(sum((testa[i, ] - testb[j, ])^2))
        }
    }

    expect_equal(compute_pairwise_lagged_obs_distances(1:5, 2:10, lag=2, dist_fn=dist, dist_fn_args=list()),
        testres)
})

test_that("ssr_predict with euclidean distance, no k (max nonzero weights) works", {
    ## slow manually computed result that I know is right, using Euclidean distance
    train_lag <- 2
    prediction_step <- 3
    theta <- 2
    
    train_data <- 1:12
    predict_data <- 2:7
    lagged_train_data <- cbind(1:10, 2:11, 3:12)
    lagged_predict_data <- cbind(2:5, 3:6, 4:7)
    
    dists <- matrix(NA, nrow=nrow(lagged_train_data), ncol=nrow(lagged_predict_data))
    for(i in seq_len(nrow(lagged_train_data))) {
        for(j in seq_len(nrow(lagged_predict_data))) {
            dists[i, j] <- sqrt(sum((lagged_train_data[i, ] - lagged_predict_data[j, ])^2))
        }
    }
    
    ## cannot get predictions for the first train_lag elements of predict_data:
    ## not enough preceeding observations to form lagged obs. vectors
    ## expect NA weights
    expected_weights <- matrix(NA, nrow=length(train_data), ncol=length(predict_data))
    
    ## to predict prediction_step steps ahead from the remaining values of predict_data,
    for(j in seq(from=train_lag + 1, to=length(predict_data))) {
        j_lagged <- j - train_lag  # j indexes predict_data, j_lagged indexes lagged_predict_data
        
        dists_for_nonzero_weights <- dists[seq(from=1, length=length(train_data) - train_lag - prediction_step), j_lagged]
        
        weights_j <- c(rep(0, train_lag), # no weight to the first train_lag elements of train_data -- there was not enough data preceeding them to form lagged obs. vectors
            rep(0, prediction_step), # no weight to next prediction_step elements of train_data -- the weights that were assigned to the corresponding lagged obs. vectors are projected forward prediction_step steps
            exp(-theta * dists_for_nonzero_weights / mean(dists_for_nonzero_weights))
        )
        
        weights_j <- weights_j / sum(weights_j)
        expected_weights[, j] <- weights_j
    }

    ssr_fit <- list(control=list(dist_fn=dist, dist_fn_args=list(method="euclidean")),
        lag=2,
        theta=theta)

    expect_equal(ssr_predict(ssr_fit, train_data, predict_data, prediction_steps=prediction_step)$weights[[1]],
        expected_weights)
})
