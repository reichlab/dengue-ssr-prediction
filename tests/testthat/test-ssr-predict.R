library(ssr)
library(magrittr)
library(plyr)
library(dplyr)


context("ssr prediction functions")

test_that("compute_lagged_obs_vecs works -- all variables used", {
    test_data <- data.frame(a = 1:10, b = rnorm(10))
    lags <- list(a = c(1, 5), b = c(0, 2))
    leading_rows_to_drop <- 6
    trailing_rows_to_drop <- 1
    
    expected <- test_data
    expected <- expected %>% mutate(a_lag1=lag(a, 1),
        a_lag5=lag(a, 5),
        b_lag0=b,
        b_lag2=lag(b, 2))
    expected <- expected[seq(from=7, to=9), 3:6]
  
    actual <- compute_lagged_obs_vecs(test_data,
        lags,
        c(seq_len(leading_rows_to_drop),
            seq(from=nrow(test_data) - trailing_rows_to_drop + 1,
                to=nrow(test_data))))
    
    expect_identical(actual, expected)
})

test_that("compute_lagged_obs_vecs works -- one variable not used", {
    test_data <- data.frame(a = 1:10, b = rnorm(10))
    lags <- list(b = c(0, 2))
    leading_rows_to_drop <- 6
    trailing_rows_to_drop <- 0
    
    expected <- test_data
    expected <- expected %>% mutate(b_lag0=b,
        b_lag2=lag(b, 2))
    expected <- expected[seq(from=7, to=10), 3:4]
    
    actual <- compute_lagged_obs_vecs(test_data,
        lags,
        seq_len(leading_rows_to_drop))
    
    expect_identical(actual, expected)
})


test_that("compute_kernel_values works", {
    test_data <- data.frame(a = 1:10, b = rnorm(10))
    lags <- list(a = c(1), b = c(0, 2, 3))
    leading_rows_to_drop <- 3
    trailing_rows_to_drop <- 1
    
    train_lagged_obs <- compute_lagged_obs_vecs(test_data,
        lags,
        c(seq_len(leading_rows_to_drop),
            seq(from=nrow(test_data) - trailing_rows_to_drop + 1,
                to=nrow(test_data))))
    
    prediction_lagged_obs <- compute_lagged_obs_vecs(test_data,
        lags,
        seq_len(9))
    
    bws <- c(1.3, 1.4, 1.5, 0.8)
    expected <- sapply(seq_len(ncol(train_lagged_obs)), function(ind) {
        squared_exp_kernel(x=train_lagged_obs[, ind],
            center=prediction_lagged_obs[, ind],
            bw=bws[ind],
            log=TRUE)
    })
    expected <- apply(expected, 1, sum)
    
    actual <- compute_kernel_values(train_lagged_obs,
        prediction_lagged_obs,
        lags,
        theta=list(a_lag1=list(bw=1.3),
            b_lag0=list(bw=1.4),
            b_lag2=list(bw=1.5),
            b_lag3=list(bw=0.8)),
        control=list(kernel_fns=list(a="squared_exp_kernel",
            b="squared_exp_kernel")),
        log = TRUE)
    
    expect_identical(actual, expected)
})

# test_that("ssr_predict with squared exponential kernel works", {
#     test_data <- data.frame(a = 1:10, b = rnorm(10))
#     lags <- list(a = c(1), b = c(0, 2, 3))
#     leading_rows_to_drop <- 3
#     trailing_rows_to_drop <- 1
#     
#     train_lagged_obs <- compute_lagged_obs_vecs(test_data,
#         lags,
#         leading_rows_to_drop,
#         trailing_rows_to_drop)
#     
#     prediction_lagged_obs <- compute_lagged_obs_vecs(test_data,
#         lags,
#         9,
#         0)
#     
#     bws <- c(1.3, 1.4, 1.5, 0.8)
#     expected <- sapply(seq_len(ncol(train_lagged_obs)), function(ind) {
#         squared_exp_kernel(x=train_lagged_obs[, ind],
#             center=prediction_lagged_obs[, ind],
#             bw=bws[ind],
#             log=TRUE)
#     })
#     expected <- apply(expected, 1, sum)
#     
#     actual <- compute_kernel_values(train_lagged_obs,
#         prediction_lagged_obs,
#         lags,
#         theta=list(a_lag1=list(bw=1.3),
#             b_lag0=list(bw=1.4),
#             b_lag2=list(bw=1.5),
#             b_lag3=list(bw=0.8)),
#         control=list(kernel_fns=list(a="squared_exp_kernel",
#             b="squared_exp_kernel")),
#         log = TRUE)
#     
#     expect_identical(actual, expected)
#     ## slow manually computed result that I know is right
#     lag <- 2
#     prediction_horizon <- 3
#     theta <- 2
#     
#     train_data <- 1:12
#     predict_data <- 2:7
#     lagged_train_data <- cbind(1:10, 2:11, 3:12)
#     lagged_predict_data <- cbind(2:5, 3:6, 4:7)
#     
#     dists <- matrix(NA, nrow=nrow(lagged_train_data), ncol=nrow(lagged_predict_data))
#     for(i in seq_len(nrow(lagged_train_data))) {
#         for(j in seq_len(nrow(lagged_predict_data))) {
#             dists[i, j] <- sqrt(sum((lagged_train_data[i, ] - lagged_predict_data[j, ])^2))
#         }
#     }
#     
#     ## cannot get predictions for the first train_lag elements of predict_data:
#     ## not enough preceeding observations to form lagged obs. vectors
#     ## expect NA weights
#     expected_weights <- matrix(NA, nrow=length(train_data), ncol=length(predict_data))
#     
#     ## to predict prediction_step steps ahead from the remaining values of predict_data,
#     for(j in seq(from=train_lag + 1, to=length(predict_data))) {
#         j_lagged <- j - train_lag  # j indexes predict_data, j_lagged indexes lagged_predict_data
#         
#         dists_for_nonzero_weights <- dists[seq(from=1, length=length(train_data) - train_lag - prediction_step), j_lagged]
#         
#         weights_j <- c(rep(0, train_lag), # no weight to the first train_lag elements of train_data -- there was not enough data preceeding them to form lagged obs. vectors
#             rep(0, prediction_step), # no weight to next prediction_step elements of train_data -- the weights that were assigned to the corresponding lagged obs. vectors are projected forward prediction_step steps
#             exp(-theta * dists_for_nonzero_weights / mean(dists_for_nonzero_weights))
#         )
#         
#         weights_j <- weights_j / sum(weights_j)
#         expected_weights[, j] <- weights_j
#     }
# 
#     ssr_fit <- list(control=list(dist_fn=dist, dist_fn_args=list(method="euclidean")),
#         lag=2,
#         theta=theta)
# 
#     expect_equal(ssr_predict(ssr_fit, train_data, predict_data, prediction_steps=prediction_step)$weights[[1]],
#         expected_weights)
# })
