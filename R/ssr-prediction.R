## Functions to perform prediction from new data and an object with a fit
##
## ssr_predict
## ssr_predict_given_lagged_obs
## ssr_kernel_centers_and_weights_predict_given_lagged_obs
## ssr_point_predict_given_lagged_obs
## ssr_dist_predict_given_lagged_lead_obs


#' Make predictions from an estimated ssr model forward prediction_horizon time
#' steps from the end of predict_data, based on the weighting variables, lags,
#' kernel functions, and bandwidths specified in the ssr_fit object.
#' 
#' @param ssr_fit is an object representing a fitted ssr model
#' @param prediction_data is a vector of data points to use in prediction
#' @param prediction_horizon is an integer specifying the number of steps ahead
#'     to perform prediction
#' @param normalize_weights boolean, should the weights be normalized?
#' @param prediction_type character; either "distribution" or "point",
#'     indicating the type of prediction to perform.
#' 
#' @return an object with prediction results; the contents depend on the value
#'     of prediction_type
ssr_predict <- function(ssr_fit,
        prediction_data,
        leading_rows_to_drop=max(ssr_fit$vars_and_lags$lag_value),
        additional_training_rows_to_drop=NULL,
        prediction_horizon,
        prediction_type="distribution",
        n) {
    ## get training and prediction examples
    training_examples <- assemble_training_examples(ssr_fit$train_data,
        ssr_fit$vars_and_lags,
        ssr_fit$y_names,
        leading_rows_to_drop,
        additional_training_rows_to_drop,
        prediction_horizon,
        drop_trailing_rows=TRUE)
    
    prediction_examples <- assemble_training_examples(prediction_data,
        ssr_fit$vars_and_lags,
        ssr_fit$y_names,
        leading_rows_to_drop,
        c(),
        prediction_horizon,
        drop_trailing_rows=FALSE)
    
    ## do prediction
    ssr_predict_given_lagged_obs(training_examples$lagged_obs,
        training_examples$lead_obs,
        prediction_examples$lagged_obs,
        prediction_examples$lead_obs,
        ssr_fit,
        prediction_type,
        n)
}

#' Make predictions from an estimated ssr model forward prediction_horizon time
#' steps from the end of predict_data, based on the kernel functions and
#' bandwidths specified in the ssr_fit object.  This function requires that the
#' lagged and lead observation vectors have already been computed.
#' 
#' @param train_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the training data.  Each row
#'     corresponds to a time point.  Each column is a (lagged) variable.
#' @param train_lead_obs is a vector with length = nrow(train_lagged_obs) with
#'     the value of the prediction target variable corresponding to each row in
#'     the train_lagged_obs matrix.
#' @param prediction_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the prediction data.  There is
#'     only one row, representing one time point.  Each column is a (lagged)
#'     variable.
#' @param ssr_fit is an object representing a fitted ssr model
#' @param prediction_type character; either "distribution" or "point",
#'     indicating the type of prediction to perform.
#' 
#' @return an object with prediction results; the contents depend on the value
#'     of prediction_type
ssr_predict_given_lagged_obs <- function(train_lagged_obs,
    train_lead_obs,
    prediction_lagged_obs,
    prediction_test_lead_obs,
    ssr_fit,
    prediction_type="distribution",
    n) {
    
    if(identical(prediction_type, "centers-and-weights")) {
        kernel_centers_and_weights <-
            ssr_kernel_centers_and_weights_predict_given_lagged_obs(train_lagged_obs,
                train_lead_obs,
                prediction_lagged_obs,
                ssr_fit,
                normalize_weights)
        
        return(kernel_centers_and_weights)
    } else if(identical(prediction_type, "distribution")) {
    	return(ssr_dist_predict_given_lagged_lead_obs(train_lagged_obs,
            train_lead_obs,
            prediction_lagged_obs,
            prediction_test_lead_obs,
            ssr_fit))
    } else if(identical(prediction_type, "point")) {
        kernel_centers_and_weights <-
            ssr_kernel_centers_and_weights_predict_given_lagged_obs(train_lagged_obs,
                train_lead_obs,
                prediction_lagged_obs,
                ssr_fit,
                normalize_weights)
        
        return(ssr_point_predict_given_kernel_centers_and_weights(
            kernel_centers_and_weights = kernel_centers_and_weights,
            ssr_fit = ssr_fit))
    } else if(identical(prediction_type, "sample")) {
        kernel_centers_and_weights <-
            ssr_kernel_centers_and_weights_predict_given_lagged_obs(train_lagged_obs,
                train_lead_obs,
                prediction_lagged_obs,
                ssr_fit,
                normalize_weights)
        
        return(ssr_sample_predict_given_kernel_centers_and_weights(
            n = n,
            kernel_centers_and_weights = kernel_centers_and_weights,
            ssr_fit = ssr_fit))
    } else {
        stop("Invalid prediction type.")
    }
}

#' Make predictions from an estimated ssr model forward prediction_horizon time
#' steps from the end of predict_data, based on the kernel functions and
#' bandwidths specified in the ssr_fit object.  This function requires that the
#' lagged and lead observation vectors have already been computed.
#' 
#' @param train_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the training data.  Each row
#'     corresponds to a time point.  Each column is a (lagged) variable.
#' @param train_lead_obs is a vector with length = nrow(train_lagged_obs) with
#'     the value of the prediction target variable corresponding to each row in
#'     the train_lagged_obs matrix.
#' @param prediction_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the prediction data.  There is
#'     only one row, representing one time point.  Each column is a (lagged)
#'     variable.
#' @param prediction_test_lead_obs is a matrix (with column names) containing
#'     prediction target vectors computed from the prediction data.  Each row
#'     represents one time point.  Each column is a (leading) target variable.
#' @param ssr_fit is an object representing a fitted ssr model
#' @param normalize_weights boolean, should the weights be normalized?
#' 
#' @return a named list with four components:
#'     log_weights: a vector of length = length(train_lagged_obs) with
#'         the log of weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     weights: a vector of length = length(train_lagged_obs) with
#'         the weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     centers: a copy of the train_lead_obs argument -- kernel centers
#'     pred_bws: a named vector of bandwidths, one per column of centers
ssr_kernel_centers_and_weights_predict_given_lagged_obs <- function(train_lagged_obs,
    train_lead_obs,
    prediction_lagged_obs,
    ssr_fit) {
    
    ## compute log of kernel function values representing similarity
    ## of lagged observation vectors from train_data and predict_data.
    ## result is a vector of length nrow(train_lagged_obs)
    unnormalized_log_weights <- compute_kernel_values(train_lagged_obs,
        prediction_lagged_obs,
        kernel_components = ssr_fit$ssr_control$kernel_components,
        theta = ssr_fit$theta_hat,
        log = TRUE)
    
    log_weights <- compute_normalized_log_weights(unnormalized_log_weights)
    
    return(list(unnormalized_log_weights=unnormalized_log_weights,
        log_weights=log_weights,
        weights=exp(log_weights),
        conditioning_vars=train_lagged_obs,
        centers=train_lead_obs))
}


#' Make predictions from an estimated ssr model forward prediction_horizon time
#' steps from the end of predict_data, based on the kernel functions and
#' bandwidths specified in the ssr_fit object.  This function requires that the
#' lagged and lead observation vectors have already been computed.
#' 
#' @param train_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the training data.  Each row
#'     corresponds to a time point.  Each column is a (lagged) variable.
#' @param train_lead_obs is a vector with length = nrow(train_lagged_obs) with
#'     the value of the prediction target variable corresponding to each row in
#'     the train_lagged_obs matrix.
#' @param prediction_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the prediction data.  There is
#'     one row, representing one time point.  Each column is a (lagged) variable.
#' @param prediction_test_lead_obs is a matrix (with column names) containing
#'     prediction target vectors computed from the prediction data.  Each row
#'     represents one time point.  Each column is a (leading) target variable.
#' @param ssr_fit is an object representing a fitted ssr model
#' 
#' @return a list with three components:
#'     log_weights: a vector of length = length(train_lagged_obs) with
#'         the log of weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     weights: a vector of length = length(train_lagged_obs) with
#'         the weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     centers: a copy of the train_lead_obs argument -- kernel centers
ssr_dist_predict_given_lagged_lead_obs <- function(train_lagged_obs,
    train_lead_obs,
    prediction_lagged_obs,
    prediction_test_lead_obs,
    ssr_fit) {
    
    ## compute log of kernel function values representing similarity
    ## of lagged observation vectors from train_data and predict_data.
    ## result is a vector of length nrow(train_lagged_obs)
    log_kernel_values_x <- compute_kernel_values(train_lagged_obs,
        prediction_lagged_obs,
        kernel_components = ssr_fit$ssr_control$kernel_components,
        theta = ssr_fit$theta_hat,
        log = TRUE)
    
    ## compute log(sum(kernel values x))
    log_sum_kernel_values_x <- logspace_sum(log_kernel_values_x)
    
    log_result <- sapply(seq_len(nrow(prediction_test_lead_obs)),
        function(prediction_test_row_ind) {
            log_kernel_values_xy <- compute_kernel_values(
                cbind(train_lagged_obs, train_lead_obs),
                cbind(prediction_lagged_obs,
                    prediction_test_lead_obs),
                kernel_components = ssr_fit$ssr_control$kernel_components,
                theta = ssr_fit$theta_hat,
                log = TRUE
            )
            
            return(logspace_sum(log_kernel_values_xy) -
                    log_sum_kernel_values_x)
        })
    
    return(log_result)
}

#' Make predictions from an estimated ssr model forward prediction_horizon time
#' steps from the end of predict_data, based on the kernel functions and
#' bandwidths specified in the ssr_fit object.  This function requires that the
#' lagged and lead observation vectors have already been computed.
#' 
#' @param train_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the training data.  Each row
#'     corresponds to a time point.  Each column is a (lagged) variable.
#' @param train_lead_obs is a vector with length = nrow(train_lagged_obs) with
#'     the value of the prediction target variable corresponding to each row in
#'     the train_lagged_obs matrix.
#' @param prediction_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the prediction data.  There is
#'     only one row, representing one time point.  Each column is a (lagged)
#'     variable.
#' @param prediction_test_lead_obs is a matrix (with column names) containing
#'     prediction target vectors computed from the prediction data.  Each row
#'     represents one time point.  Each column is a (leading) target variable.
#' @param ssr_fit is an object representing a fitted ssr model
#' @param normalize_weights boolean, should the weights be normalized?
#' 
#' @return a list with three components:
#'     log_weights: a vector of length = length(train_lagged_obs) with
#'         the log of weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     weights: a vector of length = length(train_lagged_obs) with
#'         the weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     centers: a copy of the train_lead_obs argument -- kernel centers
ssr_point_predict_given_kernel_centers_and_weights <- function(kernel_centers_and_weights,
    ssr_fit) {
    stop("point predictions for ssr are not yet implemented")
}

#' Draw a sample from the predictive distribution corresponding to an estimated
#' ssr model forward prediction_horizon time steps from the end of predict_data.
#' This function requires that the kernel weights and centers have already been
#' computed.
#' 
#' @param train_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the training data.  Each row
#'     corresponds to a time point.  Each column is a (lagged) variable.
#' @param train_lead_obs is a vector with length = nrow(train_lagged_obs) with
#'     the value of the prediction target variable corresponding to each row in
#'     the train_lagged_obs matrix.
#' @param prediction_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the prediction data.  There is
#'     only one row, representing one time point.  Each column is a (lagged)
#'     variable.
#' @param prediction_test_lead_obs is a matrix (with column names) containing
#'     prediction target vectors computed from the prediction data.  Each row
#'     represents one time point.  Each column is a (leading) target variable.
#' @param ssr_fit is an object representing a fitted ssr model
#' @param normalize_weights boolean, should the weights be normalized?
#' 
#' @return a list with three components:
#'     log_weights: a vector of length = length(train_lagged_obs) with
#'         the log of weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     weights: a vector of length = length(train_lagged_obs) with
#'         the weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     centers: a copy of the train_lead_obs argument -- kernel centers
ssr_sample_predict_given_kernel_centers_and_weights <- function(n,
    kernel_centers_and_weights,
    ssr_fit) {
    
    result <- matrix(NA, nrow = n, ncol = 1)
    sampled_kernel_inds <- sample(length(kernel_centers_and_weights$weights),
        size = n,
        replace = TRUE,
        prob = kernel_centers_and_weights$weights)
    
    for(kernel_ind in unique(sampled_kernel_inds)) {
        result_inds <- which(sampled_kernel_inds == kernel_ind)
        
        result[result_inds, ] <- simulate_values_from_product_kernel(n = length(result_inds),
            conditioning_obs = kernel_centers_and_weights$conditioning_vars[kernel_ind],
            center = kernel_centers_and_weights$centers[kernel_ind],
            kernel_components = ssr_fit$ssr_control$kernel_components,
            theta = ssr_fit$theta)
    }
    
    return(result)
}
