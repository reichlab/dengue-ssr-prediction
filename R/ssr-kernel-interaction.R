## Functions for interaction between ssr estimation and prediction routines and kernel functions
## 
## initialize_theta
## extract_vectorized_theta_est_from_theta
## update_theta_from_vectorized_theta_est
## compute_kernel_values
## simulate_values_from_product_kernel

#' Initialize parameter values.
#'
#' @param prev_theta previous theta values before update to variables and lags
#' @param updated_vars_and_lags lags after update
#' @param update_var_name character; the name of the variable added/removed from the model
#' @param update_lag_value integer; the lag that was added/removed from the model
#' @param data data matrix
#' @param ssr_control list of control parameters for ssr
#' 
#' @return list of theta parameters
initialize_theta <- function(prev_theta,
    updated_vars_and_lags,
    update_var_name,
    update_lag_value,
    data,
    ssr_control) {
    theta <- lapply(seq_along(ssr_control$kernel_components),
        function(ind) {
            updated_var_and_lag_in_current_component <-
                any(ssr_control$kernel_components[[ind]]$vars_and_lags$var_name == update_var_name &
                        ssr_control$kernel_components[[ind]]$vars_and_lags$lag_value == update_lag_value)
            if(updated_var_and_lag_in_current_component ||
                is.null(prev_theta[[ind]])) {
                potential_cols_in_component <-
                    ssr_control$kernel_components[[ind]]$vars_and_lags$combined_name
                cols_used <- colnames(data) %in% potential_cols_in_component
                if(length(cols_used) > 0) {
                    fn_name <- ssr_control$kernel_components[[ind]]$
                        initialize_kernel_params_fn
                    
                    fn_args <- ssr_control$kernel_components[[ind]]$
                        initialize_kernel_params_args
                    fn_args$update_var_name <- update_var_name
                    fn_args$update_lag_value <- update_lag_value
                    fn_args$prev_theta <- prev_theta[[ind]]
                    fn_args$ssr_control <- ssr_control
                    fn_args$x <- data[, cols_used, drop=FALSE]
                    
                    return(do.call(fn_name, fn_args))
                } else {
                    return(NULL)
                }
            } else {
                return(prev_theta[[ind]])
            }
        }
    )
    
    return(theta)
}


#' Extract a vector of parameter values that are to be estimated from theta,
#' on the estimation scale.
#' 
#' @param theta_list kernel parameters theta in list form
#' @param lags list representing combinations of variables and lags
#'     included in the model
#' @param ssr_control control parameters for the ssr fit
#' @param add_fixed_params boolean -- should parameters that are being held
#'     fixed in the estimation process be added to the return value?
#' 
#' @return numeric vector with parameter values
extract_vectorized_theta_est_from_theta <- function(theta,
	vars_and_lags,
    ssr_control,
    add_fixed_params = FALSE) {
    
    theta_vector <- c()
    
    for(ind in seq_along(ssr_control$kernel_components)) {
        ## parameters that are being estimated
        theta_vector <- c(theta_vector, do.call(
			ssr_control$kernel_components[[ind]]$vectorize_kernel_params_fn,
       		c(list(theta_list = theta[[ind]],
				vars_and_lags = vars_and_lags),
       			ssr_control$kernel_components[[ind]]$vectorize_kernel_params_args)
       	))
    }

    return(theta_vector)
}

#' Convert theta from vector form to list form
#' 
#' @param theta_est_vector numeric vector of kernel parameters theta that are
#'     being estimated, on estimation scale.
#' @param ssr_control control parameters for the ssr fit
#' 
#' @return list of lists of parameter values -- outer list has one component
#'     for each kernel function, inner list has one component
#'     for each parameter used in the corresponding kernel function
update_theta_from_vectorized_theta_est <- function(theta_est_vector,
	theta,
    ssr_control) {
    
    theta_vector_component_start_ind <- 1L
    for(ind in seq_along(ssr_control$kernel_components)) {
        ## parameters that are being estimated
        temp <- do.call(
			ssr_control$kernel_components[[ind]]$
				update_theta_from_vectorized_theta_est_fn,
        	c(list(theta_est_vector = theta_est_vector[
					seq(from = theta_vector_component_start_ind,
        				to = length(theta_est_vector))],
        		theta = theta[[ind]]),
        		ssr_control$kernel_components[[ind]]$
					update_theta_from_vectorized_theta_est_args)
        )
        
        theta_vector_component_start_ind <- theta_vector_component_start_ind +
			temp$num_theta_vals_used
        
        theta[[ind]] <- temp$theta
    }

    return(theta)
}


#' Compute kernel values measuring the similarity of each row in the
#' train_lagged_obs data frame to the prediction_lagged_obs.
#' 
#' @param train_obs a data frame with lagged observation vectors computed
#'     from the training data
#' @param prediction_obs a data frame with the lagged observation vector
#'     computed from the prediction data.  It is assumed that
#'     prediction_lagged_obs contains only one row.
#' @param theta a list with one component for each component of the kernel
#'     function.  This component is a named list with arguments to the
#'     corresponding kernel function.
#' @param ssr_control a list of ssr_control parameters for ssr
#' @param log boolean; if TRUE (default), return kernel values on the log scale
compute_kernel_values <- function(train_obs,
    prediction_obs,
    kernel_components,
    theta,
    log = TRUE) {
	if(!(identical(nrow(prediction_obs), 1L))) {
		stop("In call to compute_kernel_values, prediction_obs must have exactly 1 row.")
	}

    ## create a matrix of log kernel values by component
	## rows correspond to time points in train_obs, columns to components of
	## the kernel function
    log_kernel_component_values <- matrix(0,
        nrow=nrow(train_obs),
        ncol=length(kernel_components))
    
    for(ind in seq_along(kernel_components)) {
		combined_names_in_component <-
			kernel_components[[ind]]$vars_and_lags$combined_name
        col_names <- colnames(train_obs)[
            colnames(train_obs) %in% combined_names_in_component]
        
        if(length(col_names) > 0) {
	        ## assemble arguments to kernel function
	        kernel_fn_args <- c(theta[[ind]],
	            kernel_components[[ind]]$theta_fixed)
	        kernel_fn_args$x <- train_obs[, col_names, drop = FALSE]
	        kernel_fn_args$center <- prediction_obs[, col_names, drop = FALSE]
	        kernel_fn_args$log <- TRUE
	        
	        ## call kernel function and store results
	        log_kernel_component_values[, ind] <-
	            do.call(kernel_components[[ind]]$kernel_fn,
					kernel_fn_args)
	    }
    }
    
    ## return on scale requested by user
    ## these computations assume product kernel --
    ## if we're doing something else, change apply(..., sum)
    if(log) {
        return(apply(log_kernel_component_values, 1, sum))
    } else {
        return(exp(apply(log_kernel_component_values, 1, sum)))
    }
}



#' Simulate values from (conditional) kernel.
#' 
#' @param train_obs a data frame with lagged observation vectors computed
#'     from the training data
#' @param prediction_obs a data frame with the lagged observation vector
#'     computed from the prediction data.  It is assumed that
#'     prediction_lagged_obs contains only one row.
#' @param theta a list with one component for each component of the kernel
#'     function.  This component is a named list with arguments to the
#'     corresponding kernel function.
#' @param ssr_control a list of ssr_control parameters for ssr
#' @param log boolean; if TRUE (default), return kernel values on the log scale
simulate_values_from_product_kernel <- function(n,
    conditioning_obs,
    center,
    kernel_components,
    theta) {
    if(!(identical(nrow(conditioning_obs), 1L)) || !(identical(nrow(center), 1L))) {
        stop("In call to simulate_values_from_product_kernel, conditioning_obs and center must have exactly 1 row.")
    }
    
    ## create a matrix to contain simulated values
    ## rows correspond to simulations, columns to variables simulated
    simulated_values <- matrix(NA,
        nrow = n,
        ncol = ncol(center))
    colnames(simulated_values) <- colnames(center)
    
    for(ind in seq_along(kernel_components)) {
        combined_names_in_component <-
            kernel_components[[ind]]$vars_and_lags$combined_name
        conditioning_col_names <- colnames(conditioning_obs)[
            colnames(conditioning_obs) %in% combined_names_in_component]
        center_col_names <- colnames(center)[
            colnames(center) %in% combined_names_in_component]
        
        if(length(conditioning_col_names) > 0) {
            simulated_values[, conditioning_col_names] <-
                rep(conditioning_obs[, conditioning_col_names], each = n)
        }
        
        if(length(center_col_names) > 0) {
            ## assemble arguments to kernel function
            rkernel_args <- c(theta[[ind]],
                kernel_components[[ind]]$theta_fixed)
            rkernel_args$n <- n
            rkernel_args$conditioning_obs <-
                conditioning_obs[, conditioning_col_names, drop = FALSE]
            rkernel_args$center <- center[, center_col_names, drop = FALSE]
            
            ## call kernel function and store results
            simulated_values[, center_col_names] <-
                do.call(kernel_components[[ind]]$rkernel_fn,
                    rkernel_args)
        }
    }
    
    return(simulated_values)
}
