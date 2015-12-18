### functions to fit ssr and make predictions and forecasts

### functions for creating parameters specifying how the fit should be performed

#' Assemble a list of ssr_control parameters for the ssr function with
#'     user-specified values.
#' 
#' @param X_names a character vector of length >= 1 containing names of
#'     variables in the data data frame to use in forming the lagged
#'     observation process used for calculating weights
#' @param y_names a character vector of length 1 containing the name of the
#'     variable in the data data frame to use as the target for prediction
#' @param time_name (optional) a character vector of length 1 containing the
#'     name of the variable in the data data frame to use as the time.
#' @param prediction_horizons integer vector: the number of time steps between
#'     the last observation and the time at which we make a prediction
#' @param kernel_components a list with one component for each component of the
#'     kernel function.  Each component is a list with the following entries:
#'       - vars_and_lags: a data frame with two columns: "var_name" and
#'         "lag_value".  Each row specifies a combination of variable and lag
#'         that is included in this component of the kernel function.
#'       - kernel_fn: a function to evaluate the kernel
#'       - theta_fixed: a named list of parameters to kernel_fn whose values are
#'         held fixed (i.e., not estimated)
#'       - theta_est: a named list of parameters to kernel_fn whose values are
#'         to be estimated
#'       - initialize_theta_fn: a function to initialize both theta_fixed and
#'         theta_est
#'       - initialize_theta_args: a named list of arguments to
#'         initialize_theta_fn
#'       - vectorize_theta_est: a function that converts theta_est into an
#'         ordered vector on a scale suitable for passing as the first argument
#'         to optim.  Required to return a list with three components:
#'           (1) theta_est - vector of parameters to be estimated
#'           (2) lb - vector of lower bounds to theta_est
#'           (3) ub - vector of upper bounds to theta_est
#'       - update_theta_from_vectorized_theta_est: a function that updates
#'         theta_est (in list form) from theta_est (in vector form).
#' @param crossval_buffer during cross-validation, the number of indices before
#'     the time at which we are making a prediction to drop from the "training
#'     examples".
#' @param loss_fn_name a string giving the name of the function use to
#'     compute loss from predictions
#' @param loss_fn_args a named list giving arguments to the loss function
#' @param par_packages a character vector containing names of packages that need
#'     to be loaded in instances of R when computations are performed in
#'     parallel.
#' 
#' @return the (at this point, unvalidated) list of ssr_control parameters
create_ssr_control <- function(X_names,
        y_names,
        time_name,
        prediction_horizons,
        kernel_components,
        crossval_buffer,
        loss_fn_name,
        loss_fn_args,
        par_packages = NULL) {
    ssr_control <- list()

    ssr_control$X_names <- X_names
    ssr_control$y_names <- y_names
    ssr_control$time_name <- time_name
    
    ssr_control$prediction_horizons <- prediction_horizons
    
    ssr_control$kernel_components <- kernel_components
    
    ssr_control$crossval_buffer <- crossval_buffer
    
    ssr_control$loss_fn_name <- loss_fn_name
    ssr_control$loss_fn_args <- loss_fn_args

    return(ssr_control)
}

#' Assemble a list of ssr_control parameters for the ssr function with default
#'     values
#' 
#' @param X_names a character vector of length >= 1 containing names of
#'     variables in the data data frame to use in forming the lagged
#'     observation process used for calculating weights
#' @param y_names a character vector of length 1 containing the name of the
#'     variable in the data data frame to use as the target for prediction
#' @param time_name (optional) a character vector of length 1 containing the
#'     name of the variable in the data data frame to use as the time.
#' @param data a data frame where rows are consecutive observations
#' 
#' @return the list of ssr_control parameters
create_ssr_control_default <- function(X_names, y_names, time_name, data) {
    ssr_control <- list()
    
    ssr_control$X_names <- X_names
    ssr_control$y_names <- y_names
    ssr_control$time_name <- time_name
    
    ssr_control$kernel_components <- get_default_kernel_components(X_names,
        y_names,
        time_name,
        data)
    
    ssr_control$loss_fn_name <- "mase"
    ssr_control$loss_fn_args <- list()
    
    return(ssr_control)
}

#' Get default kernel functions based on a brief look at the data.  This is
#' unreliable.  Update to return periodic_kernel if X_names[i] == time_name?
#' 
#' @param X_names a character vector of length >= 1 containing names of
#'     variables in the data data frame to use in forming the lagged
#'     observation process used for calculating weights
#' @param y_names a character vector of length 1 containing the name of the
#'     variable in the data data frame to use as the target for prediction
#' @param time_name (optional) a character vector of length 1 containing the
#'     name of the variable in the data data frame to use as the time.
#' @param data a data frame where rows are consecutive observations
#' 
#' @return a list of default parameters for kernel components -- probably all bad
get_default_kernel_components <- function(X_names, y_names, time_name, data) {
	stop("Function get_default_kernel_components is not yet implemented")

    return(kernel_components)
}

#' Validate ssr_control parameters for ssr -- not implemented
#' 
#' @param ssr_control a list of ssr_control parameters for ssr
#' @param X_names a character vector of length >= 1 containing names of
#'     variables in the data data frame to use in forming the lagged
#'     observation process used for calculating weights
#' @param y_names a character vector of length 1 containing the name of the
#'     variable in the data data frame to use as the target for prediction
#' @param time_name (optional) a character vector of length 1 containing the
#'     name of the variable in the data data frame to use as the time.
#' @param data a data frame where rows are consecutive observations
#' 
#' @return no return value -- either stops with an error or not.
validate_ssr_control <- function(ssr_control, X_names, y_names, time_name, data) {
#    warning("ssr ssr_control parameter validation not yet implemented")
}


### functions for parameter estimation

#' Estimate the parameters for ssr.  There is redundancy here in that X_names,
#' y_names, time_name are all included in the ssr_control object as well as
#' parameters to this function.  Decide on an interface.
#' 
#' @param X_names a character vector of length >= 1 containing names of
#'     variables in the data data frame to use in forming the lagged
#'     observation process used for calculating weights
#' @param y_names a character vector of length 1 containing the name of the
#'     variable in the data data frame to use as the target for prediction
#' @param time_name (optional) a character vector of length 1 containing the
#'     name of the variable in the data data frame to use as the time.
#' @param data a data frame where rows are consecutive observations
#' @param ssr_control a list of parameters ssr_controlling how the fit is done.
#'     See the documentation for ssr_control.
#' 
#' @return an object representing an estimated ssr model.  Currently, a list
#'     with 7 components.
ssr <- function(X_names,
        y_names,
        time_name,
        data,
        ssr_control) {
    ## get/validate ssr_control argument
    if(missing(ssr_control)) {
        ssr_control <- create_ssr_control_default(X_names, y_names, time_name, data)
        warning("ssr_control argument not supplied to ssr -- using defaults, which may be bad")
    } else {
        validate_ssr_control(ssr_control, X_names, y_names, time_name, data)
    }
    
    ## estimate lags and kernel parameters via cross-validation
    param_estimates <- est_ssr_params_stepwise_crossval(data, ssr_control)
    
    return(list(ssr_control=ssr_control,
        X_names=X_names,
        y_names=y_names,
        time_name=time_name,
		vars_and_lags=param_estimates$vars_and_lags,
        theta_hat=param_estimates$theta_hat,
        train_data=data))
}

#' Use a stepwise procedure to estimate parameters and select lags by optimizing
#' a cross-validation estimate of predictive performance
#' 
#' @param data the data frame to use in performing cross validation
#' @param ssr_control a list of parameters specifying how the fitting is done
#' 
#' @return a list with two components: vars_and_lags is the estimated "optimal" lags
#'     to use for each variable, and theta_hat is the estimated "optimal"
#'     kernel parameters to use for each combination of variable and lag
est_ssr_params_stepwise_crossval <- function(data, ssr_control) {
    all_vars_and_lags <- rbind.fill(lapply(ssr_control$kernel_components,
		function(kernel_component) {
			kernel_component$vars_and_lags
		}))
    
    ## initialize cross-validation process
    ## the variable selected in previous iteration
    selected_var_lag_ind <- NULL
    ## cross-validation estimate of loss associated with current estimates
    crossval_prediction_loss <- Inf
    ## initial parameters: no variables/lags selected, no kernel parameters
    vars_and_lags <- data.frame(var_name = NA,
		lag_value = NA,
		combined_name = NA)
    theta_hat <- vector("list", length(ssr_control$kernel_components))
    
    all_evaluated_models <- list()
    all_evaluated_model_descriptors <- list()
    
    repeat {
        ## get cross-validation estimates of performance for model obtained by
        ## adding or removing each variable/lag combination (except for the one
        ## updated in the previous iteration) from the model, as well as the
        ## corresponding parameter estimates
        
        ## commented out use of foreach for debugging purposes
#        crossval_results <- foreach(i=seq_len(nrow(all_vars_and_lags)),
#            .packages=c("ssr", ssr_control$par_packages),
#            .combine="c") %dopar% {
        crossval_results <- lapply(seq_len(nrow(all_vars_and_lags)),
			function(i) {
	        	descriptor_current_model <- update_vars_and_lags(vars_and_lags,
					all_vars_and_lags[i, "var_name"],
					all_vars_and_lags[i, "lag_value"])
        		
	        	model_i_previously_evaluated <- any(sapply(
					all_evaluated_model_descriptors,
	        		function(descriptor) {
	        			identical(descriptor_current_model, descriptor)
	        		}
				))
        		
	            if(!model_i_previously_evaluated) {
	            	potential_step_result <- 
	                    est_ssr_params_stepwise_crossval_one_potential_step(
	                        prev_vars_and_lags=vars_and_lags,
	                        prev_theta=theta_hat,
	                        update_var_name=all_vars_and_lags[i, "var_name"],
	                        update_lag_value=all_vars_and_lags[i, "lag_value"],
	                        data=data,
	                        ssr_control=ssr_control)
                
	                all_evaluated_models <-
	                	c(all_evaluated_models,
	                		list(potential_step_result))
	               	all_evaluated_model_descriptors <-
	               		c(all_evaluated_model_descriptors,
	               			potential_step_result["all_vars_and_lags"])
	                
	                return(potential_step_result)
#	               return(list(potential_step_result)) # put results in a list so that combine="c" is useful
	            } else {
	            	return(NULL)
	            }
#	        }
	        }
		)
       
        ## drop elements corresponding to previously explored models
        non_null_components <- sapply(crossval_results,
        	function(component) { !is.null(component) }
        )
        crossval_results <- crossval_results[non_null_components]
        
        ## pull out loss achieved by each model, find the best value
        loss_achieved <- sapply(crossval_results, function(component) {
            component$loss
        })
        optimal_loss_ind <- which.min(loss_achieved)
       
#        print("loss achieved is:")
#        print(loss_achieved)
#        print("\n")
        
        ## either update the model and keep going or stop the search
        if(loss_achieved[optimal_loss_ind] < crossval_prediction_loss) {
            ## found a model improvement -- update and continue
            selected_var_lag_ind <- optimal_loss_ind
            crossval_prediction_loss <- loss_achieved[selected_var_lag_ind]
			vars_and_lags <-
				crossval_results[[selected_var_lag_ind]]$vars_and_lags
            theta_hat <- crossval_results[[selected_var_lag_ind]]$theta
        } else {
            ## could not find a model improvement -- stop search
            break
        }
    }

    return(list(vars_and_lags=vars_and_lags,
        theta_hat=theta_hat))
}

#' Update the list that keeps track of which combinations of variables and lags
#' are included in the model by adding or removing (as appropriate) the
#' combination specified by update_var_name and update_lag_values
#'
#' @param prev_vars_and_lags list of previous variable/lag combinations to update
#' @param update_var_name name of the variable to update
#' @param update_lag_value value of the lag to update
#'
#' @return updated list of variable/lag combinations
update_vars_and_lags <- function(prev_vars_and_lags,
		update_var_name,
		update_lag_value) {
	updated_vars_and_lags <- prev_vars_and_lags
	
	existing_ind <- which(prev_vars_and_lags$var_name == update_var_name &
		prev_vars_and_lags$lag_value == update_lag_value)
	
    if(length(existing_ind) > 0) {
        ## remove variable/lag combination from model
		updated_vars_and_lags <- updated_vars_and_lags[-existing_ind, , drop = FALSE]
    } else {
        ## add variable/lag combination to model
		updated_vars_and_lags <- rbind(updated_vars_and_lags,
			data.frame(var_name = update_var_name,
				lag_value = update_lag_value,
				combined_name = paste0(update_var_name,
					"_lag",
					update_lag_value)))
    }
    
    return(updated_vars_and_lags)
}


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
}

#' Estimate the parameters theta and corresponding cross-validated estimate of
#' loss for one possible model obtained by adding or removing a variable/lag
#' combination from the model obtained in the previous iteration of the stepwise
#' search procedure.
#' 
#' @param prev_vars_and_lags list representing combinations of variables and lags
#'     included in the model obtained at the previous step
#' @param prev_theta list representing the kernel parameter estimates obtained
#'     at the previous step
#' @param update_var_name the name of the variable to try adding or removing
#'     from the model
#' @param update_lag_value the value of the lag for the variable specified by
#'     update_var_name to try adding or removing from the model
#' @param data the data frame with observations used in estimating model
#'     parameters
#' @param ssr_control a list of parameters specifying how the fitting is done
#' 
#' @return a list with three components: loss is a cross-validation estimate of
#'     the loss associated with the estimated parameter values for the given
#'     model, lags is a list representing combinations of variables and lags
#'     included in the updated model, and theta is a list representing the
#'     kernel parameter estimates in the updated model
est_ssr_params_stepwise_crossval_one_potential_step <- function(
		prev_vars_and_lags,
    	prev_theta,
    	update_var_name,
    	update_lag_value,
    	data,
    	ssr_control) {
    ## updated variable and lag combinations included in model
    updated_vars_and_lags <- update_vars_and_lags(prev_vars_and_lags,
		update_var_name,
		update_lag_value)
    
    ## initial values for theta; a list with two components:
    ##     theta_est, vector of initial values in vector form on scale appropriate for
    ##         estimation.
    ##     theta_fixed, list of values that will not be estimated, one component
    ##         for each component kernel function
    theta_init <- initialize_theta(prev_theta,
    	update_var_name,
    	update_lag_value,
    	data,
    	ssr_control)
    	
    theta_est_init <- extract_vectorized_theta_est_from_theta(theta_init,
    	ssr_control)
    
    ## optimize parameter values
    optim_result <- optim(par=theta_est_init,
        fn=ssr_crossval_estimate_parameter_loss,
#        gr = gradient_ssr_crossval_estimate_parameter_loss,
		gr=NULL,
		theta=theta_init,
        vars_and_lags=updated_vars_and_lags,
        data=data,
        ssr_control=ssr_control,
        method="L-BFGS-B",
        lower=-50,
        #		upper=10000,
        #control=list(),
        hessian=FALSE)
    
#    print(optim_result)
    
    ## convert back to list and original parameter scale
    updated_theta <- update_theta_from_vectorized_theta_est(optim_result$par,
        theta_init,
        ssr_control)
    
    return(list(
        loss=optim_result$value,
        vars_and_lags=updated_vars_and_lags,
        theta=updated_theta
    ))
}

#' Using cross-validation, estimate the loss associated with a particular set
#' of lags and kernel function parameters.
#' 
#' @param theta_est_vector vector of kernel function parameters that are being
#'     estimated
#' @param theta list of kernel function parameters, both those that are being
#'     estimated and those that are out of date.  Possibly the values of
#'     parameters being estimated are out of date; they will be replaced with
#'     the values in theta_est_vector.
#' @param vars_and_lags list representing combinations of variables and lags
#'     included in the model
#' @param data the data frame to use in performing cross validation
#' @param ssr_control a list of parameters specifying how the fitting is done
#' 
#' @return numeric -- cross-validation estimate of loss associated with the
#'     specified parameters
ssr_crossval_estimate_parameter_loss <- function(theta_est_vector,
		theta,
		vars_and_lags,
    	data,
    	ssr_control) {
	max_lag <- max(unlist(lapply(ssr_control$kernel_components,
		function(kernel_component) {
			kernel_component$vars_and_lags$lag_value
		}
	)))
    ## set up theta list containing both the kernel parameters that are being
    ## estimated and the kernel parameters that are being held fixed
    ## also, transform back to original scale
    ## convert back to list and original parameter scale
    theta <- update_theta_from_vectorized_theta_est(theta_est_vector,
        theta,
        ssr_control)
    
    ## create data frame of "examples" -- lagged observation vectors and
    ## corresponding prediction targets
    cross_validation_examples <- assemble_training_examples(data,
		vars_and_lags,
        ssr_control$y_names,
        leading_rows_to_drop=max_lag,
        additional_rows_to_drop=NULL,
        prediction_horizons=ssr_control$prediction_horizons,
        drop_trailing_rows=TRUE)
    
    ## This could be made more computationally efficient by computing
    ## kernel values for all relevant combinations of lags for each variable,
    ## then combining as appropriate -- currently, the same kernel value may be
    ## computed multiple times in the call to ssr_predict_given_lagged_obs
    crossval_loss_by_time_ind <- sapply(
        seq_len(nrow(cross_validation_examples$lagged_obs)),
        function(t_pred) {
            ## get training indices -- those indices not within
            ## t_pred +/- ssr_control$crossval_buffer
            t_train <- seq_len(nrow(cross_validation_examples$lagged_obs))
            t_train <- t_train[!(t_train %in%
                seq(from=t_pred - ssr_control$crossval_buffer,
                    to=t_pred + ssr_control$crossval_buffer))]
            
            ## calculate kernel weights and centers for prediction at
            ## prediction_lagged_obs based on train_lagged_obs and
            ## train_lead_obs
            ## assemble lagged and lead observations -- subsets of
            ## cross_validation_examples given by t_pred and t_train
            ## we can re-use the weights at different prediction_target_inds,
            ## and just have to adjust the kernel centers
            prediction_target_ind <- 1
            
            train_lagged_obs <- cross_validation_examples$lagged_obs[
                t_train, , drop=FALSE]
            train_lead_obs <- cross_validation_examples$lead_obs[
                t_train, prediction_target_ind, drop=FALSE]
            prediction_lagged_obs <- 
                cross_validation_examples$lagged_obs[
                    t_pred, , drop=FALSE]
            prediction_lead_obs <-
                cross_validation_examples$lead_obs[
                    t_pred, prediction_target_ind, drop=FALSE]
            
            ## for each prediction target variable, compute loss
            crossval_loss_by_prediction_target <- sapply(
                seq_len(ncol(cross_validation_examples$lead_obs)),
                function(prediction_target_ind) {
                    ## calculate and return value of loss function based on prediction
                    ## and realized value
                    loss_fn_args <- ssr_control$loss_fn_args
                    loss_fn_args$prediction_result <- ssr_predict_given_lagged_obs(
		                train_lagged_obs=train_lagged_obs,
		                train_lead_obs=train_lead_obs,
		                prediction_lagged_obs=prediction_lagged_obs,
		                prediction_lead_obs=prediction_lead_obs,
		                ssr_fit=list(theta_hat=theta,
		                    ssr_control=ssr_control
		                ),
		                prediction_type=ssr_control$loss_fn_prediction_type)
                    
                    loss_fn_args$obs <- as.numeric(
                        cross_validation_examples$lead_obs[
                            t_pred, prediction_target_ind]
                    )
                    
                    return(do.call(ssr_control$loss_fn_name, loss_fn_args))
                })
            
            return(sum(crossval_loss_by_prediction_target))
        })
    
#    browser()
    if(any(is.na(crossval_loss_by_time_ind))) {
        ## parameters resulted in numerical instability
        ## return largest non-infinite value
        return(.Machine$double.xmax)
    } else {
        return(sum(crossval_loss_by_time_ind))
    }
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
			ssr_control$kernel_components[[ind]]$vectorize_kernel_param_fns,
       		c(list(theta_list = theta[[ind]],
				vars_and_lags = vars_and_lags,
       			ssr_control = ssr_control),
       			ssr_control$kernel_components[[ind]]$vectorize_kernel_param_args)
       	))
    }

    return(theta)
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
				update_theta_from_vectorized_theta_est,
        	c(list(theta_est_vector = theta_est_vector[
					seq(from = theta_vector_component_start_ind,
        				to = length(theta_vector))],
        		theta = theta[[ind]],
        		ssr_control = ssr_control),
        		ssr_control$kernel_components[[ind]]$
					update_theta_from_vectorized_theta_est_args)
        )
        
        theta_vector_component_start_ind <- theta_vector_component_start_ind +
			temp$num_theta_vals_used
        
        theta[[ind]] <- temp$params
    }

    return(theta)
}






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
        normalize_weights=TRUE,
        prediction_type="distribution") {
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
        normalize_weights,
        prediction_type)
}

#' Construct data frame of lagged observation vectors and data frame of
#' corresponding prediction targets.
#' 
#' @param data data frame with observations of all variables at all times
#' @param vars_and_lags list representing combinations of variables and lags
#'     included in the lagged observation vectors
#' @param y_names names of variables included as prediction targets
#' @param leading_rows_to_drop an integer specifying the number of leading rows
#'     to drop.  This will typically be max(unlist(lags)), but could be larger
#'     for example if we are in the process of searching lags to determine
#'     which ones to include.  Then for example if our search is for lags up to
#'     52, we would want to set leading_rows_to_drop = 52.
#' @param additional_rows_to_drop an integer vector specifying indices of
#'     additional rows to drop.  For example, if we are performing
#'     cross-validation, we might want to drop all rows within +/- 52 indices of
#'     the current prediction target.
#' @param prediction_horizons an integer vector specifying the number of time
#'     steps between the last observation and the prediction target.
#' @param drop_trailing_rows boolean:  drop the last prediction_horizon rows?
#'     These are the rows for which we can form lagged observation vectors, but
#'     we cannot obtain a corresponding prediction target.
#'     
#' @return a list with two components: lagged_obs is a data frame with lagged
#'     observation vectors, and lead_obs is a data frame with corresponding
#'     prediction targets.
assemble_training_examples <- function(data,
		vars_and_lags,
    	y_names,
    	leading_rows_to_drop,
    	additional_rows_to_drop,
    	prediction_horizons,
    	drop_trailing_rows=TRUE) {
    ## which rows should not be used as regression/density estimation examples
    ## either because the corresponding regression example cannot be formed
    ## or because we are performing cross-validation and don't want to use
    ## times adjacent to the prediction target

    ## too early
    all_train_rows_to_drop <- seq_len(leading_rows_to_drop)

    ## passed in indices -- near prediction target
    all_train_rows_to_drop <- c(all_train_rows_to_drop, additional_rows_to_drop)
    
    ## too late
    max_prediction_horizon <- max(prediction_horizons)
    if(drop_trailing_rows) {
        all_train_rows_to_drop <- c(all_train_rows_to_drop,
            seq(from=nrow(data) - max_prediction_horizon + 1,
                to=nrow(data)))
    }
    
    ## compute lagged observation vectors for train and prediction data
    train_lagged_obs <- compute_lagged_obs_vecs(data,
		vars_and_lags,
        all_train_rows_to_drop)
    
    ## compute lead observation series for train data
    train_lead_obs <- data.frame(rep(NA, nrow(data)))
    for(y_name in y_names) {
        for(prediction_horizon in prediction_horizons) {
            train_lead_obs_inds <-
                seq(from=1, to=nrow(data)) + prediction_horizon
            component_name <- paste0(y_name, "_horizon", prediction_horizon)
            train_lead_obs[[component_name]] <-
                data[train_lead_obs_inds, y_name, drop=TRUE]
        }
    }
    train_lead_obs <- train_lead_obs[-all_train_rows_to_drop, , drop=FALSE]
    
    ## drop initial column of NA's
    train_lead_obs <- train_lead_obs[, -1, drop=FALSE]
    
    return(list(lagged_obs=train_lagged_obs,
        lead_obs=train_lead_obs))
}

#' Construct data frame of lagged observation vectors.
#' 
#' @param data data frame with observations of all variables at all times
#' @param vars_and_lags list representing combinations of variables and lags
#'     included in the lagged observation vectors
#' @param leading_rows_to_drop an integer specifying the number of leading rows
#'     to drop.  This will typically be max(unlist(lags)), but could be larger
#'     for example if we are in the process of searching lags to determine
#'     which ones to include.  Then for example if our search is for lags up to
#'     52, we would want to set leading_rows_to_drop = 52.
#'     
#' @return a list with one component: lagged_obs is a data frame with lagged
#'     observation vectors.
assemble_prediction_examples <- function(data,
		vars_and_lags,
    	y_names,
    	prediction_horizons,
    	leading_rows_to_drop) {
    all_prediction_rows_to_drop <- seq_len(leading_rows_to_drop) # too early

    prediction_lagged_obs <- compute_lagged_obs_vecs(data,
		vars_and_lags,
        all_prediction_rows_to_drop)
    
    prediction_lead_obs <- data.frame(rep(NA, nrow(data)))
    for(y_name in y_names) {
        for(prediction_horizon in prediction_horizons) {
            prediction_lead_obs_inds <-
                seq(from=1, to=nrow(data)) + prediction_horizon
            component_name <- paste0(y_name, "_horizon", prediction_horizon)
            prediction_lead_obs[[component_name]] <-
                data[prediction_lead_obs_inds, y_name, drop=TRUE]
        }
    }
    prediction_lead_obs <- prediction_lead_obs[-all_prediction_rows_to_drop, , drop=FALSE]
    
    return(list(lagged_obs=prediction_lagged_obs))
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
#' @param normalize_weights boolean, should the weights be normalized?
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
    normalize_weights=TRUE,
    prediction_type="distribution") {
    if(identical(prediction_type, "distribution")) {
    	return(ssr_dist_predict_given_lagged_lead_obs(train_lagged_obs,
		    train_lead_obs,
		    prediction_lagged_obs,
		    prediction_test_lead_obs,
		    ssr_fit))
    } else if(identical(prediction_type, "point")) {
    	return(ssr_point_predict_given_lagged_obs(train_lagged_obs,
		    train_lead_obs,
		    prediction_lagged_obs,
		    ssr_fit,
		    normalize_weights))
    } else if(identical(prediction_type, "centers-and-weights")) {
    	return(ssr_kernel_centers_and_weights_predict_given_lagged_obs(train_lagged_obs,
		    train_lead_obs,
		    prediction_lagged_obs,
		    ssr_fit,
		    normalize_weights))
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
#' @return a list with three components:
#'     log_weights: a vector of length = length(train_lagged_obs) with
#'         the log of weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     weights: a vector of length = length(train_lagged_obs) with
#'         the weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     centers: a copy of the train_lead_obs argument -- kernel centers
ssr_kernel_centers_and_weights_predict_given_lagged_obs <- function(train_lagged_obs,
    train_lead_obs,
    prediction_lagged_obs,
    ssr_fit,
    normalize_weights=TRUE) {
    
    ## compute log of kernel function values representing similarity
    ## of lagged observation vectors from train_data and predict_data.
    ## result is a vector of length nrow(train_lagged_obs)
    log_weights <- compute_kernel_values(train_lagged_obs,
        prediction_lagged_obs,
        ssr_fit$theta_hat,
        ssr_fit$ssr_control,
        log = TRUE)

    ## if requested, normalize log weights
    if(normalize_weights) {
        log_weights <- compute_normalized_log_weights(log_weights)
    }
    
    return(list(log_weights=log_weights,
        weights=exp(log_weights),
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
ssr_point_predict_given_lagged_obs <- function(train_lagged_obs,
    train_lead_obs,
    prediction_lagged_obs,
    ssr_fit,
    normalize_weights=TRUE) {
    
    ## compute log of kernel function values representing similarity
    ## of lagged observation vectors from train_data and predict_data.
    ## result is a vector of length nrow(train_lagged_obs)
    log_weights <- compute_kernel_values(train_lagged_obs,
        prediction_lagged_obs,
        ssr_fit$theta_hat,
        ssr_fit$ssr_control,
        log = TRUE)

    ## if requested, normalize log weights
    if(normalize_weights) {
        log_weights <- compute_normalized_log_weights(log_weights)
    }
    
    return(list(log_weights=log_weights,
        weights=exp(log_weights),
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
        ssr_fit$theta_hat,
        ssr_fit$ssr_control,
        log = TRUE)
    
    ## compute log(sum(kernel values x))
    log_sum_kernel_values_x <- logspace_sum(log_kernel_values_x)
    
    log_result <- sapply(seq_len(nrow(prediction_test_lead_obs)),
        function(prediction_test_row_ind) {
            log_kernel_values_xy <- compute_kernel_values(
                cbind(train_lagged_obs, train_lead_obs),
                cbind(prediction_lagged_obs,
                    prediction_test_lead_obs),
		        ssr_fit$theta_hat,
		        ssr_fit$ssr_control,
		        log = TRUE
            )
            
            return(logspace_sum(log_kernel_values_xy) -
                log_sum_kernel_values_x)
        })
    
    return(log_result)
}




#' Normalize a vector of weights on log scale: given a vector of log_weights
#' such that exp(log_weights) are proportional to the final weights, update
#' so that exp(log_weights) sums to 1.
#' 
#' @param log_weights: a vector of log(w) where w is proportional to the weights
#' 
#' @return normalized log_weights so that sum(exp(log_weights)) = 1
compute_normalized_log_weights <- function(log_weights) {
    ## normalize
    norm_const <- logspace_sum_matrix_rows(matrix(log_weights, nrow = 1))
    log_weights <- log_weights - norm_const
    
    ## normalize again -- if the initial log_weights were "extreme", the norm_const
    ## computed above may be approximate.
    norm_const <- logspace_sum_matrix_rows(matrix(log_weights, nrow = 1))
    
    return(log_weights - norm_const)
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
    theta,
    ssr_control,
    log = TRUE) {
	if(!(identical(nrow(prediction_obs), 1L))) {
		stop("In call to compute_kernel_values, prediction_obs must have exactly 1 row.")
	}

    ## create a matrix of log kernel values by component
	## rows correspond to time points in train_obs, columns to components of
	## the kernel function
    log_kernel_component_values <- matrix(0,
        nrow=nrow(train_obs),
        ncol=length(ssr_control$kernel_components))
    
    for(ind in seq_along(ssr_control$kernel_components)) {
		combined_names_in_component <-
			ssr_control$kernel_components[[ind]]$vars_and_lags$combined_name
        col_names <- colnames(train_obs)[
            colnames(train_obs) %in% combined_names_in_component]
        
        if(length(col_names) > 0) {
	        ## assemble arguments to kernel function
	        kernel_fn_args <- c(theta[[ind]],
	            ssr_control$kernel_components[[ind]]$theta_fixed)
	        kernel_fn_args$x <- train_obs[, col_names, drop = FALSE]
	        kernel_fn_args$center <- prediction_obs[, col_names, drop = FALSE]
	        kernel_fn_args$log <- TRUE
	        
	        ## call kernel function and store results
	        log_kernel_component_values[, ind] <-
	            do.call(ssr_control$kernel_components[[ind]]$kernel_fn,
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


#' Compute a data frame with lagged observation vectors
#' 
#' @param data a data frame
#' @param vars_and_lags a named list: The component name matches the name of 
#'     one of the variables in data, and the component value is an integer
#'     vector of lags to include for that variable
#' @param rows_to_drop an integer vector specifying rows to drop after computing
#'     lagged observation vectors.
compute_lagged_obs_vecs <- function(data,
		vars_and_lags,
    	rows_to_drop) {
    ## set or validate leading_rows_to_drop
    max_lag <- max(vars_and_lags$lag_value)
    if(any(rows_to_drop < 0 | rows_to_drop > nrow(data))) {
        stop("all entries of rows_to_drop must integers between 1 and nrow(data)")
    } else if(!all(seq_len(max_lag) %in% rows_to_drop)) {
        stop("all integers between 1 and the maximum entry in the lags argument must be contained in rows_to_drop")
    }
    
    ## create a data frame with one column for each entry in the lags argument
    result <- as.data.frame(matrix(NA,
        nrow=nrow(data),
        ncol=nrow(vars_and_lags)))
    colnames(result) <- vars_and_lags$combined_name
    
    ## set column values in result
	for(new_var_ind in seq_len(nrow(vars_and_lags))) {
		lag_val <- vars_and_lags[new_var_ind, "lag_value"]
		var_name <- vars_and_lags[new_var_ind, "var_name"]
		combined_name <- vars_and_lags[new_var_ind, "combined_name"]
		
		result_inds <- seq(from=lag_val + 1, to=nrow(result))
		data_inds <- seq(from=1, to=nrow(result) - lag_val)
		
		result[result_inds, combined_name] <- data[data_inds, var_name]
	}
	
    ## drop specified rows
    if(length(rows_to_drop) > 0) {
        result <- result[-rows_to_drop, , drop=FALSE]
    }
    
    return(result)
}



#' This is a wrapper for the ssr_predict function to perform prediction from
#' the dengue data sets for the competition.  Computes KDE predictions for
#' the number of cases in the weeks indexed by t + prediction_horizon, where
#'   - t is specified by last_obs_season and last_obs_week
#'   - prediction_horizon varies over the values in prediction_horizons
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
#' @param lag_max is the number of lags to use in forming the state space
#'     representation via lagged observations
#' @param prediction_horizons is a vector giving the number of steps ahead
#'     to do prediction from the week before first_predict_week
#' @param data is the data set to use
#' @param season_var is a string with the name of the variable in the data data
#'     frame with the season.  I'm sure there is a better way to do this...
#' @param week_var is a string with the name of the variable in the data data
#'     frame with the week in the season.  I'm sure there is a better way to do this...
#' @param predict_vars is a string vector with the names of the variables in the data data
#'     frame to use for SSR.  I'm sure there is a better way to do this...
#' @param prediction_types is a character vector indicating the prediction types
#'     to perform: may contain one or more of "pt" and "density"
#' 
#' @return list of data frames with predictions.  one list component per entry
#'     in prediction_types argument.  Density estimates are from KDE along a
#'     grid of values for total_counts for each week that we predicted at.
ssr_predict_dengue_one_week <- function(last_obs_season,
    last_obs_week,
	vars_and_lags,
    theta,
    prediction_bw,
    prediction_horizon,
    data,
    season_var,
    week_var,
    X_names,
    y_names,
    time_name,
    ssr_fit,
    prediction_types = c("pt", "density")) {
    
  	if(length(y_names) > 1) {
  		stop("SSR Prediction for multiple variables is not yet implemented!")
  	}
  	
  	if(missing(ssr_fit)) {
  		ssr_fit <- list(ssr_control=create_ssr_control_default(X_names, y_names, time_name, data),
            X_names=X_names,
            y_names=y_names,
            time_name=time_name,
			vars_and_lags=vars_and_lags,
            theta_hat=theta,
            train_data=data)
  	}
	
    ## get indices for prediction data
    ## in order to predict at time t + prediction_horizon,
    ##  - data must start at time t - lag_max so that we can get lag_max + 1
    ##    values in forming the lagged observation vector
    ##  - we need the number of data points equal to lag_max + 1
    ## first_predict_ind is the first index in the data data frame
    ## that needs to be included in the prediction data
    lag_max <- max(vars_and_lags$lag_value)
    last_obs_ind <- which(data[[season_var]] == last_obs_season &
        data[[week_var]] == last_obs_week)
    first_predict_data_ind <- last_obs_ind - lag_max
    predict_data_inds <- seq(from=first_predict_data_ind, length=lag_max + 1)
    
    ## for training inds, drop anything within +/- 1 year of the last
    ## observation
    time_points_per_year <- 52
    predict_seasons <- data[[season_var]][last_obs_ind + prediction_horizon]
    additional_training_rows_to_drop <-
        seq(from=max(1, last_obs_ind - time_points_per_year),
            to=min(nrow(data), last_obs_ind + time_points_per_year))
    
    ## do prediction
    ssr_predictions <- ssr_predict(
        ssr_fit=ssr_fit,
        prediction_data=data[predict_data_inds, ],
        leading_rows_to_drop=lag_max,
        additional_training_rows_to_drop=additional_training_rows_to_drop,
        prediction_horizon=prediction_horizon,
        normalize_weights=TRUE,
        prediction_type="centers-and-weights")
    
    ## get point estimate or density predictions
    result <- list()
    
    if("centers-and-weights" %in% prediction_types) {
    	result$pred_centers_and_weights <- ssr_predictions
    }
    
    if("pt" %in% prediction_types) {
        result$pt_preds <- get_pt_predictions_one_week(ssr_predictions)
    }
    
    if("density" %in% prediction_types) {
        result$dist_preds <- get_dist_predictions_one_week(ssr_predictions)
    }
    
    return(result)
}

## a function to get point predictions for one week
get_pt_predictions_one_week <- function(ssr_predictions) {
    # get weighted mean
    pt_est <- weighted.mean(ssr_predictions$centers[, 1],
        ssr_predictions$weights)

    return(pt_est)
}

## a function to get kde predictions for one week
get_dist_predictions_one_week <- function(ssr_predictions, bw) {
    ## do weighted kde
    kde_est <- density(ssr_predictions$centers[, 1],
        weights = ssr_predictions$weights,
        bw = bw)

    return(data.frame(x = kde_est$x,
        est_density = kde_est$y))
}


#' Get week and season when we're predicting based on the index of the last obs
#' and the number of steps forward we're predicting
#' 
#' @param last_obs_ind index in the data data frame of the last observation
#'     before we start predicting
#' @param prediction_horizon the number of steps forward we're predicting
#' @param data the data frame, with columns named season_week and season
#' @param season_var is a string with the name of the variable in the data data
#'     frame with the season.  I'm sure there is a better way to do this...
#' @param week_var is a string with the name of the variable in the data data
#'     frame with the week in the season.  I'm sure there is a better way to do this...
#' 
#' @return a list with two components indicating the time of prediction:
#'     1) week is an integer from 1 to 52 and
#'     2) season is a string of the form "2008/2009"
get_prediction_season_week <- function(last_obs_ind, prediction_horizon, data, season_var, week_var) {
    wk <- (data[[week_var]][last_obs_ind] + prediction_horizon) %% 52
    if(wk == 0) {
        wk <- 52
    }
    seasons_advanced <- (data[[week_var]][last_obs_ind] + prediction_horizon
            - wk) / 52
    start_season_last_obs <- as.integer(substr(
            as.character(data[[season_var]][last_obs_ind]),
            start = 1,
            stop = 4))
	season <- start_season_last_obs + seasons_advanced
	
	if(nchar(as.character(data[[season_var]][last_obs_ind])) > 4) {
		season <- paste0(season, "/", season + 1)
	}
    
    return(list(week = wk, season = season))
}


### loss functions

#' Compute MASE from predictions that are in the form of kernel weights and
#' centers.  This is currently broken because we need predictions for a whole
#' time series to compute MASE, and the current "implementation" only takes
#' values for one time point....
mase_from_kernel_weights_and_centers <- function(
    kernel_weights_and_centers,
    obs) {
    stop("mase_from_kernel_weights_and_centers is broken -- do not use")
    pred <- get_pt_predictions_one_week(kernel_weights_and_centers)
    mase(obs, pred)
}

#' Compute "mean" absolute error for one time point from prediction in the form
#' of kernel weights and centers.
#' 
#' @param kernel_weights_and_centers a named list with two components: weights
#'     is a vector of kernel weights, and centers is a vector of kernel centers
#' @param obs is a numeric with length 1 containign the observed value for one
#'     time point
#' 
#' @return abs(obs - prediction) where prediction is the weighted mean of the
#'     kernel centers.
mae_from_kernel_weights_and_centers <- function(
    kernel_weights_and_centers,
    obs) {
    pred <- get_pt_predictions_one_week(kernel_weights_and_centers)
    mae(obs, pred)
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
