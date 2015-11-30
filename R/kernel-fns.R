#' Squared exponential kernel function
#' 
#' @param x a vector of values at which to evaluate the kernel function
#' @param center a real number, center point for the kernel function
#' @param bw kernel bandwidth
#' @param return log scale value?
#' 
#' @return vector of the kernel function value at each point in x
squared_exp_kernel <- function(x, center, bw, log) {
    if(!is.numeric(center) | length(center) != 1) {
        stop("center must be numeric with length 1")
    }
    
    result <- -0.5 * ((x - center) / bw)^2
    
    if(log) {
        return(result)
    } else {
        return(exp(result))
    }
}


#' Periodic kernel function
#' 
#' @param x a vector of values at which to evaluate the kernel function
#' @param center a real number, center point for the kernel function
#' @param period kernel period
#' @param bw kernel bandwidth
#' 
#' @return vector of the kernel function value at each point in x
periodic_kernel <- function(x, center, period, bw, log) {
    result <- -0.5 * (sin(period * (x - center)) / bw)^2
    
    if(log) {
        return(result)
    } else {
        return(exp(result))
    }
}

#' Compute the parameters bw, bw_continuous, conditional_bw_discrete, and
#' conditional_center_discrete_offset_multiplier for the pdtmvn_kernel function
#' from the eigen-decomposition of the bandwidth matrix.
#' 
#' @param bw_evecs matrix whose columns are the eigenvectors of the bandwidth
#'     matrix.
#' @param bw_evals vector of eigenvalues of the bandwidth matrix.
#' @param continuous_vars Vector containing column indices for continuous
#'     variables.
#' @param discrete_vars Vector containing column indices for discrete variables.
#' 
#' @return Named list with four components: bw, bw_continuous,
#'     conditional_bw_discrete, and conditional_center_discrete_offset_multiplier
compute_pdtmvn_kernel_bw_params_from_bw_eigen <- function(bw_evecs,
    bw_evals,
    continuous_vars,
    discrete_vars) {
    
    bw <- bw_evecs %*% diag(bw_evals) %*% t(bw_evecs)
    
    bw_params <- c(list(bw = sigma),
        compute_sigma_subcomponents(sigma = bw,
            continuous_vars = continuous_vars,
            discrete_vars = discrete_vars,
            validate_level = 0)
    )
    
    names(bw_params)[names(bw_params) == "sigma_continuous"] <- "bw_continuous"
    names(bw_params)[names(bw_params) == "conditional_sigma_discrete"] <- 
        "conditional_bw_discrete"
    names(bw_params)[names(bw_params) == "conditional_mean_discrete_offset_multiplier"] <- 
        "conditional_center_discrete_offset_multiplier"
    
    return(bw_params)
}

#' Evaluate the kernel function given by the pdtmvn distribution.
#' 
#' @param a matrix of values at which to evaluate the kernel function, with
#'     column names specified.  Each row is an observation, each column is an
#'     observed variable.
#' @param center a real vector, center point for the kernel function
#' @param bw bandwidth matrix
#' @param bw_continuous the portion of bw corresponding to continuous variables
#' @param conditional_bw_discrete the Schur complement of the portion of bw
#'     corresponding to discrete variables
#' @param conditional_center_discrete_offset_multiplier Sigma_dc Sigma_c^{-1}.
#'     This is used in computing the mean of the underlying multivariate normal
#'     distribution for the discrete variables conditioning on the continuous
#'     variables.
#' @param continuous_vars character vector with names of variables that are to
#'     be treated as continuous.  May contain entries that do not appear in
#'     colnames(x).
#' @param discrete_vars character vector with names of variables that are to be
#'     treated as discrete.  May contain entries that do not appear in
#'     colnames(x)
#' @param discrete_var_range_functions a list with one entry for each element
#'     of discrete_vars.  Each entry is a named list of length 2; the element
#'     named "a" is a character string with the name of a function that returns
#'     a(x) for any real x, and the element named "b" is a character string with
#'     the name of a function that returns b(x) for any real x.
#' @param lower Vector of lower truncation points
#' @param upper Vector of upper truncation points
#' @param log logical; if TRUE, return the log of the kernel function value
#' 
#' @return the value of the kernel function given by the pdtmvn distribution
#'     at x.
pdtmvn_kernel <- function(x,
    center,
    bw,
    bw_continuous,
    conditional_bw_discrete,
    conditional_center_discrete_offset_multiplier,
    continuous_vars,
    discrete_vars,
    discrete_var_range_functions,
    lower,
    upper,
    log) {
    vars_used <- colnames(x)
    if(is.null(vars_used)) {
        stop("x must have column names")
    }
    continuous_vars_used <- vars_used[vars_used %in% continuous_vars]
    discrete_vars_used <- vars_used[vars_used %in% discrete_vars]
    
    return(dpdtmvn(x = x,
        mean = center,
        sigma = bw,
        sigma_continuous = bw_continuous,
        conditional_sigma_discrete = conditional_bw_discrete,
        conditional_mean_discrete_offset_multiplier = 
            conditional_center_discrete_offset_multiplier,
        lower = lower,
        upper = upper,
        continuous_vars = continuous_vars_used,
        discrete_vars = discrete_vars_used,
        discrete_var_range_functions = discrete_var_range_functions,
        log = TRUE,
        validate_level = 1))
}

#' A function to vectorize the parameters of the pdtmvn_kernel and convert
#' to estimation scale.
#' @param theta_list parameters for the pdtmvn_kernel in list format
#' @param parameterization character; currently, only supported value is
#'     "bw-diagonalized-est-eigenvalues"
#' @param ssr_control list of control parameters for ssr
#'
#' @return vector containing parameters that are estimated on a scale
#'     suitable for numerical optimization
vectorize_params_pdtmvn_kernel <- function(theta_list, parameterization, ssr_control) {
	if(identical(parameterization, "bw-diagonalized-est-eigenvalues")) {
		return(theta_list$log_bw_evals)
	} else {
		stop("Invalid parameterization for pdtmvn kernel function")
	}
}

#' A function to unvectorize the parameters of the pdtmvn_kernel and convert
#' from estimation scale to scale actually used.
#' 
#' @param theta_vector a numeric vector of parameters that are being optimized,
#'     on a scale suitable for use in optim.
#' @param parameterization character; currently, only supported value is
#'     "bw-diagonalized-est-eigenvalues"
#' @param bw_evecs square matrix with eigenvectors of the bandwidth matrix in
#'     columns
#' @param 
unvectorize_params_pdtmvn_kernel <- function(theta_vector, parameterization, bw_evecs, continuous_vars, discrete_vars, ssr_control) {
	if(identical(parameterization, "bw-diagonalized-est-eigenvalues")) {
		num_bw_evals <- ncol(bw_evecs)
		return(list(
			params = c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs,
			    	bw_evals = exp(theta_vector[seq_len(num_bw_evals)]),
			    	continuous_vars,
			    	discrete_vars),
			    list(log_bw_evals = theta_vector[seq_len(num_bw_evals)])),
			num_theta_vals_used = num_bw_evals
		))
	} else {
		stop("Invalid parameterization for pdtmvn kernel function")
	}
}


initialize_params_pdtmvn_kernel <- function(parameterization,
	x,
	update_var_name,
	update_lag_value,
	prev_theta,
	continuous_vars,
	discrete_vars,
	ssr_control) {
	require("robust")
	
	if(identical(parameterization, "bw-diagonalized-est-eigenvalues")) {
		sample_cov_hat <- covRob(x)$cov
		sample_cov_eigen <- eigen(sample_cov_hat)
		num_bw_evals <- ncol(bw_evecs)
		return(c(
			compute_pdtmvn_kernel_bw_params_from_bw_eigen(sample_cov_eigen$vectors,
			    bw_evals = sample_cov_eigen$values,
			    continuous_vars,
			    discrete_vars),
			list(
				bw_evecs = sample_cov_eigen$vectors,
				bw_evals = sample_cov_eigen$values,
				log_bw_evals = log(sample_cov_eigen$values)),
			ssr_control
			
		))
	} else {
		stop("Invalid parameterization for pdtmvn kernel function")
	}
}
