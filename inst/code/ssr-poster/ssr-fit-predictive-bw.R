###
### In this file, we use leave-one-chunk-of-time-out cross-validation to
### estimate the bandwidths associated with the prediction target.  We
### use the previously estimated bandwidths for variables that are
### conditioned on, which were obtained by minimizing cross-validated
### estimates of MAE for the point estimates.  These previous results
### were obtained by code in dengue-ssr-prediction/inst/ssr-competition/
###

library(lubridate)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape)
library(ssr)

crossval_estimate_prediction_target_bandwidth <- function(ssr_fit, prediction_horizon, season_var, week_var) {
#	prediction_horizon <- 4L
#	season_var <- "season"
#	week_var <- "season_week"
#	
#	debug(crossval_est_ssr_predictive_dist_log_loss)
#	
#	crossval_est_ssr_predictive_dist_log_loss(bw = 1,
#        lags = ssr_fit$lags_hat,
#        theta = ssr_fit$theta_hat,
#        prediction_horizon = ssr_fit$ssr_control$prediction_horizons,
#        data = ssr_fit$train_data,
#        season_var = season_var,
#        week_var = week_var,
#        X_names = ssr_fit$ssr_control$X_names,
#        y_names = ssr_fit$ssr_control$y_names,
#        time_name = NULL,
#        ssr_fit = ssr_fit
#		)
	
	prediction_target_bw_est <- optim(par = 1,
        fn = crossval_est_ssr_predictive_dist_neg_log_loss,
        gr = NULL,
        lags = ssr_fit$lags_hat,
        theta = ssr_fit$theta_hat,
        prediction_horizon = ssr_fit$ssr_control$prediction_horizons,
        data = ssr_fit$train_data,
        season_var = season_var,
        week_var = week_var,
        X_names = ssr_fit$ssr_control$X_names,
        y_names = ssr_fit$ssr_control$y_names,
        time_name = NULL,
        ssr_fit = ssr_fit,
        method = "L-BFGS-B",
#        method = "Brent",
#		lower = 0,
#        upper = 100000,
        control = list(),
        hessian = FALSE)
    
    ssr_fit$predictive_bw <- exp(prediction_target_bw_est$par)
    
    return(ssr_fit)
}

crossval_est_ssr_predictive_dist_neg_log_loss <- function(log_bw, lags, theta, prediction_horizon, data, season_var, week_var, X_names, y_names, time_name, ssr_fit) {
	bw <- exp(log_bw)
	last_obs_season_possibilities <- sapply(1990:2008, function(first_year) { paste0(first_year, "/", first_year+1) })
	last_obs_week_possibilities <- seq_len(52)
	season_week_combos <- data.frame(last_obs_season = rep(last_obs_season_possibilities, each = 52),
		last_obs_week = rep(last_obs_week_possibilities, times = length(last_obs_season_possibilities)))
	season_week_combos <- season_week_combos[!(season_week_combos$last_obs_season == "2008/2009" &
		season_week_combos$last_obs_week > 52 - prediction_horizon), ]
	
	# drop first indices in cases when lag 1 variables were used.
	max_lag <- max(unlist(lags))
	if(max_lag > 0) {
		season_week_combos <- season_week_combos[-seq_len(max_lag), , drop=FALSE]
	}
	
	## the next few lines adjust for the fact that I got halfway through rewriting the ssr package before
	## producing results for the poster.....
	new_theta <- list()
	new_theta_fixed <- list()
	new_kernel_fns <- list()
	new_kernel_variable_groups <- list()
	for(var_name in names(ssr_fit$lags_hat)) {
		for(lag_val in ssr_fit$lags_hat[[var_name]]) {
#			if(lag_val == 0) {
#				full_var_name <- var_name
#			} else {
				full_var_name <- paste0(var_name, "_lag", lag_val)
#			}
#			new_theta[[full_var_name]] <- c(ssr_fit$theta_est[full_var_name], ssr_fit$ssr_control$theta_fixed[full_var_name])
			new_theta[[full_var_name]] <- c(ssr_fit$theta_hat[full_var_name], ssr_fit$ssr_control$theta_fixed[var_name])
			if(var_name == "time_ind") {
				new_kernel_fns[[full_var_name]] <- "periodic_kernel"
			} else {
				new_kernel_fns[[full_var_name]] <- "squared_exp_kernel"
			}
			if(var_name %in% names(ssr_fit$ssr_control$theta_fixed)) {
				new_theta_fixed[[full_var_name]] <- ssr_fit$ssr_control$theta_fixed[[var_name]]
			} else {
				new_theta_fixed[[full_var_name]] <- list()
			}
			new_kernel_variable_groups[[full_var_name]] <- full_var_name
		}
	}
	ssr_fit$ssr_control$theta <- new_theta
	ssr_fit$ssr_control$kernel_fns <- new_kernel_fns
	ssr_fit$ssr_control$theta_fixed <- new_theta_fixed
	ssr_fit$ssr_control$kernel_variable_groups <- new_kernel_variable_groups
	
	results_by_week <- sapply(seq_len(nrow(season_week_combos)), function(row_ind) {
		last_obs_season <- season_week_combos$last_obs_season[row_ind]
		last_obs_week <- season_week_combos$last_obs_week[row_ind]
	
		ssr_prediction <- ssr_predict_dengue_one_week(last_obs_season = last_obs_season,
		    last_obs_week = last_obs_week,
		    lags = lags,
		    theta = theta,
		    prediction_bw = bw,
		    prediction_horizon = prediction_horizon,
		    data = data,
		    season_var = season_var,
		    week_var = week_var,
		    X_names = X_names,
		    y_names = y_names,
		    time_name = time_name,
		    ssr_fit = ssr_fit,
		    prediction_types = "centers-and-weights")
	
		last_obs_ind <- which(ssr_fit$train_data$season == last_obs_season &
			ssr_fit$train_data$season_week == last_obs_week)
		obs <- ssr_fit$train_data[last_obs_ind + prediction_horizon, y_names]
		return(logspace_sum(ssr_prediction$pred_centers_and_weights$log_weights + 
			sapply(ssr_prediction$pred_centers_and_weights$centers, function(center) {dnorm(log(obs + 1), mean = log(center + 1), sd = bw, log=TRUE)})) -
			log(obs + 1))
	})
#	browser()
	
	return(-1 * sum(results_by_week))
}



location <- "sanjuan"
location_for_ssr_fit_file <- "San_Juan"


junk <- lapply(seq(from=1, to=52, by=1),
    function(phl) {
        file_name <- file.path("/media", "evan", "data", "Reich", "dengue-ssr-prediction", "inst", "intermediate-results", "ssr-poster", "poster-fits",
            paste0("fit-poster-ssr-ph",
            phl,
            "-",
            location_for_ssr_fit_file,
            ".Rdata"))
        read_env <- new.env()
        load(file_name, envir=read_env)
        ssr_fit <- read_env$ssr_fit
        
        ssr_fit <- crossval_estimate_prediction_target_bandwidth(ssr_fit, prediction_horizon = phl, season_var = "season", week_var = "season_week")
        
        save(ssr_fit, file = file.path("/media", "evan", "data", "Reich", "dengue-ssr-prediction", "inst", "intermediate-results", "ssr-poster", "poster-fits",
            paste0("updated-fit-poster-ssr-ph",
            phl,
            "-",
            location_for_ssr_fit_file,
            ".Rdata")))
        
        return(NULL)
    })

names(ssr_fits_by_prediction_horizon_limit) <-
    paste0("phl", seq(from=4, to=52, by=4))

## assemble data set
data <- San_Juan_train

ssr_predict_dengue_one_week <- function(last_obs_season,
    last_obs_week,
    lags,
    theta,
    prediction_bw,
    prediction_horizon,
    data,
    season_var,
    week_var,
    X_names,
    y_names,
    time_name,
    prediction_types = c("pt", "density"))

