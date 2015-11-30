library(lubridate)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape)
library(ssr)


location <- "sanjuan"
location_for_ssr_fit_file <- "San_Juan"

ssr_fits_by_ph <- lapply(seq(from=1, to=52, by=1),
#ssr_fits_by_ph <- lapply(seq(from=1, to=2, by=1),
    function(phl) {
        file_name <- file.path("/media", "evan", "data", "Reich", "dengue-ssr-prediction", "inst", "intermediate-results", "ssr-poster", "poster-fits",
            paste0("updated-fit-competition-ssr-ph",
            phl,
            "-",
            location_for_ssr_fit_file,
            ".Rdata"))
        read_env <- new.env()
        load(file_name, envir=read_env)
        ssr_fit <- read_env$ssr_fit
        
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
        
        return(ssr_fit)
    })


predictions_df <- data.frame(ph=rep(seq_len(52), times = 5 * 52),
	last_obs_season=rep(c("2008/2009", "2009/2010", "2010/2011", "2011/2012", "2012/2013"), each = 52, times = 52),
	last_obs_week=rep(seq_len(52) - 1, each = 52 * 5),
	model="sarima",
	stringsAsFactors=FALSE)
predictions_df$prediction_season <- predictions_df$last_obs_season
predictions_df$prediction_week <- predictions_df$last_obs_week + predictions_df$ph

inds_last_obs_season_prev_year <- which(predictions_df$last_obs_week == 0)
predictions_df$last_obs_season[inds_last_obs_season_prev_year] <- 
	sapply(predictions_df$last_obs_season[inds_last_obs_season_prev_year],
		function(next_season) {
			start_year <- as.integer(substr(next_season, 1, 4)) - 1L
			paste0(start_year, "/", start_year + 1)
		}
	)
predictions_df$last_obs_week[inds_last_obs_season_prev_year] <- 52L

inds_prediction_season_next_year <- which(predictions_df$prediction_week > 52)
predictions_df$prediction_season[inds_prediction_season_next_year] <- 
	sapply(predictions_df$prediction_season[inds_prediction_season_next_year],
		function(next_season) {
			start_year <- as.integer(substr(next_season, 1, 4)) + 1L
			paste0(start_year, "/", start_year + 1)
		}
	)
predictions_df$prediction_week[inds_prediction_season_next_year] <-
	predictions_df$prediction_week[inds_prediction_season_next_year] - 52L
	
predictions_df$log_score <- NA
predictions_df$prediction <- NA
predictions_df$prediction_median <- NA
predictions_df$AE <- NA
predictions_df$predictive_50pct_lb <- NA
predictions_df$predictive_50pct_ub <- NA
predictions_df$predictive_90pct_lb <- NA
predictions_df$predictive_90pct_ub <- NA
predictions_df$week_start_date <- San_Juan_test$week_start_date[1]

q_weighted_mixt_norm <- function(q, weights, centers, bw) {
	sapply(q, function(q_ind) {
		temp <- optim(par = sum(weights * centers),
	        fn = function(cutpt) {
	        	prob_lt_cutpt <- sum(weights * sapply(seq_along(centers), function(ind) {
	        		pnorm(cutpt, mean = centers[ind], sd = bw)
	        	}))
	        	return((q_ind - prob_lt_cutpt)^2)
	        },
	        gr = NULL,
	        method = "L-BFGS-B",
#	        method = "Brent",
#			lower = 0,
#	        upper = 100000,
	        control = list(),
	        hessian = FALSE)
	    
	    return(temp$par)
	})
}

q_weighted_mixt_lnorm <- function(q, weights, centers, bw) {
	sapply(q, function(q_ind) {
		temp <- optim(par = sum(weights * centers),
	        fn = function(cutpt) {
	        	prob_lt_cutpt <- sum(weights * sapply(seq_along(centers), function(ind) {
	        		pnorm(cutpt, mean = centers[ind], sd = bw)
	        	}))
	        	return((q_ind - prob_lt_cutpt)^2)
	        },
	        gr = NULL,
	        method = "L-BFGS-B",
#	        method = "Brent",
#			lower = 0,
#	        upper = 100000,
	        control = list(),
	        hessian = FALSE)
	    
	    return(exp(temp$par))
	})
}


for(predictions_df_row_ind in seq_len(nrow(predictions_df))) {
	last_obs_season <- predictions_df[predictions_df_row_ind, "last_obs_season"]
	last_obs_week <- predictions_df$last_obs_week[predictions_df_row_ind]
	predictions_df$last_obs_ind[predictions_df_row_ind] <- which(San_Juan_test$season == last_obs_season &
		San_Juan_test$season_week == last_obs_week)
}
predictions_df <- predictions_df[predictions_df$last_obs_ind + predictions_df$ph <= nrow(San_Juan_test), ]

expected_value_exp_norm_mixt <- function(log_weights, centers, sds) {
	if(length(sds) == 1) {
		sds <- rep(sds, length(log_weights))
	}
	
	integrand <- function(x, log_weights, centers, sds) {
#		browser()
		return(sapply(x, function(x_ind) {
			z <- log(x_ind + 1)
			norm_mixt_component <- logspace_sum(sapply(seq_along(log_weights), function(ind) { log_weights[ind] + dnorm(x_ind, centers[ind], sds[ind], log = TRUE) }))
			return(exp(norm_mixt_component + log(x_ind) - log(x_ind + 1)))
		}))
	}
	
	return(integrate(integrand, lower = 0, upper = Inf, log_weights = log_weights, centers = centers, sds = sds, subdivisions = 10000000L)$value)
}

for(predictions_df_row_ind in seq_len(nrow(predictions_df))) {
	ph <- as.numeric(predictions_df$ph[predictions_df_row_ind])
	last_obs_season <- predictions_df[predictions_df_row_ind, "last_obs_season"]
	last_obs_week <- predictions_df$last_obs_week[predictions_df_row_ind]
	last_obs_ind <- which(San_Juan_test$season == last_obs_season &
		San_Juan_test$season_week == last_obs_week)
	
	predictions_df$week_start_date[predictions_df_row_ind] <- San_Juan_test$week_start_date[last_obs_ind + as.numeric(ph)]
	
	ssr_fit <- ssr_fits_by_ph[[ph]]
	
	new_data <- San_Juan_test[seq_len(last_obs_ind + ph), , drop=FALSE]
	## add log column
	new_data$log_total_cases <- log(new_data$total_cases + 1)

	## add smooth log column
	sm <- loess(log_total_cases ~ as.numeric(week_start_date), data=new_data, span=12 / nrow(new_data))
	new_data$smooth_log_cases <- sm$fitted

	## add time column
	new_data$time_ind <- seq_len(nrow(new_data))
	
	ssr_fit$train_data <- new_data
	
	ssr_prediction <- ssr_predict_dengue_one_week(last_obs_season = last_obs_season,
		last_obs_week = last_obs_week,
        lags = ssr_fit$lags_hat,
        theta = ssr_fit$theta_hat,
        prediction_bw = ssr_fit$predictive_bw,
        prediction_horizon = ssr_fit$ssr_control$prediction_horizons,
        data = new_data,
        season_var = "season",
        week_var = "season_week",
        X_names = ssr_fit$ssr_control$X_names,
        y_names = ssr_fit$ssr_control$y_names,
        time_name = NULL,
        ssr_fit = ssr_fit,
	    prediction_types = "centers-and-weights")
	
	obs <- new_data[last_obs_ind + ph, ssr_fit$y_names]
	log_obs <- log(obs + 1)
	log_centers <- log(ssr_prediction$pred_centers_and_weights$centers[, 1] + 1)
	
	
	
	predictions_df$prediction_median[predictions_df_row_ind] <- exp(sum(ssr_prediction$pred_centers_and_weights$weights * log_centers)) - 1
	predictions_df$prediction[predictions_df_row_ind] <- expected_value_exp_norm_mixt(ssr_prediction$pred_centers_and_weights$log_weights, log_centers, ssr_fit$predictive_bw)
	predictions_df$AE[predictions_df_row_ind] <- abs(predictions_df$prediction[predictions_df_row_ind] - obs)

#	predictions_df$log_score[predictions_df_row_ind] <- logspace_sum(ssr_prediction$pred_centers_and_weights$log_weights + 
#		sapply(ssr_prediction$pred_centers_and_weights$centers, function(center) {dnorm(obs, mean = center, sd = ssr_fit$predictive_bw, log=TRUE)}))
	predictions_df$log_score[predictions_df_row_ind] <- logspace_sum(ssr_prediction$pred_centers_and_weights$log_weights + 
		sapply(log_centers, function(log_center) {dnorm(log_obs, mean = log_center, sd = ssr_fit$predictive_bw, log=TRUE)})) -
		log_obs

	temp <- q_weighted_mixt_lnorm(c(0.05, 0.25, 0.75, 0.95),
		weights = ssr_prediction$pred_centers_and_weights$weights,
		centers = log_centers,
		bw = ssr_fit$predictive_bw)
	predictions_df[predictions_df_row_ind, c("predictive_90pct_lb", "predictive_50pct_lb", "predictive_50pct_ub", "predictive_90pct_ub")] <-
		temp - 1
}








predictions_df <- predictions_df[predictions_df$week_start_date %in% San_Juan_test$week_start_date[San_Juan_test$season %in% c("2009/2010", "2010/2011", "2011/2012", "2012/2013")], ]

predictions_df$ph <- as.factor(predictions_df$ph)

ssr_predictions_df <- predictions_df
save(ssr_predictions_df, file = "/media/evan/data/Reich/dengue-ssr-prediction/inst/intermediate-results/ssr-poster/predictions/ssr-predictions.Rdata")

library(ggplot2)
ggplot() +
	geom_line(aes(x = week_start_date, y = total_cases), data = San_Juan_test) +
	geom_line(aes(x = week_start_date, y = prediction, colour = ph), data = predictions_df[predictions_df$ph %in% c(1, 13, 26, 39, 52), , drop = FALSE]) +
	theme_bw()


