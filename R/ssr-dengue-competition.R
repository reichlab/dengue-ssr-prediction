## Functions that are specific to the dengue competition
## It's likely that none of these work anymore
## 
## ssr_predict_dengue_one_week
## get_pt_predictions_one_week
## get_dist_predictions_one_week
## get_prediction_season_week
## make_competition_forecasts_by_trajectory
## make_competition_forecasts_one_season_week
## make_competition_forecasts_one_season_week_by_trajectory
## simulate_counts_by_week
## simulate_counts_by_week_by_trajectory
## simulate_from_weighted_kde_given_ind
## simulate_from_weighted_kde


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



### functions to make predictions for the Dengue Forecasting competition

make_competition_forecasts_by_trajectory <- function(
    ssr_fits_by_prediction_horizon_limit,
    n_sims,
    data,
    outfile_path,
    location,
    phase="train") {
    
    ## make a data frame with seasons and weeks at which we make predictions
    
    ## values used by competition organizers -- I think week 0 means week 52
    ## of the previous season.
    if(identical(phase, "train")) {
        last_obs_seasons <- c("2005/2006",
            "2006/2007",
            "2007/2008",
            "2008/2009")
    } else {
        last_obs_seasons <- c("2009/2010",
            "2010/2011",
            "2011/2012",
            "2012/2013")
    }
    
    last_obs_weeks <- seq(from = 0, to = 48, by = 4)
    
    results <- expand.grid(last_obs_week=last_obs_weeks,
        last_obs_season=last_obs_seasons,
        stringsAsFactors=FALSE)
    
    ## values we use -- translate week 0 to week 52 of the previous season
    results$last_obs_season_ind1 <- results$last_obs_season
    results$last_obs_week_ind1 <- results$last_obs_week
    results$last_obs_season_ind1[results$last_obs_week == 0] <-
        sapply(results$last_obs_season_ind1[results$last_obs_week == 0],
            function(season) {
                init_season <- as.numeric(substr(season, 1, 4))
                return(paste0(init_season - 1, "/", init_season))
            })
    results$last_obs_week_ind1[results$last_obs_week == 0] <- 52
    
    ## set row names
    rownames(results) <- paste0(results$last_obs_season,
        "_wk",
        results$last_obs_week)
    
    ## create a separate results data frame for each quantity we are predicting
    ## and add extra columns representing the quantities being predicted
    peak_incidence_results <- results
    peak_incidence_results$point <- NA
    
    if(identical(location, "sanjuan")) {
        peak_incidence_dist_cutoffs <- c(50 * seq(from = 0, to = 10), Inf)
    } else {
        peak_incidence_dist_cutoffs <- c(15 * seq(from = 0, to = 10), Inf)
    }
    for(ind in seq_len(length(peak_incidence_dist_cutoffs) - 1)) {
        if(peak_incidence_dist_cutoffs[ind + 1] < Inf) {
            result_name <- paste0("p(",
                peak_incidence_dist_cutoffs[ind],
                "<=peak_incidence<",
                peak_incidence_dist_cutoffs[ind + 1],
                ")")
        } else {
            result_name <- paste0("p(",
                peak_incidence_dist_cutoffs[ind],
                "<=peak_incidence)")
        }
        peak_incidence_results[[result_name]] <- NA
    }
    
    peak_week_results <- results
    peak_week_results$point <- NA
    
    peak_week_dist_values <- seq_len(52)
    for(ind in seq_along(peak_week_dist_values)) {
        result_name <- paste0("p(peak_week=",
            peak_week_dist_values[ind],
            ")")
        peak_week_results[[result_name]] <- NA
    }
    
    season_incidence_results <- results
    season_incidence_results$point <- NA
    
    if(identical(location, "sanjuan")) {
        season_incidence_dist_cutoffs <- c(seq(from = 0, to = 10000, by = 1000), Inf)
    } else {
        season_incidence_dist_cutoffs <- c(seq(from = 0, to = 1000, by = 100), Inf)
    }
    for(ind in seq_len(length(season_incidence_dist_cutoffs) - 1)) {
        if(season_incidence_dist_cutoffs[ind + 1] < Inf) {
            result_name <- paste0("p(",
                season_incidence_dist_cutoffs[ind],
                "<=season_incidence<",
                season_incidence_dist_cutoffs[ind + 1],
                ")")
        } else {
            result_name <- paste0("p(",
                season_incidence_dist_cutoffs[ind],
                "<=season_incidence)")
        }
        season_incidence_results[[result_name]] <- NA
    }
    
    
    ## columns of each data frame where results are stored
    peak_incidence_results_prediction_cols <- seq(from=ncol(results) + 1,
        to=ncol(peak_incidence_results))
    peak_week_results_prediction_cols <- seq(from=ncol(results) + 1,
        to=ncol(peak_week_results))
    season_incidence_results_prediction_cols <- seq(from=ncol(results) + 1,
        to=ncol(season_incidence_results))
    
    
    ## get predictions for each combination of last_obs_season and last_obs_week
    for(results_row in seq_len(nrow(results))) {
        last_obs_week_ind0 <- results[results_row, "last_obs_week"]
        last_obs_week <- results[results_row, "last_obs_week_ind1"]
        last_obs_season <- results[results_row, "last_obs_season_ind1"]
        
        
        ## assemble weekly fits for the current season --
        ## either the observed counts if we are at or past that week, or
        ## the weights and kernel centers if we are not yet at that week.
        weekly_fits <- lapply(seq_len(52), function(prediction_week) {
                data_ind <- which(data$season == last_obs_season &
                        data$season_week == prediction_week)
                if(prediction_week <= last_obs_week_ind0) {
                    ## get observed value for the given week
                    return(list(observed_value=data[data_ind, "total_cases"]))
                } else {
                    prediction_horizon <- prediction_week - last_obs_week_ind0
                    prediction_horizon_limit <- 52 - last_obs_week_ind0
                    
                    ## get inds for prediction and training data -- ensure they're the
                    ## same for all weeks in the season so that we can perform prediction
                    ## by trajectory
                    ssr_fit_ind <- which(seq(from=4, to=52, by=4) == prediction_horizon_limit)
                    ssr_fit <- ssr_fits_by_prediction_horizon_limit[[ssr_fit_ind]]
                    
                    ## prediction data are those within max_lag of the last observed week
                    max_lag <- max(unlist(ssr_fit$lags_hat))
                    last_obs_week_ind <- which(data$season == last_obs_season &
                            data$season_week == last_obs_week)
                    prediction_data_inds <- seq(
                        from=last_obs_week_ind - max_lag,
                        to=last_obs_week_ind)
                    
                    if(identical(phase, "train")) {
                        ## in the training phase, forecasts use all the data not within
                        ## +/- 1 year of the last observed week as regression examples.
                        ## drop indices that are within 1 year of the last observed week or
                        ## within the last (prediction_horizon = 52 - last_obs_week_ind0)
                        ## weeks of the end of the data
                        training_data_inds_to_drop <- c(
                            seq(from=last_obs_week_ind - 52,
                                to=last_obs_week_ind + 52),
                            seq(from=nrow(data) - (52 - last_obs_week_ind0),
                                to=nrow(data)))
                        training_data_inds_to_drop <- training_data_inds_to_drop[
                            training_data_inds_to_drop >= 1 &
                                training_data_inds_to_drop <= nrow(data)
                        ]
                    } else {
                        ## in the testing phase, we need to update the ssr_fit with
                        ## additional regression examples that are available in the new
                        ## data that occur before the last observed week.
                        ## we compute the data smooths here so that data after the last
                        ## observed week are not used.
                        ## note that the fit parameter estimates are not updated.
                        training_data_inds_to_drop <- c()
                        
                        train_data <- data[seq_len(last_obs_week_ind), , drop=FALSE]
                        
                        ## add log column
                        train_data$log_total_cases <- log(train_data$total_cases + 1)
                        
                        ## add smooth log column
                        sm <- loess(log_total_cases ~ as.numeric(week_start_date),
                            data=train_data,
                            span=12 / nrow(train_data))
                        train_data$smooth_log_cases <- sm$fitted
                        
                        data$smooth_log_cases <- NA
                        data$smooth_log_cases[seq_len(last_obs_week_ind)] <- sm$fitted
                        
                        ssr_fit$train_data <- train_data
                    }
                    
                    ## get kernel weights and centers for prediction in the given
                    ## week
                    return(ssr_predict(
                            ssr_fit,
                            prediction_data=data[prediction_data_inds, , drop=FALSE],
                            leading_rows_to_drop=max(unlist(ssr_fit$lags_hat)),
                            additional_training_rows_to_drop=training_data_inds_to_drop,
                            prediction_horizon=prediction_horizon,
                            normalize_weights=TRUE))
                }
            })
        
        ## get predictions
        results_one_season_week <-
            make_competition_forecasts_one_season_week_by_trajectory(
                weekly_fits=weekly_fits,
                n_sims=n_sims,
                location=location)
        
        ## store predictions in corresponding results data frames
        peak_incidence_results[results_row,
            peak_incidence_results_prediction_cols] <- c(
                results_one_season_week$peak_incidence_pt_est,
                results_one_season_week$peak_incidence_dist_est)
        peak_week_results[results_row,
            peak_week_results_prediction_cols] <- c(
                results_one_season_week$peak_week_pt_est,
                results_one_season_week$peak_week_dist_est)
        season_incidence_results[results_row,
            season_incidence_results_prediction_cols] <- c(
                results_one_season_week$season_incidence_pt_est,
                results_one_season_week$season_incidence_dist_est)
    }
    
    ## transpose so that results are in format required for export
    results <- t(results)
    peak_incidence_results <- t(peak_incidence_results)
    peak_week_results <- t(peak_week_results)
    season_incidence_results <- t(season_incidence_results)
    
    ## save results data frames
    peak_incidence_results_for_output <- format(peak_incidence_results, digits=15)
    peak_week_results_for_output <- format(peak_week_results, digits=15)
    season_incidence_results_for_output <- format(season_incidence_results, digits=15)
    
    if(identical(phase, "train")) {
        write.csv(peak_incidence_results_for_output[-(1:4), ],
            file=file.path(outfile_path, paste0("kerneloftruth_peakinc_", location, ".csv")))
        write.csv(peak_week_results_for_output[-(1:4), ],
            file=file.path(outfile_path, paste0("kerneloftruth_peakweek_", location, ".csv")))
        write.csv(season_incidence_results_for_output[-(1:4), ],
            file=file.path(outfile_path, paste0("kerneloftruth_seasoninc_", location, ".csv")))
    } else {
        write.csv(peak_incidence_results_for_output[-(1:4), ],
            file=file.path(outfile_path, paste0("kerneloftruth_peakinc_", location, "_test.csv")))
        write.csv(peak_week_results_for_output[-(1:4), ],
            file=file.path(outfile_path, paste0("kerneloftruth_peakweek_", location, "_test.csv")))
        write.csv(season_incidence_results_for_output[-(1:4), ],
            file=file.path(outfile_path, paste0("kerneloftruth_seasoninc_", location, "_test.csv")))
    }
}

make_competition_forecasts_one_season_week <- function(weekly_fits,
    n_sims) {
    
    ## simulate n_sims realizations from the distributions implied by the
    ## weekly_kernel_fits
    ## result is a data frame with simulation index in the row,
    ## week of the season in the column
    simulated_counts_by_week <-
        simulate_counts_by_week(weekly_fits, n_sims)
    
    ## extract estimates of peak incidence
    max_counts_by_sim_ind <- apply(simulated_counts_by_week, 1, max)
    
    peak_incidence_pt_est <- mean(max_counts_by_sim_ind)
    
    peak_incidence_dist_cutoffs <- c(15 * seq(from = 0, to = 10), Inf)
    peak_incidence_dist_est <- sapply(
        seq_len(length(peak_incidence_dist_cutoffs) - 1),
        function(lb_ind) {
            mean(max_counts_by_sim_ind >= peak_incidence_dist_cutoffs[lb_ind] &
                    max_counts_by_sim_ind < peak_incidence_dist_cutoffs[lb_ind + 1])
        })
    
    ## extract estimates of peak week
    peak_week_by_sim_ind <- apply(simulated_counts_by_week, 1, which.max)
    
    peak_week_pt_est <- mean(peak_week_by_sim_ind)
    
    peak_week_dist_est <- sapply(seq_len(52), function(week) {
            mean(peak_week_by_sim_ind == week)
        })
    
    ## extract estimates of season incidence
    season_incidence_by_sim_ind <- apply(simulated_counts_by_week, 1, sum)
    
    season_incidence_pt_est <- mean(season_incidence_by_sim_ind)
    
    season_incidence_dist_cutoffs <- c(seq(from = 0, to = 1000, by = 100), Inf)
    season_incidence_dist_est <- sapply(
        seq_len(length(season_incidence_dist_cutoffs) - 1),
        function(lb_ind) {
            mean(season_incidence_by_sim_ind >= season_incidence_dist_cutoffs[lb_ind] &
                    season_incidence_by_sim_ind < season_incidence_dist_cutoffs[lb_ind + 1])
        })
    
    ## In theory, our predictive distributions assign non-zero probability to
    ## each bin.  However, our simulation-based procedure for estimating these
    ## probabilities may result in estimated probabilities of 0.  In order to
    ## ensure that we achieve a finite score, we inflate the estimated probability
    ## of each class where the simulation-based estimate is 0.
    peak_incidence_dist_est <- peak_incidence_dist_est +
        (0.0001 / length(peak_incidence_dist_est))
    peak_incidence_dist_est <- peak_incidence_dist_est / sum(peak_incidence_dist_est)
    
    peak_week_dist_est <- peak_week_dist_est +
        (0.0001 / length(peak_week_dist_est))
    peak_week_dist_est <- peak_week_dist_est / sum(peak_week_dist_est)
    
    peak_incidence_dist_est <- peak_incidence_dist_est +
        (0.0001 / length(peak_incidence_dist_est))
    peak_incidence_dist_est <- peak_incidence_dist_est / sum(peak_incidence_dist_est)
    
    return(list(
            peak_incidence_pt_est=peak_incidence_pt_est,
            peak_incidence_dist_est=peak_incidence_dist_est,
            peak_week_pt_est=peak_week_pt_est,
            peak_week_dist_est=peak_week_dist_est,
            season_incidence_pt_est=season_incidence_pt_est,
            season_incidence_dist_est=season_incidence_dist_est
        ))
}



make_competition_forecasts_one_season_week_by_trajectory <-
    function(weekly_fits,
        n_sims,
        location) {
    
    ## simulate n_sims realizations from the distributions implied by the
    ## weekly_kernel_fits
    ## result is a data frame with simulation index in the row,
    ## week of the season in the column
    simulated_counts_by_week <-
        simulate_counts_by_week_by_trajectory(weekly_fits, n_sims)
    
    ## extract estimates of peak incidence
    max_counts_by_sim_ind <- apply(simulated_counts_by_week, 1, max)
    
    peak_incidence_pt_est <- mean(max_counts_by_sim_ind)
    
    if(identical(location, "sanjuan")) {
        peak_incidence_dist_cutoffs <- c(50 * seq(from = 0, to = 10), Inf)
    } else {
        peak_incidence_dist_cutoffs <- c(15 * seq(from = 0, to = 10), Inf)
    }
    peak_incidence_dist_est <- sapply(
        seq_len(length(peak_incidence_dist_cutoffs) - 1),
        function(lb_ind) {
            mean(max_counts_by_sim_ind >= peak_incidence_dist_cutoffs[lb_ind] &
                    max_counts_by_sim_ind < peak_incidence_dist_cutoffs[lb_ind + 1])
        })
    
    ## extract estimates of peak week
    peak_week_by_sim_ind <- apply(simulated_counts_by_week, 1, which.max)
    
    peak_week_pt_est <- mean(peak_week_by_sim_ind)
    
    peak_week_dist_est <- sapply(seq_len(52), function(week) {
            mean(peak_week_by_sim_ind == week)
        })
    
    ## extract estimates of season incidence
    season_incidence_by_sim_ind <- apply(simulated_counts_by_week, 1, sum)
    
    season_incidence_pt_est <- mean(season_incidence_by_sim_ind)
    
    if(identical(location, "sanjuan")) {
        season_incidence_dist_cutoffs <- c(seq(from = 0, to = 10000, by = 1000), Inf)
    } else {
        season_incidence_dist_cutoffs <- c(seq(from = 0, to = 1000, by = 100), Inf)
    }
    season_incidence_dist_est <- sapply(
        seq_len(length(season_incidence_dist_cutoffs) - 1),
        function(lb_ind) {
            mean(season_incidence_by_sim_ind >= season_incidence_dist_cutoffs[lb_ind] &
                    season_incidence_by_sim_ind < season_incidence_dist_cutoffs[lb_ind + 1])
        })
    
    ## In theory, our predictive distributions assign non-zero probability to
    ## each bin.  However, our simulation-based procedure for estimating these
    ## probabilities may result in estimated probabilities of 0.  In order to
    ## ensure that we achieve a finite score, we inflate the estimated probability
    ## of each class where the simulation-based estimate is 0.
    peak_incidence_dist_est <- peak_incidence_dist_est +
        (0.0001 / length(peak_incidence_dist_est))
    peak_incidence_dist_est <- peak_incidence_dist_est / sum(peak_incidence_dist_est)
    
    peak_week_dist_est <- peak_week_dist_est +
        (0.0001 / length(peak_week_dist_est))
    peak_week_dist_est <- peak_week_dist_est / sum(peak_week_dist_est)
    
    season_incidence_dist_est <- season_incidence_dist_est +
        (0.0001 / length(season_incidence_dist_est))
    season_incidence_dist_est <- season_incidence_dist_est / sum(season_incidence_dist_est)
    
    return(list(
            peak_incidence_pt_est=peak_incidence_pt_est,
            peak_incidence_dist_est=peak_incidence_dist_est,
            peak_week_pt_est=peak_week_pt_est,
            peak_week_dist_est=peak_week_dist_est,
            season_incidence_pt_est=season_incidence_pt_est,
            season_incidence_dist_est=season_incidence_dist_est
        ))
}

simulate_counts_by_week <- function(weekly_fits, n_sims) {
    results <- sapply(seq_len(52), function(week) {
            if(length(weekly_fits[[week]]$observed_value) == 1) {
                return(rep(weekly_fits[[week]]$observed_value, n_sims))
            } else {
                return(simulate_from_weighted_kde(n_sims, weekly_fits[[week]]))
            }
        })
    
    return(results)
}

simulate_counts_by_week_by_trajectory <- function(weekly_fits, n_sims) {
    ## get weights -- same for all weekly fits.
    weights <- weekly_fits[[52]]$weights
    
    ## select index of the trajectory to sample from for each simulation
    trajectory_index <- sample(length(weights),
        size=n_sims,
        replace=TRUE,
        prob=weights)
    
    ## return observed value if available or simulated value
    results <- sapply(seq_len(52), function(week) {
            if(length(weekly_fits[[week]]$observed_value) == 1) {
                return(rep(weekly_fits[[week]]$observed_value, n_sims))
            } else {
                return(simulate_from_weighted_kde_given_ind(trajectory_index,
                        weekly_fits[[week]]))
            }
        })
    
    return(results)
}

simulate_from_weighted_kde_given_ind <- function(inds, weighted_kde_fit) {
    ## get bandwidth
    density_fit <- density(x=weighted_kde_fit$centers[, 1],
        weights=weighted_kde_fit$weights,
        bw="SJ")
    
    ## get simulated values -- from a mixture of normals with sd=density_fit$bw
    component_means <- weighted_kde_fit$centers[inds, 1]
    
    return(rnorm(length(inds), component_means, density_fit$bw))
}

simulate_from_weighted_kde <- function(n, weighted_kde_fit) {
    ## get bandwidth
    density_fit <- density(x=weighted_kde_fit$centers[, 1],
        weights=weighted_kde_fit$weights,
        bw="SJ")
    
    ## get simulated values -- from a mixture of normals with sd=density_fit$bw
    component_means <- sample(weighted_kde_fit$centers[, 1],
        size=n,
        replace=TRUE,
        prob=weighted_kde_fit$weights)
    
    return(rnorm(n, component_means, density_fit$bw))
}


