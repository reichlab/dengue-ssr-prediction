## read in and plot data

library(lubridate)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape)
library(ssr)

sj <- San_Juan_train

## add lag columns
sj <- sj %>% mutate(total_cases_lag3 = lag(total_cases, 3),
                    total_cases_lag1 = lag(total_cases))

## add log column
sj$log_total_cases <- log(sj$total_cases + 1)

## add smooth log column
sm <- loess(log_total_cases ~ as.numeric(week_start_date), data=sj, span=1/80)
sj$smooth_log_cases <- sm$fitted

                
        
## simple time series plot
ggplot(sj, aes(x=week_start_date, y=total_cases)) + geom_line()
ggplot(sj, aes(x=week_start_date)) + geom_linerange(aes(ymin=0, ymax=total_cases))

## phase plots
(p <- ggplot(sj, aes(x=total_cases_lag1, y=total_cases, color=season)) + 
    geom_point() + geom_path() +
    scale_x_log10() + scale_y_log10() +
    geom_abline(b=1, a=0, color="grey"))

p + facet_wrap(~ season)

## smooth data
sm <- loess(total_cases ~ as.numeric(week_start_date), data=sj, span=1/80)
sj$smooth_cases <- sm$fitted
sj$smooth_cases_lag1 <- lag(sj$smooth_cases)

ggplot(sj, aes(x=week_start_date)) + 
    geom_linerange(aes(ymin=0, ymax=total_cases)) +
    geom_line(aes(y=smooth_cases), color="red", size=1)

## phase plot with smooth data
(p1 <- ggplot(sj, aes(x=smooth_cases_lag1, y=smooth_cases, color=season)) + 
    geom_point() + geom_path() +
    scale_x_log10() + scale_y_log10() +
    geom_abline(b=1, a=0, color="grey"))

p1 + facet_wrap(~ season)







# 1, 13, 26, 39, and 52 weeks ahead density predictions, for plots
options(error = recover)
target_inds <- which(sj$season == "2008/2009")
#prediction_horizon <- 1

ssr_ests <- rbind.fill(lapply(c(1, 13, 26, 39, 52), function(prediction_horizon) {
    last_obs_season_week_and_prediction_horizon_combos <- lapply(target_inds,
        function(ti) {
            list(last_obs_season = as.character(sj$season[ti - prediction_horizon]),
                last_obs_week = sj$season_week[ti - prediction_horizon],
                prediction_horizon = prediction_horizon)
        }
    )
#    params <- last_obs_season_week_and_prediction_horizon_combos[[1]]
    
    rbind.fill(lapply(last_obs_season_week_and_prediction_horizon_combos,
        function(params) {
            preds <- ssr_predict_dengue_one_week(
                last_obs_season = params$last_obs_season,
                last_obs_week = params$last_obs_week,
                lags = list(smooth_log_cases = c(0, 1)),
                theta=list(smooth_log_cases_lag0=list(bw=.1),
                    smooth_log_cases_lag1=list(bw=.1)),
                prediction_horizon = params$prediction_horizon,
                data = sj,
				season_var = "season",
				week_var = "season_week",
                X_names = c("smooth_log_cases"),
                y_name = c("total_cases"),
                time_name = "week_start_date",
                prediction_types = "density")$dist_preds
            preds$last_obs_season <- params$last_obs_season
            preds$last_obs_week <- params$last_obs_week
            last_obs_ind <- which(sj[["season"]] == params$last_obs_season &
                    sj[["season_week"]] == params$last_obs_week)
            temp <- get_prediction_season_week(last_obs_ind,
                params$prediction_horizon,
                sj,
                season_var = "season",
                week_var = "season_week")
            preds$prediction_season <- temp$season
            preds$prediction_week <- temp$week
            preds$prediction_horizon <- params$prediction_horizon
            
            return(preds)
        }
    ))
}))


sj$season_season_week <- paste(as.character(sj$season), sj$season_week, sep = "-")
ssr_ests$season_season_week <- paste(ssr_ests$prediction_season, ssr_ests$prediction_week,
    sep = "-")

ssr_ests$model <- paste0("ssr_p_", ssr_ests$prediction_horizon)

## get naive model estimates -- only handles holdout_seasons with length 1?
get_season_week_kde_ests_one_week <- function(wk, holdout_seasons, data) {
    inds <- which(!(data$season %in% holdout_seasons) & data$season_week == wk)
    
    ## do weighted kde
#    kde_log_scale <- density(data[inds, "log_total_cases"],
#        bw = "SJ")
    kde_orig_scale <- density(data[inds, "total_cases"],
        bw = "SJ")
    
    return(data.frame(#log_total_cases = kde_log_scale$x,
            total_cases = kde_orig_scale$x,
#            est_density_log_scale = kde_log_scale$y,
            est_density_orig_scale = kde_orig_scale$y,
            season = rep(holdout_seasons, each = length(kde_orig_scale$x)),
            season_week = wk
        ))
}

season_week_kde_ests <- rbind.fill(
    lapply(seq_len(52), get_season_week_kde_ests_one_week,
        holdout_seasons = "2008/2009",
        data = sj)
)

season_week_kde_ests$model <- "naive"



# combine estimates from ssr with various lags and naive model
ssr_ests$total_cases <- ssr_ests$x
ssr_ests$est_density_orig_scale <- ssr_ests$est_density
ssr_ests$season_week <- ssr_ests$prediction_week
combined_ests <- rbind.fill(ssr_ests, season_week_kde_ests)


# get observed totals
sj_inds_keep <- seq_len(nrow(sj))[sj$season_season_week %in%
    ssr_ests$season_season_week]
obs_counts <- sj[sj_inds_keep,
    c("smooth_log_cases", "log_total_cases", "total_cases", "season_season_week", "season", "season_week")]
obs_counts$exp_smooth_log_cases <- exp(obs_counts$smooth_log_cases)
obs_counts <- melt(obs_counts, id = c("season", "season_week", "season_season_week"))


## plot estimates, compared with observed cases and smoothed cases
## log scale
# p <- ggplot() +
#     geom_line(aes(x = log_total_cases, y = est_density_log_scale, colour = model), data = combined_ests) +
#     geom_vline(aes(xintercept = value, linetype = variable), data = obs_counts[obs_counts$variable %in% c("smooth_log_cases", "log_total_cases"), ], show_guide = TRUE) +
#     facet_wrap(~ season_week) +
#     ggtitle("Predictions by week for 2008/2009 season (log scale)") +
#     theme_bw()
# 
# print(p)

## original scale
p <- ggplot() +
    geom_line(aes(x = total_cases, y = est_density_orig_scale, colour = model), data = combined_ests) +
    geom_line(aes(x = total_cases, y = est_density, colour = model), data = combined_ests) +
    geom_vline(aes(xintercept = value, linetype = variable), data = obs_counts[obs_counts$variable %in% c("exp_smooth_log_cases", "total_cases"), ], show_guide = TRUE) +
    coord_cartesian(xlim = c(0, 50)) +
    facet_wrap(~ season_week) +
    ggtitle("Predictions by week for 2008/2009 season") +
    theme_bw()

print(p)




## compute point estimates and get MASE for each lag and naive model
holdout_seasons <- c("2007/2008", "2008/2009")
target_inds <- which(sj$season %in% holdout_seasons)
prediction_horizons_vec <- 1:52
ssr_pt_ests <- rbind.fill(lapply(prediction_horizons_vec, function(prediction_horizon) {
    last_obs_season_week_and_prediction_horizon_combos <- lapply(target_inds,
        function(ti) {
            list(last_obs_season = as.character(sj$season[ti - prediction_horizon]),
                last_obs_week = sj$season_week[ti - prediction_horizon],
                prediction_horizon = prediction_horizon)
        }
    )

    rbind.fill(lapply(last_obs_season_week_and_prediction_horizon_combos,
        function(params) {
            preds <- ssr_predict_dengue_one_week(
                last_obs_season = params$last_obs_season,
                last_obs_week = params$last_obs_week,
                lags = list(smooth_log_cases = c(0, 1)),
                theta=list(smooth_log_cases_lag0=list(bw=.1),
                    smooth_log_cases_lag1=list(bw=.1)),
                prediction_horizon = params$prediction_horizon,
                data = sj,
                season_var = "season",
                week_var = "season_week",
                X_names = c("smooth_log_cases"),
                y_name = c("total_cases"),
                time_name = "week_start_date",
                prediction_types = "pt")$pt_preds
            preds$last_obs_season <- params$last_obs_season
            preds$last_obs_week <- params$last_obs_week
            last_obs_ind <- which(sj[["season"]] == params$last_obs_season &
                    sj[["season_week"]] == params$last_obs_week)
            temp <- get_prediction_season_week(last_obs_ind,
                params$prediction_horizon,
                sj,
                season_var = "season",
                week_var = "season_week")
            preds$prediction_season <- temp$season
            preds$prediction_week <- temp$week
            preds$prediction_horizon <- params$prediction_horizon
            
            return(preds)
        }
    ))
}))


sj$season_season_week <- paste(as.character(sj$season), sj$season_week, sep = "-")
ssr_pt_ests$est_total_cases_orig_scale <- ssr_pt_ests$pt_est
ssr_pt_ests$season_week <- ssr_pt_ests$prediction_week
ssr_pt_ests$season <- ssr_pt_ests$prediction_season
ssr_pt_ests$season_season_week <- paste(ssr_pt_ests$season, ssr_pt_ests$season_week,
    sep = "-")

ssr_pt_ests$model <- paste0("ssr_p_", ssr_pt_ests$prediction_horizon)

## get naive model point estimates -- only handles holdout_seasons with length 1?
## turn this and naive density estimates function above into more real functions
get_season_week_pt_ests_one_week <- function(wk, holdout_seasons, data) {
    inds <- which(!(data$season %in% holdout_seasons) & data$season_week == wk)
    
    ## do (unweighted) mean
    pt_est_log_scale <- mean(data[inds, "log_total_cases"])
    pt_est_orig_scale <- mean(data[inds, "total_cases"])
    
    return(data.frame(est_total_cases_log_scale = pt_est_log_scale,
            est_total_cases_orig_scale = pt_est_orig_scale,
            est_total_cases_orig_scale_from_log_scale = exp(pt_est_log_scale),
            season = holdout_seasons,
            season_week = wk
        ))
}

season_week_pt_ests <- rbind.fill(
    lapply(seq_len(52), get_season_week_pt_ests_one_week,
        holdout_seasons = holdout_seasons,
        data = sj)
)

season_week_pt_ests$model <- "naive"



# combine estimates from ssr with various lags and naive model
combined_pt_ests <- rbind.fill(ssr_pt_ests, season_week_pt_ests)


sj_inds_keep <- seq_len(nrow(sj))[sj$season_season_week %in%
				ssr_pt_ests$season_season_week]
obs_counts <- sj[sj_inds_keep,
		c("smooth_log_cases", "log_total_cases", "total_cases", "season_season_week", "season", "season_week")]
obs_counts$exp_smooth_log_cases <- exp(obs_counts$smooth_log_cases)
obs_counts <- melt(obs_counts, id = c("season", "season_week", "season_season_week"))


# get mase for all models and plot
mase_by_model <- tapply(combined_pt_ests$est_total_cases_orig_scale,
    combined_pt_ests$model,
    function(preds) {
        mase(obs_counts$value[obs_counts$variable == "total_cases"], preds)
    })

mase_by_pred_step <- data.frame(
    pred_step = as.integer(sapply(strsplit(names(mase_by_model[2:53]), "_"), function(lc) lc[length(lc)])),
    mase = mase_by_model[2:53]
)

ggplot() +
    geom_point(aes(x = pred_step, y = mase), data = mase_by_pred_step) +
    geom_line(aes(x = pred_step, y = mase), data = mase_by_pred_step) +
    geom_hline(yintercept = mase_by_model[names(mase_by_model) == "naive"]) +
    coord_cartesian(ylim = c(0, 5.1)) +
    ggtitle("MASE for SSR at prediction steps = 1, ..., 52\ncompared with naive weekly model\n in 2008/2009 season") +
    theme_bw()


















### get ssr fit
library(doMC)

registerDoMC(cores=5)        ## number of cores on this machine

set.seed(998468235L,kind="L'Ecuyer")
mcopts <- list(preschedule=FALSE,set.seed=TRUE)


ssr_control <- create_ssr_control(X_names="smooth_log_cases",
    y_names="total_cases",
    time_name=NULL,
    max_lag=list(smooth_log_cases=1),
    prediction_horizon=1,
    kernel_fns=list(smooth_log_cases="squared_exp_kernel"),
    theta_est=list(smooth_log_cases="bw"),
    theta_fixed=list(),
    theta_transform_fns=list(
        squared_exp_kernel=list(
            bw=list(transform="log",
                detransform="exp")
        )
    ),
    crossval_buffer=52,
    loss_fn_name="mae_from_kernel_weights_and_centers",
    loss_fn_args=list())

options(error=recover)
debug(ssr_crossval_estimate_parameter_loss)
ssr_fit <- ssr(X_names="smooth_log_cases",
    y_names="total_cases",
    time_name=NULL,
    data=sj,
    ssr_control=ssr_control)
