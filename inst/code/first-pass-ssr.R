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
target_inds <- which(sj$season == "2008/2009")
ssr_ests <- rbind.fill(lapply(c(1, 13, 26, 39, 52), function(prediction_steps) {
    last_obs_season_week_and_prediction_steps_combos <- lapply(target_inds,
        function(ti) {
            list(last_obs_season = as.character(sj$season[ti - prediction_steps]),
                last_obs_week = sj$season_week[ti - prediction_steps],
                prediction_steps = prediction_steps)
        }
    )
        
    rbind.fill(lapply(last_obs_season_week_and_prediction_steps_combos,
        function(params) {
            ssr_predict_dengue_stepsahead_one_week(
                last_obs_season = params$last_obs_season,
                last_obs_week = params$last_obs_week,
                theta = 10,
                tr_lag = 1,
                prediction_steps = params$prediction_steps,
                data = sj,
                prediction_types = "density")$density_preds
        }
    ))
}))


sj$season_season_week <- paste(as.character(sj$season), sj$season_week, sep = "-")
ssr_ests$season_season_week <- paste(ssr_ests$season, ssr_ests$season_week,
    sep = "-")

ssr_ests$model <- paste0("ssr_p_", ssr_ests$prediction_step)

## get naive model estimates -- only handles holdout_seasons with length 1?
get_season_week_kde_ests_one_week <- function(wk, holdout_seasons, data) {
    inds <- which(!(data$season %in% holdout_seasons) & data$season_week == wk)
    
    ## do weighted kde
    kde_log_scale <- density(data[inds, "log_total_cases"],
        bw = "SJ")
    kde_orig_scale <- density(data[inds, "total_cases"],
        bw = "SJ")
    
    return(data.frame(log_total_cases = kde_log_scale$x,
            total_cases = kde_orig_scale$x,
            est_density_log_scale = kde_log_scale$y,
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
p <- ggplot() +
    geom_line(aes(x = log_total_cases, y = est_density_log_scale, colour = model), data = combined_ests) +
    geom_vline(aes(xintercept = value, linetype = variable), data = obs_counts[obs_counts$variable %in% c("smooth_log_cases", "log_total_cases"), ], show_guide = TRUE) +
    facet_wrap(~ season_week) +
    ggtitle("Predictions by week for 2008/2009 season (log scale)") +
    theme_bw()

print(p)

## original scale
p <- ggplot() +
    geom_line(aes(x = total_cases, y = est_density_orig_scale, colour = model), data = combined_ests) +
    geom_vline(aes(xintercept = value, linetype = variable), data = obs_counts[obs_counts$variable %in% c("exp_smooth_log_cases", "total_cases"), ], show_guide = TRUE) +
    coord_cartesian(xlim = c(0, 50)) +
    facet_wrap(~ season_week) +
    ggtitle("Predictions by week for 2008/2009 season") +
    theme_bw()

print(p)




## compute point estimates and get MASE for each lag and naive model
target_inds <- which(sj$season == "2008/2009")
prediction_steps_vec <- 1:52
#prediction_steps_vec <- c(1, 13, 26, 39, 52)
ssr_pt_ests <- rbind.fill(lapply(prediction_steps_vec, function(prediction_steps) {
    last_obs_season_week_and_prediction_steps_combos <- lapply(target_inds,
        function(ti) {
            list(last_obs_season = as.character(sj$season[ti - prediction_steps]),
                last_obs_week = sj$season_week[ti - prediction_steps],
                prediction_steps = prediction_steps)
         }
    )

    rbind.fill(lapply(last_obs_season_week_and_prediction_steps_combos,
        function(params) {
            ssr_predict_dengue_stepsahead_one_week(
                last_obs_season = params$last_obs_season,
                last_obs_week = params$last_obs_week,
                theta = 10,
                tr_lag = 1,
                prediction_steps = params$prediction_steps,
                data = sj,
                prediction_types = "pt")$pt_preds
        }
    ))
}))


sj$season_season_week <- paste(as.character(sj$season), sj$season_week, sep = "-")
ssr_pt_ests$season_season_week <- paste(ssr_pt_ests$season, ssr_pt_ests$season_week,
    sep = "-")

ssr_pt_ests$model <- paste0("ssr_p_", ssr_pt_ests$prediction_step)

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
        holdout_seasons = "2008/2009",
        data = sj)
)

season_week_pt_ests$model <- "naive"



# combine estimates from ssr with various lags and naive model
combined_pt_ests <- rbind.fill(ssr_pt_ests, season_week_pt_ests)


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


