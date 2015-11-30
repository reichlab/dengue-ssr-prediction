library(ggplot2)
library(ssr)
library(plyr)
library(dplyr)
library(tidyr)

load("/media/evan/data/Reich/dengue-ssr-prediction/inst/intermediate-results/ssr-poster/predictions/ssr-predictions.Rdata")
load("/media/evan/data/Reich/dengue-ssr-prediction/inst/intermediate-results/ssr-poster/predictions/sarima-predictions.Rdata")

ssr_predictions_df$model <- "ssr"

predictions_df <- rbind.fill(ssr_predictions_df, sarima_predictions_df)
predictions_df$total_cases <- predictions_df$prediction

realized_values_df <- data.frame(
	ph = rep(seq_len(52), each = length(unique(San_Juan_test$week_start_date))),
	week_start_date = rep(unique(San_Juan_test$week_start_date), times = 52L)
)
sj_lookup_inds <- sapply(realized_values_df$week_start_date, function(rvwsd) { which(San_Juan_test$week_start_date == rvwsd) })
realized_values_df$total_cases <- San_Juan_test$total_cases[sj_lookup_inds]
realized_values_df$prediction_season <- San_Juan_test$season[sj_lookup_inds]
realized_values_df$prediction_week <- San_Juan_test$season_week[sj_lookup_inds]

load("/media/evan/data/Reich/dengue-ssr-prediction/inst/intermediate-results/ssr-poster/predictions/ssr-predictions-orig-scale.Rdata")

predictions_df$total_cases[predictions_df$model == "ssr"] <- 
	predictions_df$prediction[predictions_df$model == "ssr"] <-
	ssr_predictions_df$prediction
	
ribbons_df <- predictions_df %>%
	select(week_start_date, prediction_season, prediction_week, last_obs_season, last_obs_week, ph, model, predictive_50pct_lb:predictive_90pct_ub) %>%
	gather("bound_type", "predictive_value", predictive_50pct_lb:predictive_90pct_ub) %>%
	mutate(interval_type = ifelse(grepl("50", bound_type), "50", "90"),
		bound_type = ifelse(grepl("lb", bound_type), "lower", "upper")) %>%
	spread(bound_type, predictive_value)

plot_df <- rbind.fill(predictions_df, realized_values_df)

ggplot() +
	geom_line(aes(x = week_start_date, y = total_cases), data = San_Juan_test[San_Juan_test$season %in% c("2007/2008", "2008/2009", "2009/2010", "2010/2011", "2011/2012", "2012/2013"), ]) +
	geom_line(aes(x = week_start_date, y = prediction, colour = model),
		data = predictions_df[predictions_df$ph %in% c(1, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52),]) +
	facet_wrap( ~ ph) +
	theme_bw()


#phs_used <- c(1, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52)
#all_years_used <- c("2007/2008", "2008/2009", "2009/2010", "2010/2011", "2011/2012", "2012/2013")
all_years_used <- c("2008/2009", "2009/2010", "2010/2011", "2011/2012", "2012/2013")
phs_used <- c(1, 4, 13, 26, 39, 52)
ggplot() +
	geom_ribbon(aes(x = week_start_date, ymin = lower, ymax = upper, colour = model, fill = model, alpha = interval_type),
		size = 0,
		data = ribbons_df[ribbons_df$ph %in% phs_used, ]) +
	geom_line(aes(x = week_start_date, y = total_cases), data = San_Juan_test[San_Juan_test$season %in% all_years_used, ]) +
	geom_line(aes(x = week_start_date, y = prediction, colour = model),
		data = predictions_df[predictions_df$ph %in% phs_used, ]) +
	scale_colour_manual("Model",
		breaks = c("sarima", "ssr"),
		labels = c("SARIMA", "SSR"),
		values = c("#E69F00", "#56B4E9")
		) +
	scale_fill_manual("Model",
		breaks = c("sarima", "ssr"),
		labels = c("SARIMA", "SSR"),
		values = c("#E69F00", "#56B4E9")
		) +
	scale_alpha_discrete("Prediction\nInterval\nCoverage",
		labels = c("50 Percent", "90 Percent"),
		limits = c("50", "90"),
		range = c(0.4, 0.2)) +
	facet_wrap( ~ ph, ncol = 2) +
	xlab("Time") +
	ylab("Total Cases") +
	ggtitle("Point and Interval Predictions\nFaceted by Prediction Horizon") +
	theme_bw(base_size = 32)





all_years_used <- c("2009/2010", "2010/2011", "2011/2012", "2012/2013")
all_last_obs_weeks_used <- c(1, 13, 26)
San_Juan_test$prediction_season <- San_Juan_test$season
San_Juan_test$prediction_week <- San_Juan_test$season_week

ribbons_df$prediction_season_for_last_obs_ribbon_plot <- ribbons_df$prediction_season
ribbons_df$prediction_season_for_last_obs_ribbon_plot[
	ribbons_df$last_obs_season == ribbons_df$prediction_season &
	ribbons_df$prediction_week < ribbons_df$last_obs_week] <- NA

predictions_df$prediction_season_for_last_obs_ribbon_plot <- ribbons_df$prediction_season
predictions_df$prediction_season_for_last_obs_ribbon_plot[
	predictions_df$last_obs_season == predictions_df$prediction_season &
	predictions_df$prediction_week < predictions_df$last_obs_week] <- NA

facet_labels <- function(variable, value) {
	if(identical(variable, "last_obs_week")) {
		return(paste0("Last Observed Week: ", value))
	} else {
		return(as.character(value))
	}
}

p <- ggplot() +
#	geom_ribbon(aes(x = prediction_week, ymin = lower, ymax = upper, colour = model, fill = model, alpha = interval_type),
	geom_ribbon(aes(x = prediction_week, ymin = lower, ymax = upper, fill = model),
		alpha = 0.3,
		size = 0,
		data = ribbons_df[ribbons_df$interval_type == "50" &
			ribbons_df$last_obs_week %in% all_last_obs_weeks_used &
			(ribbons_df$last_obs_season == ribbons_df$prediction_season &
				ribbons_df$prediction_week > ribbons_df$last_obs_week), ]) +
	geom_line(aes(x = prediction_week, y = total_cases), data = San_Juan_test[San_Juan_test$season %in% all_years_used, ]) +
	geom_line(aes(x = prediction_week, y = prediction, colour = model, linetype = model),
		data = predictions_df[predictions_df$last_obs_week %in% all_last_obs_weeks_used &
			(predictions_df$last_obs_season == predictions_df$prediction_season &
				predictions_df$prediction_week > predictions_df$last_obs_week), ]) +
	scale_linetype_manual("Model",
		breaks = c("sarima", "ssr"),
		labels = c("SARIMA", "SSR"),
		values = c(2, 6)
	) +
	scale_colour_manual("Model",
		breaks = c("sarima", "ssr"),
		labels = c("SARIMA", "SSR"),
		values = c("#E69F00", "#56B4E9")
		) +
	scale_fill_manual("Model",
		breaks = c("sarima", "ssr"),
		labels = c("SARIMA", "SSR"),
		values = c("#E69F00", "#56B4E9")
		) +
#	scale_alpha_discrete("Prediction\nInterval\nCoverage",
#		labels = c("50 Percent", "90 Percent"),
#		limits = c("50", "90"),
#		range = c(0.4, 0.2)) +
	facet_grid(prediction_season ~ last_obs_week,
		labeller = "facet_labels") +
	xlab("Week of Season") +
	ylab("Total Cases") +
	ggtitle("Point and 50% Interval Predictions for Each Test Season\nFaceted by Season and Last Observed Week") +
	theme_bw(base_size = 22)

pdf("/media/evan/data/Reich/dengue-ssr-prediction/inst/intermediate-results/ssr-poster/plots/Predictions-by-season-and-last-obs-week.pdf", width = 20, height = 15)
print(p)
dev.off()
















ggplot() +
	geom_line(aes(x = ph, y = AE, colour = model), stat="summary", fun.y = "mean", data = predictions_df) +
	theme_bw()

	
ggplot() +
	geom_point(aes(x = ph, y = AE, colour = model), data = predictions_df) +
	theme_bw()

	
pointwise_result_differences <- data.frame(
	ph = rep(seq_len(52), each = length(unique(predictions_df$week_start_date))),
	week_start_date = rep(unique(predictions_df$week_start_date), times = 52L)
)
sj_lookup_inds <- sapply(pointwise_result_differences$week_start_date, function(rvwsd) { which(San_Juan_test$week_start_date == rvwsd) })
pointwise_result_differences$total_cases <- San_Juan_test$total_cases[sj_lookup_inds]
pointwise_result_differences$prediction_season <- San_Juan_test$season[sj_lookup_inds]
pointwise_result_differences$prediction_week <- San_Juan_test$season_week[sj_lookup_inds]

sarima_lookup_inds <- sapply(seq_len(nrow(pointwise_result_differences)),
		function(row_ind) {
			which(predictions_df$ph == pointwise_result_differences$ph[row_ind] & predictions_df$week_start_date == pointwise_result_differences$week_start_date[row_ind] & predictions_df$model == "sarima")
		})
ssr_lookup_inds <- sapply(seq_len(nrow(pointwise_result_differences)),
		function(row_ind) {
			which(predictions_df$ph == pointwise_result_differences$ph[row_ind] & predictions_df$week_start_date == pointwise_result_differences$week_start_date[row_ind] & predictions_df$model == "ssr")
		})
	
pointwise_result_differences$AE_diff <-
	predictions_df$total_cases[sarima_lookup_inds] -
	predictions_df$total_cases[ssr_lookup_inds]

ggplot() +
	geom_violin(aes(x = ph, y = AE_diff), data = pointwise_result_differences) +
	theme_bw()


pointwise_result_differences$log_score_diff <-
	predictions_df$log_score[ssr_lookup_inds] -
	predictions_df$log_score[sarima_lookup_inds]


pointwise_result_differences$ph <- as.factor(pointwise_result_differences$ph)

dist_summaries_df <- data.frame(
	ph = factor(1:52),
	mean = sapply(1:52, function(ph) {
		mean(pointwise_result_differences$log_score_diff[pointwise_result_differences$ph == ph])
	}),
	median = sapply(1:52, function(ph) {
		median(pointwise_result_differences$log_score_diff[pointwise_result_differences$ph == ph])
	})
)

dist_summaries_df <- melt(dist_summaries_df, "ph")

p <- ggplot() +
	geom_violin(aes(x = ph, y = log_score_diff), data = pointwise_result_differences) +
	geom_hline(yintercept=0) +
	geom_point(aes(x = ph, y = value, colour = variable), shape = "-", size = 8, data = dist_summaries_df) +
	scale_colour_manual("Summary",
		values = c("#E69F00", "#56B4E9")) +
	ylab("Difference in Log Scores\n(Positive Values Indicate SSR outperformed SARIMA)") +
	xlab("Prediction Horizon") +
	ggtitle("Difference in Log Score vs. Prediction Horizon") +
	theme_bw(base_size = 22)


pdf("/media/evan/data/Reich/dengue-ssr-prediction/inst/intermediate-results/ssr-poster/plots/log-score-diff-by-ph.pdf", width=20, height=9.5)
print(p)
dev.off()

#ggplot() +
#	geom_violin(aes(x = ph, y = log_score_diff), data = pointwise_result_differences) +
#	geom_point(aes(x = ph, y = log_score_diff), stat="summary", fun.y="mean", data = pointwise_result_differences) +
#	geom_point(aes(x = ph, y = log_score_diff), stat="summary", colour="red", fun.y="median", data = pointwise_result_differences) +
#	facet_wrap(~ prediction_season, ncol = 1) +
#	theme_bw()


#ggplot() +
#	geom_violin(aes(x = prediction_season, y = log_score_diff), data = pointwise_result_differences) +
#	geom_point(aes(x = prediction_season, y = log_score_diff), stat="summary", fun.y="mean", data = pointwise_result_differences) +
#	geom_point(aes(x = prediction_season, y = log_score_diff), stat="summary", colour="red", fun.y="median", data = pointwise_result_differences) +
#	theme_bw()



ggplot() +
#	geom_point(aes(x = total_cases, y = log_score_diff, colour = prediction_season, shape = prediction_season), stat="summary", fun.y="mean", data = pointwise_result_differences) +
	geom_point(aes(x = total_cases, y = log_score_diff, colour = prediction_season, shape = prediction_season, alpha = as.numeric(ph)), data = pointwise_result_differences) +
	theme_bw()




p <- ggplot() +
	geom_hline(yintercept = 0) +
	geom_point(aes(x = total_cases, y = log_score_diff, group = prediction_week, colour = prediction_season, shape = prediction_season),
		size = 4,
		stat="summary", fun.y="median",
		data = pointwise_result_differences) +
	scale_colour_manual("Prediction Season",
		values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")) +
	scale_shape("Prediction Season") +
	xlab("Observed Total Cases in Prediction Target Week") +
	ylab("Median Difference in Log Scores\nAcross all Prediction Horizons for Target Week\n(Positive Values Indicate SSR outperformed SARIMA)") +
	ggtitle("Difference in Log Score vs. Observed Total Cases in Prediction Target Week") +
	theme_bw(base_size = 22)


pdf("/media/evan/data/Reich/dengue-ssr-prediction/inst/intermediate-results/ssr-poster/plots/log-score-diff-by-obs-total-counts.pdf", width=20, height=9.5)
print(p)
dev.off()







