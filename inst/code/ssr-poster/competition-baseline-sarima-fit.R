library(ssr)
library(forecast)

sj_log_incidence_ts <- ts(San_Juan_train$total_cases, frequency = 52)

sj_log_sarima_fit_manual1 <- arima(sj_log_incidence_ts,
	order = c(1, 0, 0),
	seasonal = list(order = c(4, 1, 0), period = 52))

#sj_sarima_fit <- auto.arima(sj_incidence_ts)
#
#sj_log_incidence_ts <- ts(log(San_Juan_train$total_cases + 1), frequency = 52)
#sj_log_sarima_fit <- auto.arima(sj_log_incidence_ts)
#
#sj_log_sarima_fit_manual1 <- arima(sj_log_incidence_ts,
#	order = c(1, 1, 3),
#	seasonal = list(order = c(1, 0, 0), period = 52))
#
#sj_log_sarima_fit_manual2 <- arima(sj_log_incidence_ts,
#	order = c(1, 1, 3),
#	seasonal = list(order = c(0, 1, 0), period = 52))
#
#sj_log_sarima_fit_manual3 <- arima(sj_log_incidence_ts,
#	order = c(1, 1, 3),
#	seasonal = list(order = c(0, 0, 1), period = 52))
#
#sj_log_incidence <- log(San_Juan_train$total_cases + 1)
#
#sj_seasonally_differenced_log_ts <- ts(sj_log_incidence[seq(from = 53, to = length(sj_log_incidence))] -
#    sj_log_incidence[seq(from = 1, to = length(sj_log_incidence) - 52)], frequency = 52)
#
#temp <- sj_log_incidence[seq(from = 53, to = length(sj_log_incidence))] -
#    sj_log_incidence[seq(from = 1, to = length(sj_log_incidence) - 52)]
#
#sj_twice_seasonally_differenced_log_ts <- ts(temp[seq(from = 1 + 3 * 52, to = length(temp))] -
#	temp[seq(from = 1, to = length(temp) - 3 * 52)])
#
#sj_seasonally_differenced_log_sarima_fit <- auto.arima(sj_seasonally_differenced_log_ts)
#
#coef(sj_seasonally_differenced_log_sarima_fit)
#
## adjust model coefficients and data for effects of seasonal differencing
## this means we can use R's built in forecast.Arima and predict.Arima methods
#sj_seasonally_differenced_log_sarima_fit$coef["sar1"] <- sj_seasonally_differenced_log_sarima_fit$coef["sar1"] + 1
#sj_seasonally_differenced_log_sarima_fit$series <- "sj_log_incidence_ts"
#sj_seasonally_differenced_log_sarima_fit$x <- sj_log_incidence_ts

#plot(forecast(sj_seasonally_differenced_log_sarima_fit,h=30))

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
predictions_df$AE <- NA
predictions_df$predictive_50pct_lb <- NA
predictions_df$predictive_50pct_ub <- NA
predictions_df$predictive_90pct_lb <- NA
predictions_df$predictive_90pct_ub <- NA
predictions_df$week_start_date <- San_Juan_test$week_start_date[1]



#sj_log_sarima_fit <- sj_log_sarima_fit_manual1  ## very bad
#sj_log_sarima_fit <- sj_log_sarima_fit_manual2  ## quite bad
#sj_log_sarima_fit <- sj_log_sarima_fit_manual3  ## very bad

sarima_inds <- which(predictions_df$model == "sarima")
#sarima_inds <- which(predictions_df$model == "sarima" & predictions_df$ph %in% c(1, 13, 26, 39, 52))
last_obs_inds <- sapply(sarima_inds, function(predictions_df_row_ind) {
	which(San_Juan_test$season == predictions_df$last_obs_season[predictions_df_row_ind] &
		San_Juan_test$season_week == predictions_df$last_obs_week[predictions_df_row_ind])
})

for(last_obs_ind in unique(last_obs_inds)) {
	new_data <- ts(San_Juan_test$total_cases[seq_len(last_obs_ind)], frequency = 52)
	updated_sj_log_sarima_fit <- Arima(new_data, model = sj_log_sarima_fit_manual1)

	predict_result <- predict(updated_sj_log_sarima_fit, n.ahead = 52)
	for(ph in seq_len(52)) {
		predictions_df_row_ind <- which(last_obs_inds == last_obs_ind &
			predictions_df$ph == ph)
		predictions_df$week_start_date[predictions_df_row_ind] <- San_Juan_test$week_start_date[last_obs_ind + ph]
	
		predictive_log_mean <- as.numeric(predict_result$pred[ph])
		predictions_df$prediction[predictions_df_row_ind] <- as.numeric(predict_result$pred[ph])

		predictions_df$log_score[predictions_df_row_ind] <- dnorm(San_Juan_test$total_cases[last_obs_ind + ph],
			mean = predictive_log_mean,
			sd = as.numeric(predict_result$se[ph]),
			log = TRUE)
		temp <- qnorm(c(0.05, 0.25, 0.75, 0.95),
			mean = predictive_log_mean,
			sd = as.numeric(predict_result$se[ph]))
		predictions_df[predictions_df_row_ind, c("predictive_90pct_lb", "predictive_50pct_lb", "predictive_50pct_ub", "predictive_90pct_ub")] <-
			temp - 1
	}
}

predictions_df <- predictions_df[predictions_df$week_start_date %in% San_Juan_test$week_start_date[San_Juan_test$season %in% c("2009/2010", "2010/2011", "2011/2012", "2012/2013")], ]

predictions_df$ph <- as.factor(predictions_df$ph)

competition_baseline_sarima_predictions_df <- predictions_df
save(competition_baseline_sarima_predictions_df, file = "/media/evan/data/Reich/dengue-ssr-prediction/inst/intermediate-results/ssr-poster/predictions/competition-baseline-sarima-predictions.Rdata")

library(ggplot2)
ggplot() +
	geom_line(aes(x = week_start_date, y = total_cases), data = San_Juan_test) +
	geom_line(aes(x = week_start_date, y = prediction, colour = ph), data = predictions_df[predictions_df$ph %in% c(1, 13, 26, 39, 52), , drop = FALSE]) +
	theme_bw()


