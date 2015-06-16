# functions to fit ssr and make predictions

ssr_control <- function(lag = 1, dist_fn = dist, dist_fn_args = list(method = "euclidean")) {
	# assemble a list of control parameters for the ssr function
	#  - lag is the number of lags to include
	#  - distance_fn is a distance method 

	control <- list()

	control$lag <- lag

	control$dist_fn <- dist_fn
	control$dist_fn_args <- dist_fn_args

	return(control)
}

ssr <- function(train_data, predict_data, control = ssr_control()) {
	# estimate the parameters for ssr
	#  - each of train_data, predict_vector is a vector (for now we only do univariate time series)
	#  - control is a list of parameters.  see the documentation for ssr_control for a description of these parameters

	# do some estimation process in here and return the results?

	return(list(control = control))
}

ssr_predict <- function(ssr_fit, train_data, predict_data, prediction_lags, k) {
	# ssr_fit is a fit object from an ssr model
	# train_data and predict_data are vectors
	# prediction_lags is a vector of steps ahead to perform prediction.
	# k is the maximum number of non-zero observation weights.  
	#
	# returns a list of lists.  the list has one component for each lag in prediction_lags; component p has two sub-components:
	#  weights: a W_p by J matrix, where W_p is J = nrow(predict_data).  component p of the list is a W by matrix
	#
	# three steps to prediction:
	#  1) compute distances between all points of the form train_data[i:(i + ssr_fit$lag)] and predict_data[j:(j + ssr_fit$lag)]
	#  2) translate distances to weights
	#  3) trace forward to get observations corresponding to each weight.

	# step 1 -- compute distances between lagged observation vectors from train_data and predict_data.
	dists <- compute_pairwise_lagged_obs_distances(train_data, predict_data, ssr_fit$lag, ssr_fit$dist_fn, ssr_fit$dist_fn_args)

	# step 2 -- compute weights
	# TO DO: enforce the k argument -- at most k weights are non-zero.  note that the results of this restriction will be different for each prediction lag

	# I'm following the formulas in Perretti et al. here, but I think there's cancellation -- we should work it out.
	# w_ij = exp(-theta * d_ij / a_j), where a_j = (1/n) sum_i d_ij is the average distance between prediction case j and all training points.
	# in terms of the dists matrix, the a_j are column means
	a <- apply(dists, 2, mean)
	weights <- exp(sweep(dists, 2, -1 * ssr_fit$theta * a, FUN = "/"))


	# step 3 -- get look-ahead observations corresponding to each weight.  again, I'm returning much more than necessary.
	predictions <- lapply(prediction_lags, function(l) {
		# get inds such that we can look ahead the required number l of lags:
		# we require that i + ssr_fit$lag + l <= length(train_data), or equivalently
		# i <= length(train_data) - ssr_fit$lag - l
		inds <- seq_len(length(train_data) - ssr_fit$lag - l)
		return(list(weights = weights[inds, , drop = FALSE],
			obs = train_data[inds + ssr_fit$lag + l]))
	})
	names(predictions) <- paste0("lag_", prediction_lags)

	return(predictions)
}

get_lagged_obs_matrix <- function(orig_data, lag) {
	# orig_data is a vector, lag is an integer >= 0
	n <- length(orig_data)

	lag <- as.integer(lag)
	if(is.na(lag) || lag < 0 || lag > n - 1) {
		stop("invalid lag")
	}

	return(sapply(seq(from = 0, to = lag),
		function(l) {
			orig_data[seq(from = 1 + l, to = n - lag + l)]
		}))
}

compute_pairwise_lagged_obs_distances <- function(v1, v2, lag, dist_fn, dist_fn_args) {
	# this is a first pass -- wastes memory and processing power, but easy to understand and change

	# form matrices with lagged observations
	m1 <- get_lagged_obs_matrix(v1, lag)
	m2 <- get_lagged_obs_matrix(v2, lag)

	# m1 and m2 are matrices with the same number of columns, dist_fn is a function to compute distances, dist_fn_args are arguments to dist_fn
	# row (i, j) of the result has the distance between row i of train and row j of predict.
	return(outer(seq_len(nrow(m1)), seq_len(nrow(m2)),
		FUN = Vectorize(function(i, j) {
			dist_fn_args$x <- rbind(m1[i, ], m2[j, ])
			do.call(dist_fn, dist_fn_args)
		})))
}
