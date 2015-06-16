library(ssr)

context("ssr prediction functions")

test_that("compute_pairwise_lagged_obs_distances works", {
	# slow manually computed result that I know is right, using Euclidean distance
	testa <- cbind(1:3, 2:4, 3:5)
	testb <- cbind(2:8, 3:9, 4:10)
	testres <- matrix(NA, nrow = nrow(testa), ncol = nrow(testb))
	for(i in seq_len(nrow(testa))) {
		for(j in seq_len(nrow(testb))) {
			testres[i, j] <- sqrt(sum((testa[i, ] - testb[j, ])^2))
		}
	}

	expect_equal(compute_pairwise_lagged_obs_distances(1:5, 2:10, lag = 2, dist_fn = dist, dist_fn_args = list()),
		testres)
})

#test_that("ssr_predict with euclidean distance works", {
#	ssr_fit <- list(control = list(lag = 2, dist_fn = dist, dist_fn_args = list(method = "euclidean")),
#		theta = 2)
#})
