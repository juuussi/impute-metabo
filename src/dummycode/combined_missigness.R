y <- reference_data

Generate <- function(condition, fixed_objects = NULL) {
  fun <- function(y) ifelse(y > 1, .5, 0)
  cormat <- matrix(.5, 3, 3)
  diag(cormat) <- 1
  dat <- rmvnorm(condition$N, sigma = cormat)
  dat <- apply(dat, 2, add_missing, fun=fun)
  dat
}

results <- runSimulation(Design, replications=1000, verbose=FALSE, parallel = TRUE)
