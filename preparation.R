library(doParallel)
library(foreach)

# density function
dens <- function(x, delta) {
  if (x < 0) {
    return(0 * delta)
  } # so dimension fit
  (delta / sqrt(2 * pi)) * exp(delta) * (x^(-1.5)) *
    exp((-0.5) * (((delta^2) / x) + x))
}

#for plotting
dens_plot <- function(x, delta, p) dens(x, delta) %*% t(p)



# EM algorithm to find best delta
EM <- function(x, k = 2, initial = 1, sta = seq(min(x):max(x)),
               tolerance = 0.000000001, silent = TRUE) {
  n <- length(x)
  # initialize delta with quntiles in equal distance between 0.1 and 0.9
  if (initial == 1) {
    delta <- quantile(x, seq(0.1, 0.9, length = k))
  }
  # otherwise initialize with given starting values
  if (initial != 1) delta <- sta
  # initially equal probability for each class
  p <- rep(1 / k, k)
  # initialize w_ij matrix
  w <- matrix(NA, n, k)
  # initialize result matrix for 5000 steps
  # TODO smarter way to inizialize size
  result <- matrix(NA, 5000, (2 * k + 2))
  liknew <- 0
  iter <- 0
  finished <- 0
  startingpoint <- delta
  while (finished == 0) {
    iter <- iter + 1
    lik <- 0
    # E Step
    for (i in 1:n) {
      # todo is hier ergebnis vektor?
      temp <- dens(x[i], delta) * p
      w[i, ] <- temp / sum(temp)
      lik <- lik + log(sum(temp))
    }
    p <- (rep(1, n) %*% w) / n
    delta <- (t(x) %*% w) / (rep(1, n) %*% w)
    criterion <- abs((liknew - lik) / lik)
    if (criterion < tolerance) finished <- 1
    liknew <- lik
    if (!silent) print(c(iter, p, delta, lik))
    # add another 5000 rows if result matrix not enough
    if (iter == nrow(result)) result <- rbind(result, matrix(NA, 5000, (2 * k + 2)))
    result[iter, ] <- c(p, delta, lik, criterion)
  }
  list(
    history = result[1:iter, ], weights = w, starting = startingpoint,
    likelihood = lik, precision = criterion, delta = delta, p = p
  )
}


# non parametric bootstrap for delta conf. intervals
bs_non_param <- function(x, k = 2, initial = 1, sta = seq(min(x):max(x)),
                         tolerance = 0.000000001, seed = 1234, B = 1000, alpha = 0.05, silent = TRUE) {
  set.seed(seed)
  res <- matrix(NA, B, k)
  n <- length(x)
  for (i in 1:B) {
    x_bs <- sample(x, n, replace = TRUE)
    temp <- EM(x_bs,
      k = k, initial = initial, sta = sta,
      tolerance = tolerance, silent = silent
    )
    res[i, ] <- temp$delta
  }
  apply(res, 2, quantile, probs = c(alpha / 2, 1 - alpha / 2))
}


# sample values in a markov chain
MCMC <- function(theta_0, delta, p, n) {
  res <- matrix(NA, n, 1)
  theta <- theta_0

  for (i in 1:n) {
    theta_i <- rnorm(1, theta, 1)
    alpha <- min(1, (p %*% t(dens(theta_i, delta)) / p %*% t(dens(theta, delta))))
    if (rbinom(1, 1, prob = alpha) == 1) theta <- theta_i
    res[i, ] <- theta
  }
  return(res)
}


# get iid samples from markov chain
MCMC_iid <- function(res, n = 300, lag = 35, burn_in = 100) {
  len <- nrow(res)
  keep <- res[burn_in:len, ][seq(1, len - burn_in, by = lag)]
  keep <- keep[1:n]
  keep
}


# bootstrap CI from MCMC_iid
bs_param <- function(B = 100, delta, p, theta_0 = 1, n = 15000, lag = 40, burn_in = 100, parallel = TRUE) {
  k <- length(delta)
  result <- matrix(NA, B, 2 * k)
  if(parallel) {
    foreach (i=1:B, .combine=rbind) %dopar% {
      source("preparation.R")
      res <- MCMC(theta_0 = theta_0, delta = delta, p = p, n = n)
      samp <- MCMC_iid(res, lag = lag, burn_in = burn_in)
      param <- EM(samp, k = k)
      cbind(param$delta, param$p)
    }
  } else {
    for (i in 1:B) {
      res <- MCMC(theta_0 = theta_0, delta = delta, p = p, n = n)
      samp <- MCMC_iid(res, lag = lag, burn_in = burn_in)
      param <- EM(samp, k = k)
      result[i, ] <- cbind(param$delta, param$p)
      result
    }
  }
}




# debugging ---------------------------------------------------------------



B <- 100
theta_0 <- 1
n <- 15000
lag <- 40
burn_in <- 100
