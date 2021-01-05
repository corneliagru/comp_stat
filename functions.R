dens_cond <- function(x, delta_j) {
  return <- vector(mode = "numeric", length = length(x))
  for (i in 1:length(x)) {
    if (x[i] > 0) {
      return[i] <- (delta_j / sqrt(2 * pi)) * exp(delta_j) * (x[i]^(-3 / 2)) * exp((-1 / 2) *
        (((delta_j^2) / x[i]) + x[i]))
    } else {
      return[i] <- 0
    }
  }
  return(return)
}

p_times_dens <- function(x, delta, p) {
  p * sapply(delta, dens_cond, x = x)
}

dens <- function(x, delta, p) {
  if (length(x) == 1) {
    return(p_times_dens(x, delta, p))
  } else {
    return(apply(p_times_dens(x, delta, p), 1, sum))
  }
}

EM_fun <- function(x, k, delta, p, pres = 0.000000001) {
  n <- length(x)
  w <- matrix(NA, n, k)
  result <- matrix(NA, 5000, (2 * k + 2))
  lik_new <- 0
  iter <- 0
  continue <- TRUE
  starting_point <- delta

  while (continue == TRUE) {
    iter <- iter + 1
    lik <- 0

    # E-Step
    temp <- p_times_dens(x, delta, p)
    w <- temp / apply(temp, 1, sum)
    lik <- sum(log(apply(temp, 1, sum)))

    # M-Step
    p <- apply(w, 2, sum) / n
    delta <- apply(w * x, 2, sum) / apply(w, 2, sum)

    # Check Criterion
    criterion <- abs((lik_new - lik) / lik)
    if (criterion < pres) continue <- FALSE

    lik_new <- lik
    result[iter, ] <- c(delta, p, lik, criterion)
  }
  return(list(
    weights = w,
    delta = delta,
    p = p,
    likelihood = lik,
    history = result[1:iter, ],
    starting = starting_point,
    precision = criterion
  ))
}

MCMC_fun <- function(start, delta, p, n, burnin = 1000, discard = 100) {
  nsteps <- burnin + n * discard
  all_values <- matrix(NA, nsteps, length(start))
  x <- start

  for (i in 1:nsteps) {
    y <- runif(1, x - 3, x + 3)
    accept <- min(1, dens(y, delta, p) / dens(x, delta, p))
    u <- runif(1, 0, 1)
    if (u < accept) x <- y
    all_values[i, ] <- x
  }

  keep_values <- all_values[-(1:burnin), ][seq(1, n * discard, 100)]

  return(list(
    all_values = all_values,
    keep_values = keep_values
  ))
}
