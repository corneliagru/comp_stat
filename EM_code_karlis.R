#### truncated poisson example


#### the loglikelihood function
likel <- function(x, lambda) {
  n <- length(x)
  sumx <- sum(x)
  log(lambda) * sumx - n * log(exp(lambda) - 1) -
    sum(lfactorial(x))
}


##### Newton-Raphson
likfit <- function(x, l1 = 1) {
  i <- 1
  histl <- l1
  t <- 0
  while (t == 0) {
    i <- i + 1
    l1new <- l1 - (l1 - mean(x) + mean(x) * exp(-l1)) / (1 - mean(x) * exp(-l1))
    histl <- c(histl, l1new)
    if (abs(l1 - l1new) < 0.0000001) t <- 10
    l1 <- l1new
  }
  histl
}

### the data
x <- c(rep(1, 797), rep(2, 301), rep(3, 77), rep(4, 17), rep(5, 6), 6, 7)


### function for the EM algorithm
likfit2 <- function(x, l1 = 1) {
  n <- length(x)
  sumx <- sum(x)
  i <- 1
  s0 <- 797 / l1
  histl <- NULL
  t <- 0
  while (t == 0) {
    i <- i + 1
    l1new <- (sumx) / (s0 + n)
    s0 <- (s0 + n) * exp(-l1new)
    newlik <- likel(x, l1new)
    histl <- rbind(histl, c(l1new, newlik))
    if (abs(l1 - l1new) < 0.00000000001) t <- 10
    l1 <- l1new
  }
  histl
}


#############################################
### EM for Poisson mixture
#################################################

data <- read.table("e:\\office pc\\karlis\\compstat\\eglimata.txt")
# x num of crimes
# tr population
x <- data$V1
tr <- data$V2

### initial value
k <- 4 # number of components
sta <- c(1, 10, 20, 30) ## initial values for lambdas
initial <- 0
pres <- 0.000000001 # precisions requeired to stop
p <- rep(1 / k, k) # initial values for mixing proportions

EMpois <- function(x, tr, k = 2, initial = 1, sta = seq(1:max(x), length = k),
                   pres = 0.000000001) {
  n <- length(x)

  if (initial == 1) {
    te <- seq(0.1, 0.9, length = k)
    lambda <- quantile(x / tr, te) + 0.01
    lambda[k] <- max(x / tr)
  }
  if (initial != 1) lambda <- sta
  p <- rep(1 / k, k)
  w <- matrix(NA, n, k)
  result <- matrix(NA, 5000, (2 * k + 2))
  liknew <- 0
  iter <- 0
  telos <- 0
  startingpoint <- lambda
  while (telos == 0) {
    iter <- iter + 1
    lik <- 0
    for (i in 1:n) {
      temp <- dpois(x[i], lambda * tr[i]) * p
      w[i, ] <- temp / sum(temp)
      lik <- lik + log(sum(temp))
    }
    p <- (rep(1, n) %*% w) / n
    lambda <- (t(x) %*% w) / (t(tr) %*% w)
    criterion <- abs((liknew - lik) / lik)
    if (criterion < pres) telos <- 10
    liknew <- lik
    print(c(iter, p, lambda, lik))
    result[iter, ] <- c(p, lambda, lik, criterion)
  }
  list(
    weights = w, lam = lambda, p = p, likelihood = lik, history = result[1:iter, ],
    starting = startingpoint, presicion = criterion
  )
}


# Note that I have supplied this as function as well by removing the
# comments
