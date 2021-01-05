

# import ------------------------------------------------------------------


data <- read.table("data.txt")
x <- data$V1

source("preparation.R")

# find best parameter settings. apparently there is only one group, all deltas
# are basically the same
param <- EM(x, k = 2)
param$delta
param <- EM(x, k = 1)
delta <- param$delta
p <- param$p 


res <- MCMC(theta_0 = 1, delta = delta, p =  p, n = 15000)

acf(res,lag.max=50)
# -> n = 2000, lag 35
plot(ts(res))
# -> n = 2000, burn in 100

samp <- MCMC_iid(res, 40, 100)


plot(density(samp))
lines(density(x), col = "red")


res_bs <- bs_param(B = 10, delta = delta, p = p)

quantile(res_bs[,1], c(0.025, 0.975))
