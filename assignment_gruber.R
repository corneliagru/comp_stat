
# import ------------------------------------------------------------------


data <- read.table("data.txt")
x <- data$V1

source("preparation.R")


# description -------------------------------------------------------------

round(summary(x), digits = 2 )




# find best parameter settings. apparently there is only one group, all deltas
# are basically the same
param <- EM(x, k = 5)
param$delta

# only one class
param <- EM(x, k = 1)
delta <- param$delta
p <- param$p


xx = seq(0,8, by = 0.1)
yy <- sapply(xx, dens_plot, delta = delta, p = p)

plot(density(x), ylim = c(0,0.4), main = "")
lines(xx, yy, col = "red")
legend("topright", legend=c("observed data", "delta = 2.37"),
       col=c("black", "red"),lty = 1,  cex=1.1)






res <- MCMC(theta_0 = 1, delta = delta, p = p, n = 15000)

acf(res, lag.max = 50)
# -> n = 2000, lag 35
plot(ts(res))
# -> n = 2000, burn in 100

samp <- MCMC_iid(res, 40, 100)


plot(density(samp))
lines(density(x), col = "red")


res_bs <- bs_param(B = 10, delta = delta, p = p)

apply(res_bs, 2, quantile, probs = c(0.025, 0.975))



#prepare
param <- EM(x, k = 1)
delta <- param$delta
p <- param$p

#running on server/cluster so I can use all 40 cores
numCores <- detectCores()
registerDoParallel(numCores)
source("preparation.R")
start <- Sys.time()
res_bs <- bs_param(B = 1000, delta = delta, p = p, parallel = TRUE)
end <- Sys.time()
print(end - start)
stopImplicitCluster()





ci <- apply(res_bs, 2, quantile, probs = c(0.025, 0.975))
d_lower <- ci[1,1]
d_upper <- ci[2,1]

xx = seq(0,8, by = 0.1)
yy <- sapply(xx, dens_plot, delta = d_lower, p = p)
yyy <- sapply(xx, dens_plot, delta = d_upper, p = p)


plot(density(x), ylim = c(0,0.4), main = "")
lines(xx, yy, col = "gray20", lty = "dotted")
lines(xx, yyy, col = "gray20", lty = "longdash")
legend("topright", legend=c("delta = 2.19", "delta = 2.59"),
       col=c("gray20", "gray20"), lty=c("dotted", "longdash"), cex=1.1)
